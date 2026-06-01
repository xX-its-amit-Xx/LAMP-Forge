"""Tests for swine_enteric_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.swine_enteric_risk import SwineEntericAlertLevel, SwineEntericFlags
from lamp_forge.swine_enteric_trend import (
    SwineEntericRecord,
    SwineEntericTrend,
    SwineEntericTrendDirection,
    analyse_swine_enteric_trend,
    records_from_csv,
    records_from_enteric_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    pedv: bool = False,
    pdcov: bool = False,
    rota_a: bool = False,
) -> SwineEntericFlags:
    return SwineEntericFlags(pedv=pedv, pdcov=pdcov, rota_a=rota_a)


def _rec(
    sample_id: str,
    flags: SwineEntericFlags | None = None,
    date: str | None = None,
) -> SwineEntericRecord:
    return SwineEntericRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: SwineEntericRecord) -> SwineEntericTrend:
    return analyse_swine_enteric_trend(list(recs))


# ---------------------------------------------------------------------------
# Insufficient data / edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_swine_enteric_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(pedv=True)))
        assert t.direction is SwineEntericTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_pedv_events_empty(self) -> None:
        t = _analyse(_rec("S1", _flags(pedv=True)))
        assert t.pedv_newly_detected is False
        assert t.pedv_cleared is False


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is SwineEntericTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is SwineEntericTrendDirection.STABLE_CLEAR

    def test_stable_clear_interpretation_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        text = t.interpretation.lower()
        assert "negative" in text or "no active" in text

    def test_worst_alert_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is SwineEntericAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# EMERGING
# ---------------------------------------------------------------------------


class TestEmerging:
    def test_pedv_newly_detected(self) -> None:
        recs = [
            _rec("S1", _flags(pdcov=True)),
            _rec("S2", _flags(pdcov=True)),
            _rec("S3", _flags(pdcov=True, pedv=True)),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.EMERGING
        assert t.pedv_newly_detected is True

    def test_pedv_newly_detected_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(pedv=True)))
        assert t.direction is SwineEntericTrendDirection.EMERGING
        assert t.pedv_newly_detected is True

    def test_burden_increasing_triggers_emerging(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3", _flags(rota_a=True)),
            _rec("S4", _flags(rota_a=True, pdcov=True)),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.EMERGING

    def test_emerging_action_mentions_quarantine_or_notify(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(pedv=True)))
        text = t.recommended_action.lower()
        assert "quarantine" in text or "notify" in text or "lockdown" in text

    def test_pedv_newly_detected_interpretation_mentions_pedv(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(pedv=True)))
        assert "PEDV" in t.interpretation or "pedv" in t.interpretation.lower()


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_first_sample_positive_last_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        assert t.direction is SwineEntericTrendDirection.RESOLVING

    def test_pedv_cleared(self) -> None:
        recs = [
            _rec("S1", _flags(pedv=True)),
            _rec("S2", _flags(pedv=True)),
            _rec("S3"),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.RESOLVING
        assert t.pedv_cleared is True

    def test_burden_decreasing_triggers_resolving(self) -> None:
        recs = [
            _rec("S1", _flags(rota_a=True, pdcov=True)),
            _rec("S2", _flags(rota_a=True, pdcov=True)),
            _rec("S3"),
            _rec("S4"),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.RESOLVING

    def test_resolving_action_mentions_maintain_or_continue(self) -> None:
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "maintain" in text or "continue" in text or "surveillance" in text

    def test_pedv_cleared_action_mentions_consecutive_negative(self) -> None:
        recs = [_rec("S1", _flags(pedv=True)), _rec("S2")]
        t = analyse_swine_enteric_trend(recs)
        text = t.recommended_action.lower()
        assert "consecutive" in text or "negative" in text or "confirm" in text


# ---------------------------------------------------------------------------
# STABLE_ENDEMIC
# ---------------------------------------------------------------------------


class TestStableEndemic:
    def test_rota_a_consistently_positive(self) -> None:
        recs = [
            _rec("S1", _flags(rota_a=True)),
            _rec("S2", _flags(rota_a=True)),
            _rec("S3", _flags(rota_a=True)),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.STABLE_ENDEMIC

    def test_pdcov_consistently_positive(self) -> None:
        recs = [
            _rec("S1", _flags(pdcov=True)),
            _rec("S2", _flags(pdcov=True)),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.STABLE_ENDEMIC

    def test_stable_endemic_interpretation_non_empty(self) -> None:
        recs = [_rec("S1", _flags(rota_a=True)), _rec("S2", _flags(rota_a=True))]
        t = analyse_swine_enteric_trend(recs)
        assert t.interpretation.strip()

    def test_stable_endemic_with_high_alert_mentions_coronavirus(self) -> None:
        recs = [
            _rec("S1", _flags(pdcov=True)),
            _rec("S2", _flags(pdcov=True)),
            _rec("S3", _flags(pdcov=True)),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert "coronavirus" in t.interpretation.lower() or "pdcov" in t.interpretation.lower()


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is SwineEntericAlertLevel.NEGATIVE

    def test_pedv_detection_gives_critical_worst(self) -> None:
        recs = [_rec("S1"), _rec("S2", _flags(pedv=True)), _rec("S3")]
        t = analyse_swine_enteric_trend(recs)
        assert t.worst_alert_level is SwineEntericAlertLevel.CRITICAL

    def test_rota_a_only_gives_moderate_worst(self) -> None:
        recs = [_rec("S1", _flags(rota_a=True)), _rec("S2", _flags(rota_a=True))]
        t = analyse_swine_enteric_trend(recs)
        assert t.worst_alert_level is SwineEntericAlertLevel.MODERATE

    def test_worst_preserved_after_clearance(self) -> None:
        recs = [
            _rec("S1", _flags(pedv=True)),
            _rec("S2"),
        ]
        t = analyse_swine_enteric_trend(recs)
        assert t.worst_alert_level is SwineEntericAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(rota_a=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "pedv_newly_detected",
            "pedv_cleared",
            "worst_alert_level",
            "interpretation",
            "recommended_action",
            "timeline",
        ):
            assert key in d, f"Missing key: {key}"

    def test_timeline_length_matches_n_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        d = t.to_dict()
        tl = d["timeline"]
        assert isinstance(tl, list)
        assert len(tl) == 3

    def test_timeline_contains_sample_ids(self) -> None:
        t = _analyse(_rec("Alpha"), _rec("Beta"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert tl[0]["sample_id"] == "Alpha"
        assert tl[1]["sample_id"] == "Beta"

    def test_json_serialisable(self) -> None:
        t = _analyse(_rec("S1", _flags(pedv=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_pedv_newly_detected_is_bool(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(pedv=True)))
        d = t.to_dict()
        assert isinstance(d["pedv_newly_detected"], bool)
        assert d["pedv_newly_detected"] is True


# ---------------------------------------------------------------------------
# records_from_csv
# ---------------------------------------------------------------------------


class TestRecordsFromCsv:
    def _write_csv(self, tmp_path: Path, rows: list[dict[str, str]]) -> Path:
        p = tmp_path / "test.csv"
        fieldnames = list(rows[0].keys())
        with p.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        return p

    def test_basic_parse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "BarnA-Jan",
                    "date": "2026-01-15",
                    "pedv": "1",
                    "pdcov": "0",
                    "rota_a": "1",
                },
                {
                    "sample_id": "BarnA-Feb",
                    "date": "2026-02-15",
                    "pedv": "0",
                    "pdcov": "0",
                    "rota_a": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "BarnA-Jan"
        assert recs[0].flags.pedv is True
        assert recs[0].flags.rota_a is True
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.pedv is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "pedv": "true",
                    "pdcov": "false",
                    "rota_a": "TRUE",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.pedv is True
        assert recs[0].flags.pdcov is False
        assert recs[0].flags.rota_a is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "pedv": "0", "pdcov": "0", "rota_a": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "pedv": "0", "pdcov": "0"}],
        )
        with pytest.raises(ValueError, match="rota_a"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "pedv": "0", "pdcov": "0", "rota_a": "1"},
                {"sample_id": "B", "pedv": "0", "pdcov": "0", "rota_a": "1"},
                {"sample_id": "C", "pedv": "1", "pdcov": "0", "rota_a": "1"},
            ],
        )
        recs = records_from_csv(p)
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.EMERGING
        assert t.pedv_newly_detected is True


# ---------------------------------------------------------------------------
# records_from_enteric_json
# ---------------------------------------------------------------------------


class TestRecordsFromEntericJson:
    def _write_enteric_json(
        self,
        path: Path,
        pedv: bool = False,
        pdcov: bool = False,
        rota_a: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "pedv": pedv,
                "pdcov": pdcov,
                "rota_a": rota_a,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_enteric_json(p1, rota_a=True)
        self._write_enteric_json(p2)
        recs = records_from_enteric_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.rota_a is True
        assert recs[1].flags.rota_a is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_enteric_json(p, date="2026-03-15")
        recs = records_from_enteric_json([p])
        assert recs[0].date == "2026-03-15"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_enteric_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, (rota_a, pedv) in enumerate([(True, False), (True, False), (True, True)]):
            p = tmp_path / f"s{i}.json"
            self._write_enteric_json(p, rota_a=rota_a, pedv=pedv)
            files.append(p)
        recs = records_from_enteric_json(files)
        t = analyse_swine_enteric_trend(recs)
        assert t.direction is SwineEntericTrendDirection.EMERGING
        assert t.pedv_newly_detected is True


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(rota_a=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(pedv=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"
        assert data["pedv_cleared"] is True

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "trend.json"
        write_trend_json(_analyse(_rec("S1"), _rec("S2")), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# write_trend_csv
# ---------------------------------------------------------------------------


class TestWriteTrendCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        assert out.exists()

    def test_header_row(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == [
            "sample_id",
            "date",
            "pedv",
            "pdcov",
            "rota_a",
            "alert_level",
        ]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(pedv=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_pedv_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(pedv=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["pedv"] == "1"
        assert dr[1]["pedv"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestSwineEntericTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["swine-enteric-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        pedv: bool = False,
        pdcov: bool = False,
        rota_a: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "pedv": pedv,
                "pdcov": pdcov,
                "rota_a": rota_a,
            },
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_enteric_result_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, rota_a=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, pedv=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        assert "EMERGING" in result.output  # type: ignore[union-attr]

    def test_resolving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, rota_a=True)
        self._write_json(p2)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2, p3 = tmp_path / "s1.json", tmp_path / "s2.json", tmp_path / "s3.json"
        self._write_json(p1)
        self._write_json(p2)
        self._write_json(p3)
        result = self._run(
            [
                "--enteric-result",
                str(p1),
                "--enteric-result",
                str(p2),
                "--enteric-result",
                str(p3),
            ]
        )
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, rota_a=True)
        result = self._run(["--enteric-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "pedv": "0", "pdcov": "0", "rota_a": "0"},
                {"sample_id": "B", "pedv": "0", "pdcov": "0", "rota_a": "1"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, pedv=True)
        out = tmp_path / "trend.json"
        result = self._run(
            ["--enteric-result", str(p1), "--enteric-result", str(p2), "--out-json", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "EMERGING"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, rota_a=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            ["--enteric-result", str(p1), "--enteric-result", str(p2), "--out-csv", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_enteric_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "pedv": "0", "pdcov": "0", "rota_a": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--enteric-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "barn_a_jan.json"
        p2 = tmp_path / "barn_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, pedv=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "barn_a_jan" in output
        assert "barn_a_feb" in output

    def test_pedv_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, pedv=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "PEDV" in output
        assert "NEWLY DETECTED" in output or "newly" in output.lower()
