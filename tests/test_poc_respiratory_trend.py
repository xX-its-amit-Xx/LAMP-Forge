"""Tests for poc_respiratory_trend -- no external binaries required.

Direction logic is deterministic so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.poc_respiratory_risk import PocRespAlertLevel, PocRespFlags
from lamp_forge.poc_respiratory_trend import (
    PocRespRecord,
    PocRespTrend,
    PocRespTrendDirection,
    analyse_poc_resp_trend,
    records_from_csv,
    records_from_resp_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    hrsv: bool = False,
    iav: bool = False,
    sars2: bool = False,
) -> PocRespFlags:
    return PocRespFlags(hrsv=hrsv, iav=iav, sars2=sars2)


def _rec(
    sample_id: str,
    flags: PocRespFlags | None = None,
    date: str | None = None,
) -> PocRespRecord:
    return PocRespRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: PocRespRecord) -> PocRespTrend:
    return analyse_poc_resp_trend(list(recs))


# ---------------------------------------------------------------------------
# Edge cases / insufficient data
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_poc_resp_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)))
        assert t.direction is PocRespTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_no_event_flags(self) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)))
        assert t.iav_newly_detected is False
        assert t.dual_wave_active is False


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is PocRespTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is PocRespTrendDirection.STABLE_CLEAR

    def test_stable_clear_interpretation_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        assert "negative" in t.interpretation.lower()

    def test_worst_alert_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is PocRespAlertLevel.NEGATIVE

    def test_event_flags_false_when_all_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.iav_newly_detected is False
        assert t.dual_wave_active is False


# ---------------------------------------------------------------------------
# EMERGING
# ---------------------------------------------------------------------------


class TestEmerging:
    def test_iav_newly_detected_triggers_emerging(self) -> None:
        recs = [
            _rec("S1", _flags(sars2=True)),
            _rec("S2", _flags(sars2=True)),
            _rec("S3", _flags(sars2=True, iav=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.direction is PocRespTrendDirection.EMERGING
        assert t.iav_newly_detected is True

    def test_iav_newly_detected_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(iav=True)))
        assert t.direction is PocRespTrendDirection.EMERGING
        assert t.iav_newly_detected is True

    def test_iav_newly_detected_action_mentions_season(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(iav=True)))
        text = (t.interpretation + " " + t.recommended_action).lower()
        assert "influenza" in text or "iav" in text or "season" in text

    def test_burden_increasing_triggers_emerging(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3", _flags(iav=True)),
            _rec("S4", _flags(iav=True, sars2=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.direction is PocRespTrendDirection.EMERGING

    def test_emerging_action_mentions_protocol_or_outbreak(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(iav=True, sars2=True, hrsv=True)))
        text = t.recommended_action.lower()
        assert "protocol" in text or "outbreak" in text or "antiviral" in text


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_last_sample_negative_after_positives(self) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        assert t.direction is PocRespTrendDirection.RESOLVING

    def test_burden_decreasing_triggers_resolving(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True, sars2=True)),
            _rec("S2", _flags(iav=True, sars2=True)),
            _rec("S3"),
            _rec("S4"),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.direction is PocRespTrendDirection.RESOLVING

    def test_resolving_action_mentions_maintain_or_surveillance(self) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "maintain" in text or "continue" in text or "surveillance" in text


# ---------------------------------------------------------------------------
# STABLE_ENDEMIC
# ---------------------------------------------------------------------------


class TestStableEndemic:
    def test_consistently_iav_positive(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True)),
            _rec("S2", _flags(iav=True)),
            _rec("S3", _flags(iav=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.direction is PocRespTrendDirection.STABLE_ENDEMIC

    def test_stable_endemic_interpretation_non_empty(self) -> None:
        recs = [_rec("S1", _flags(hrsv=True)), _rec("S2", _flags(hrsv=True))]
        t = analyse_poc_resp_trend(recs)
        assert t.interpretation.strip()

    def test_iav_already_detected_not_flagged_as_newly_detected(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True)),
            _rec("S2", _flags(iav=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.iav_newly_detected is False


# ---------------------------------------------------------------------------
# dual_wave_active flag
# ---------------------------------------------------------------------------


class TestDualWaveActive:
    def test_dual_wave_when_both_in_recent_majority(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True, sars2=True)),
            _rec("S2", _flags(iav=True, sars2=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.dual_wave_active is True

    def test_dual_wave_false_when_not_co_dominant(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True)),
            _rec("S2", _flags(sars2=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.dual_wave_active is False

    def test_dual_wave_note_in_stable_endemic_interpretation(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True, sars2=True)),
            _rec("S2", _flags(iav=True, sars2=True)),
            _rec("S3", _flags(iav=True, sars2=True)),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.dual_wave_active is True
        text = (t.interpretation + " " + t.recommended_action).lower()
        assert "iav" in text or "sars" in text or "dual" in text or "twindemic" in text


# ---------------------------------------------------------------------------
# worst_alert_level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is PocRespAlertLevel.NEGATIVE

    def test_dual_infection_gives_high_worst(self) -> None:
        recs = [_rec("S1"), _rec("S2", _flags(iav=True, sars2=True)), _rec("S3")]
        t = analyse_poc_resp_trend(recs)
        assert t.worst_alert_level is PocRespAlertLevel.HIGH

    def test_hrsv_only_gives_low_worst(self) -> None:
        recs = [_rec("S1", _flags(hrsv=True)), _rec("S2", _flags(hrsv=True))]
        t = analyse_poc_resp_trend(recs)
        assert t.worst_alert_level is PocRespAlertLevel.LOW

    def test_worst_preserved_after_clearance(self) -> None:
        recs = [
            _rec("S1", _flags(iav=True)),
            _rec("S2"),
        ]
        t = analyse_poc_resp_trend(recs)
        assert t.worst_alert_level is PocRespAlertLevel.MODERATE


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(iav=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "iav_newly_detected",
            "dual_wave_active",
            "worst_alert_level",
            "interpretation",
            "recommended_action",
            "timeline",
        ):
            assert key in d, f"Missing key: {key}"

    def test_timeline_length_matches_n_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        d = t.to_dict()
        assert len(d["timeline"]) == 3

    def test_timeline_contains_sample_ids(self) -> None:
        t = _analyse(_rec("Alpha"), _rec("Beta"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert tl[0]["sample_id"] == "Alpha"
        assert tl[1]["sample_id"] == "Beta"

    def test_json_serialisable(self) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_alert_level_in_timeline(self) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        tl = t.to_dict()["timeline"]
        valid = {"NEGATIVE", "LOW", "MODERATE", "HIGH"}
        assert tl[0]["alert_level"] in valid


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
                    "sample_id": "Clinic-Jan",
                    "date": "2026-01-15",
                    "hrsv": "0",
                    "iav": "1",
                    "sars2": "0",
                },
                {
                    "sample_id": "Clinic-Feb",
                    "date": "2026-02-15",
                    "hrsv": "0",
                    "iav": "0",
                    "sars2": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "Clinic-Jan"
        assert recs[0].flags.iav is True
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.iav is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "hrsv": "true", "iav": "FALSE", "sars2": "1"}],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.hrsv is True
        assert recs[0].flags.iav is False
        assert recs[0].flags.sars2 is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "hrsv": "0", "iav": "0", "sars2": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "hrsv": "0", "iav": "0"}],
        )
        with pytest.raises(ValueError, match="sars2"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "hrsv": "0", "iav": "1", "sars2": "0"},
                {"sample_id": "B", "hrsv": "0", "iav": "1", "sars2": "0"},
                {"sample_id": "C", "hrsv": "0", "iav": "1", "sars2": "0"},
            ],
        )
        recs = records_from_csv(p)
        t = analyse_poc_resp_trend(recs)
        assert t.direction is PocRespTrendDirection.STABLE_ENDEMIC
        assert t.iav_newly_detected is False


# ---------------------------------------------------------------------------
# records_from_resp_json
# ---------------------------------------------------------------------------


class TestRecordsFromRespJson:
    def _write_resp_json(
        self,
        path: Path,
        hrsv: bool = False,
        iav: bool = False,
        sars2: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "NEGATIVE",
            "panel_flags": {"hrsv": hrsv, "iav": iav, "sars2": sars2},
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_resp_json(p1, iav=True)
        self._write_resp_json(p2)
        recs = records_from_resp_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.iav is True
        assert recs[1].flags.iav is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_resp_json(p, date="2026-03-15")
        recs = records_from_resp_json([p])
        assert recs[0].date == "2026-03-15"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_resp_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, (iav, sars2) in enumerate([(False, False), (False, False), (True, False)]):
            p = tmp_path / f"s{i}.json"
            self._write_resp_json(p, iav=iav, sars2=sars2)
            files.append(p)
        recs = records_from_resp_json(files)
        t = analyse_poc_resp_trend(recs)
        assert t.direction is PocRespTrendDirection.EMERGING
        assert t.iav_newly_detected is True


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(iav=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"

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
        assert rows[0] == ["sample_id", "date", "hrsv", "iav", "sars2", "alert_level"]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH"}
        assert dr[0]["alert_level"] in valid_levels

    def test_iav_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(iav=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["iav"] == "1"
        assert dr[1]["iav"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestPocRespTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["poc-respiratory-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        hrsv: bool = False,
        iav: bool = False,
        sars2: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "NEGATIVE",
            "panel_flags": {"hrsv": hrsv, "iav": iav, "sars2": sars2},
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_resp_result_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, iav=True)
        result = self._run(["--resp-result", str(p1), "--resp-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, iav=True)
        result = self._run(["--resp-result", str(p1), "--resp-result", str(p2)])
        assert "EMERGING" in result.output  # type: ignore[union-attr]

    def test_iav_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, iav=True)
        result = self._run(["--resp-result", str(p1), "--resp-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "IAV" in output or "season" in output.lower() or "influenza" in output.lower()

    def test_resolving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, iav=True)
        self._write_json(p2)
        result = self._run(["--resp-result", str(p1), "--resp-result", str(p2)])
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2, p3 = tmp_path / "s1.json", tmp_path / "s2.json", tmp_path / "s3.json"
        self._write_json(p1)
        self._write_json(p2)
        self._write_json(p3)
        result = self._run(
            [
                "--resp-result",
                str(p1),
                "--resp-result",
                str(p2),
                "--resp-result",
                str(p3),
            ]
        )
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, iav=True)
        result = self._run(["--resp-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "hrsv": "0", "iav": "0", "sars2": "0"},
                {"sample_id": "B", "hrsv": "0", "iav": "1", "sars2": "1"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, iav=True)
        out = tmp_path / "trend.json"
        result = self._run(
            ["--resp-result", str(p1), "--resp-result", str(p2), "--out-json", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "EMERGING"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, iav=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            ["--resp-result", str(p1), "--resp-result", str(p2), "--out-csv", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_resp_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "hrsv": "0", "iav": "0", "sars2": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--resp-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "clinic_a_jan.json"
        p2 = tmp_path / "clinic_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, iav=True)
        result = self._run(["--resp-result", str(p1), "--resp-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "clinic_a_jan" in output
        assert "clinic_a_feb" in output
