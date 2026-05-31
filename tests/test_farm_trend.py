"""Tests for farm_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.farm_risk import FarmAlertLevel, FarmPanelFlags
from lamp_forge.farm_trend import (
    FarmRecord,
    FarmTrend,
    FarmTrendDirection,
    analyse_farm_trend,
    records_from_csv,
    records_from_farm_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    asfv: bool = False,
    fmdv: bool = False,
    aiv: bool = False,
    ndv: bool = False,
    prrsv: bool = False,
) -> FarmPanelFlags:
    return FarmPanelFlags(asfv=asfv, fmdv=fmdv, aiv=aiv, ndv=ndv, prrsv=prrsv)


def _rec(
    sample_id: str,
    flags: FarmPanelFlags | None = None,
    date: str | None = None,
) -> FarmRecord:
    return FarmRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: FarmRecord) -> FarmTrend:
    return analyse_farm_trend(list(recs))


# ---------------------------------------------------------------------------
# Insufficient data / edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_farm_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(prrsv=True)))
        assert t.direction is FarmTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_notifiable_empty(self) -> None:
        t = _analyse(_rec("S1", _flags(asfv=True)))
        assert t.notifiable_newly_detected == ()
        assert t.notifiable_cleared == ()


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is FarmTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is FarmTrendDirection.STABLE_CLEAR

    def test_stable_clear_interpretation_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        text = t.interpretation.lower()
        assert "negative" in text or "no active" in text

    def test_worst_alert_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is FarmAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# EMERGING
# ---------------------------------------------------------------------------


class TestEmerging:
    def test_notifiable_newly_detected(self) -> None:
        # AIV not in any previous sample, then appears in the last.
        recs = [
            _rec("S1", _flags(prrsv=True)),
            _rec("S2", _flags(prrsv=True)),
            _rec("S3", _flags(prrsv=True, aiv=True)),
        ]
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.EMERGING
        assert "AIV" in t.notifiable_newly_detected

    def test_asfv_newly_detected_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(asfv=True)))
        assert t.direction is FarmTrendDirection.EMERGING
        assert "ASFV" in t.notifiable_newly_detected

    def test_burden_increasing_triggers_emerging(self) -> None:
        # First half has 0, second half has 2 -> delta = 2 >= 0.5
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3", _flags(prrsv=True)),
            _rec("S4", _flags(prrsv=True, ndv=True)),
        ]
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.EMERGING

    def test_emerging_action_mentions_quarantine_or_notify(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(aiv=True)))
        text = t.recommended_action.lower()
        assert "notify" in text or "quarantine" in text or "national" in text

    def test_multiple_notifiable_newly_detected(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(asfv=True, fmdv=True)))
        assert "ASFV" in t.notifiable_newly_detected
        assert "FMDV" in t.notifiable_newly_detected

    def test_cleared_note_appended_when_mixed_events(self) -> None:
        # NDV was present before, AIV is new -- cleared note should appear.
        recs = [
            _rec("S1", _flags(ndv=True)),
            _rec("S2", _flags(aiv=True)),  # NDV gone, AIV new
        ]
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.EMERGING
        assert "AIV" in t.notifiable_newly_detected
        assert "NDV" in t.notifiable_cleared
        assert "NDV" in t.interpretation  # cleared note in interpretation


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_first_sample_positive_last_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(prrsv=True)), _rec("S2"))
        assert t.direction is FarmTrendDirection.RESOLVING

    def test_notifiable_cleared(self) -> None:
        recs = [
            _rec("S1", _flags(aiv=True)),
            _rec("S2", _flags(aiv=True)),
            _rec("S3"),  # AIV cleared
        ]
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.RESOLVING
        assert "AIV" in t.notifiable_cleared

    def test_burden_decreasing_triggers_resolving(self) -> None:
        # early half: 2 pathogens; recent half: 0 pathogens
        recs = [
            _rec("S1", _flags(prrsv=True, ndv=True)),
            _rec("S2", _flags(prrsv=True, ndv=True)),
            _rec("S3"),
            _rec("S4"),
        ]
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.RESOLVING

    def test_resolving_action_mentions_maintain_or_continue(self) -> None:
        t = _analyse(_rec("S1", _flags(prrsv=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "maintain" in text or "continue" in text or "surveillance" in text


# ---------------------------------------------------------------------------
# STABLE_ENDEMIC
# ---------------------------------------------------------------------------


class TestStableEndemic:
    def test_prrsv_consistently_positive(self) -> None:
        recs = [
            _rec("S1", _flags(prrsv=True)),
            _rec("S2", _flags(prrsv=True)),
            _rec("S3", _flags(prrsv=True)),
        ]
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.STABLE_ENDEMIC

    def test_prrsv_endemic_flag_true_when_all_positive(self) -> None:
        recs = [
            _rec("S1", _flags(prrsv=True)),
            _rec("S2", _flags(prrsv=True)),
        ]
        t = analyse_farm_trend(recs)
        assert t.prrsv_endemic is True

    def test_prrsv_endemic_false_when_some_negative(self) -> None:
        recs = [
            _rec("S1", _flags(prrsv=True)),
            _rec("S2"),
            _rec("S3", _flags(prrsv=True)),
        ]
        t = analyse_farm_trend(recs)
        assert t.prrsv_endemic is False

    def test_stable_endemic_interpretation_non_empty(self) -> None:
        recs = [_rec("S1", _flags(prrsv=True)), _rec("S2", _flags(prrsv=True))]
        t = analyse_farm_trend(recs)
        assert t.interpretation.strip()

    def test_endemic_note_in_interpretation_when_prrsv_endemic(self) -> None:
        recs = [
            _rec("S1", _flags(prrsv=True)),
            _rec("S2", _flags(prrsv=True)),
            _rec("S3", _flags(prrsv=True)),
        ]
        t = analyse_farm_trend(recs)
        assert "endemic" in t.interpretation.lower() or "prrsv" in t.interpretation.lower()


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is FarmAlertLevel.NEGATIVE

    def test_asfv_detection_gives_critical_worst(self) -> None:
        recs = [_rec("S1"), _rec("S2", _flags(asfv=True)), _rec("S3")]
        t = analyse_farm_trend(recs)
        assert t.worst_alert_level is FarmAlertLevel.CRITICAL

    def test_prrsv_only_gives_low_worst(self) -> None:
        recs = [_rec("S1", _flags(prrsv=True)), _rec("S2", _flags(prrsv=True))]
        t = analyse_farm_trend(recs)
        assert t.worst_alert_level is FarmAlertLevel.MODERATE

    def test_worst_preserved_after_clearance(self) -> None:
        recs = [
            _rec("S1", _flags(aiv=True)),
            _rec("S2"),  # cleared, but worst should still reflect the earlier high
        ]
        t = analyse_farm_trend(recs)
        assert t.worst_alert_level is FarmAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(prrsv=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "notifiable_newly_detected",
            "notifiable_cleared",
            "prrsv_endemic",
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
        t = _analyse(_rec("S1", _flags(aiv=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_notifiable_newly_detected_is_list(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(asfv=True)))
        d = t.to_dict()
        assert isinstance(d["notifiable_newly_detected"], list)
        assert "ASFV" in d["notifiable_newly_detected"]


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
                    "sample_id": "FarmA-Jan",
                    "date": "2026-01-15",
                    "asfv": "0",
                    "fmdv": "0",
                    "aiv": "1",
                    "ndv": "0",
                    "prrsv": "1",
                },
                {
                    "sample_id": "FarmA-Feb",
                    "date": "2026-02-15",
                    "asfv": "0",
                    "fmdv": "0",
                    "aiv": "0",
                    "ndv": "0",
                    "prrsv": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "FarmA-Jan"
        assert recs[0].flags.aiv is True
        assert recs[0].flags.prrsv is True
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.aiv is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "asfv": "true",
                    "fmdv": "false",
                    "aiv": "TRUE",
                    "ndv": "FALSE",
                    "prrsv": "true",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.asfv is True
        assert recs[0].flags.fmdv is False
        assert recs[0].flags.aiv is True
        assert recs[0].flags.ndv is False
        assert recs[0].flags.prrsv is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0", "prrsv": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0"}],
        )
        with pytest.raises(ValueError, match="prrsv"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0", "prrsv": "1"},
                {"sample_id": "B", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0", "prrsv": "1"},
                {
                    "sample_id": "C",
                    "asfv": "0",
                    "fmdv": "0",
                    "aiv": "1",
                    "ndv": "0",
                    "prrsv": "1",
                },
            ],
        )
        recs = records_from_csv(p)
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.EMERGING
        assert "AIV" in t.notifiable_newly_detected


# ---------------------------------------------------------------------------
# records_from_farm_json
# ---------------------------------------------------------------------------


class TestRecordsFromFarmJson:
    def _write_farm_json(
        self,
        path: Path,
        asfv: bool = False,
        fmdv: bool = False,
        aiv: bool = False,
        ndv: bool = False,
        prrsv: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "asfv": asfv,
                "fmdv": fmdv,
                "aiv": aiv,
                "ndv": ndv,
                "prrsv": prrsv,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_farm_json(p1, prrsv=True)
        self._write_farm_json(p2)
        recs = records_from_farm_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.prrsv is True
        assert recs[1].flags.prrsv is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_farm_json(p, date="2026-03-15")
        recs = records_from_farm_json([p])
        assert recs[0].date == "2026-03-15"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_farm_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, (aiv, prrsv) in enumerate([(False, True), (False, True), (True, True)]):
            p = tmp_path / f"s{i}.json"
            self._write_farm_json(p, aiv=aiv, prrsv=prrsv)
            files.append(p)
        recs = records_from_farm_json(files)
        t = analyse_farm_trend(recs)
        assert t.direction is FarmTrendDirection.EMERGING
        assert "AIV" in t.notifiable_newly_detected


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(prrsv=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(aiv=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"
        assert "AIV" in data["notifiable_cleared"]

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
            "asfv",
            "fmdv",
            "aiv",
            "ndv",
            "prrsv",
            "alert_level",
        ]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(asfv=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_aiv_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(aiv=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["aiv"] == "1"
        assert dr[1]["aiv"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestFarmTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["farm-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        aiv: bool = False,
        prrsv: bool = False,
        asfv: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "asfv": asfv,
                "fmdv": False,
                "aiv": aiv,
                "ndv": False,
                "prrsv": prrsv,
            },
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_farm_result_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, prrsv=True)
        result = self._run(["--farm-result", str(p1), "--farm-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, aiv=True)
        result = self._run(["--farm-result", str(p1), "--farm-result", str(p2)])
        assert "EMERGING" in result.output  # type: ignore[union-attr]

    def test_resolving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, prrsv=True)
        self._write_json(p2)
        result = self._run(["--farm-result", str(p1), "--farm-result", str(p2)])
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2, p3 = tmp_path / "s1.json", tmp_path / "s2.json", tmp_path / "s3.json"
        self._write_json(p1)
        self._write_json(p2)
        self._write_json(p3)
        result = self._run(
            ["--farm-result", str(p1), "--farm-result", str(p2), "--farm-result", str(p3)]
        )
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, prrsv=True)
        result = self._run(["--farm-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0", "prrsv": "0"},
                {"sample_id": "B", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0", "prrsv": "1"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, aiv=True)
        out = tmp_path / "trend.json"
        result = self._run(
            ["--farm-result", str(p1), "--farm-result", str(p2), "--out-json", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "EMERGING"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, prrsv=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            ["--farm-result", str(p1), "--farm-result", str(p2), "--out-csv", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_farm_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "asfv": "0", "fmdv": "0", "aiv": "0", "ndv": "0", "prrsv": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--farm-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "farm_a_jan.json"
        p2 = tmp_path / "farm_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, aiv=True)
        result = self._run(["--farm-result", str(p1), "--farm-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "farm_a_jan" in output
        assert "farm_a_feb" in output

    def test_notifiable_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, asfv=True)
        result = self._run(["--farm-result", str(p1), "--farm-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "ASFV" in output
        assert "NEWLY DETECTED" in output or "newly" in output.lower()
