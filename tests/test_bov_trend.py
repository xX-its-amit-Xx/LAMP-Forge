"""Tests for bov_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.bov_risk import BovAlertLevel, BovPanelFlags
from lamp_forge.bov_trend import (
    BovRecord,
    BovTrend,
    BovTrendDirection,
    analyse_bov_trend,
    records_from_bov_json,
    records_from_csv,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    brsv: bool = False,
    bcov: bool = False,
    bvdv: bool = False,
    ibr: bool = False,
    mhae: bool = False,
) -> BovPanelFlags:
    return BovPanelFlags(brsv=brsv, bcov=bcov, bvdv=bvdv, ibr=ibr, mhae=mhae)


def _rec(
    sample_id: str,
    flags: BovPanelFlags | None = None,
    date: str | None = None,
) -> BovRecord:
    return BovRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: BovRecord) -> BovTrend:
    return analyse_bov_trend(list(recs))


# ---------------------------------------------------------------------------
# Insufficient data / edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_bov_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(brsv=True)))
        assert t.direction is BovTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_ibr_flags_empty(self) -> None:
        t = _analyse(_rec("S1", _flags(ibr=True)))
        assert t.ibr_newly_detected is False
        assert t.ibr_cleared is False


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is BovTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is BovTrendDirection.STABLE_CLEAR

    def test_stable_clear_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        text = t.interpretation.lower()
        assert "negative" in text or "no active" in text

    def test_worst_alert_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is BovAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# EMERGING
# ---------------------------------------------------------------------------


class TestEmerging:
    def test_ibr_newly_detected(self) -> None:
        recs = [
            _rec("S1", _flags(brsv=True)),
            _rec("S2", _flags(brsv=True)),
            _rec("S3", _flags(brsv=True, ibr=True)),
        ]
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.EMERGING
        assert t.ibr_newly_detected is True

    def test_ibr_newly_detected_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ibr=True)))
        assert t.direction is BovTrendDirection.EMERGING
        assert t.ibr_newly_detected is True

    def test_burden_increasing_triggers_emerging(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3", _flags(brsv=True)),
            _rec("S4", _flags(brsv=True, bcov=True)),
        ]
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.EMERGING

    def test_emerging_action_mentions_veterinarian_or_monitoring(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ibr=True)))
        text = t.recommended_action.lower()
        assert "veterinarian" in text or "monitoring" in text or "notify" in text

    def test_ibr_emerging_action_mentions_notification(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ibr=True)))
        text = t.recommended_action.lower()
        assert "notify" in text or "notification" in text or "regulatory" in text

    def test_cleared_note_appended_when_ibr_was_positive_before(self) -> None:
        recs = [
            _rec("S1", _flags(ibr=True)),
            _rec("S2", _flags(brsv=True)),  # IBR gone, burden increased
        ]
        t = analyse_bov_trend(recs)
        assert t.ibr_cleared is True


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_first_sample_positive_last_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(brsv=True)), _rec("S2"))
        assert t.direction is BovTrendDirection.RESOLVING

    def test_ibr_cleared(self) -> None:
        recs = [
            _rec("S1", _flags(ibr=True)),
            _rec("S2", _flags(ibr=True)),
            _rec("S3"),
        ]
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.RESOLVING
        assert t.ibr_cleared is True

    def test_ibr_cleared_action_mentions_serology(self) -> None:
        recs = [_rec("S1", _flags(ibr=True)), _rec("S2")]
        t = analyse_bov_trend(recs)
        text = t.recommended_action.lower()
        assert "serology" in text or "confirm" in text or "clearance" in text

    def test_burden_decreasing_triggers_resolving(self) -> None:
        recs = [
            _rec("S1", _flags(brsv=True, mhae=True)),
            _rec("S2", _flags(brsv=True, mhae=True)),
            _rec("S3"),
            _rec("S4"),
        ]
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.RESOLVING

    def test_resolving_action_mentions_maintain_or_surveillance(self) -> None:
        t = _analyse(_rec("S1", _flags(brsv=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "maintain" in text or "surveillance" in text or "continue" in text


# ---------------------------------------------------------------------------
# STABLE_ENDEMIC
# ---------------------------------------------------------------------------


class TestStableEndemic:
    def test_brsv_consistently_positive(self) -> None:
        recs = [
            _rec("S1", _flags(brsv=True)),
            _rec("S2", _flags(brsv=True)),
            _rec("S3", _flags(brsv=True)),
        ]
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.STABLE_ENDEMIC

    def test_bvdv_endemic_flag_true_when_all_positive(self) -> None:
        recs = [
            _rec("S1", _flags(bvdv=True)),
            _rec("S2", _flags(bvdv=True)),
        ]
        t = analyse_bov_trend(recs)
        assert t.bvdv_endemic is True

    def test_bvdv_endemic_false_when_some_negative(self) -> None:
        recs = [
            _rec("S1", _flags(bvdv=True)),
            _rec("S2"),
            _rec("S3", _flags(bvdv=True)),
        ]
        t = analyse_bov_trend(recs)
        assert t.bvdv_endemic is False

    def test_bvdv_endemic_note_in_interpretation(self) -> None:
        recs = [_rec("S1", _flags(bvdv=True)), _rec("S2", _flags(bvdv=True))]
        t = analyse_bov_trend(recs)
        assert "bvdv" in t.interpretation.lower() or "pi" in t.interpretation.lower()

    def test_stable_endemic_interpretation_non_empty(self) -> None:
        recs = [_rec("S1", _flags(brsv=True)), _rec("S2", _flags(brsv=True))]
        t = analyse_bov_trend(recs)
        assert t.interpretation.strip()


# ---------------------------------------------------------------------------
# Bacterial co-infection tracking
# ---------------------------------------------------------------------------


class TestBacterialCoinfection:
    def test_count_zero_when_no_coinfection(self) -> None:
        recs = [_rec("S1", _flags(brsv=True)), _rec("S2", _flags(mhae=True))]
        t = analyse_bov_trend(recs)
        assert t.bacterial_coinfection_count == 0

    def test_count_one_for_single_coinfection_interval(self) -> None:
        recs = [
            _rec("S1", _flags(brsv=True, mhae=True)),
            _rec("S2"),
        ]
        t = analyse_bov_trend(recs)
        assert t.bacterial_coinfection_count == 1

    def test_count_two_for_two_intervals(self) -> None:
        recs = [
            _rec("S1", _flags(brsv=True, mhae=True)),
            _rec("S2", _flags(bcov=True, mhae=True)),
            _rec("S3"),
        ]
        t = analyse_bov_trend(recs)
        assert t.bacterial_coinfection_count == 2

    def test_mhae_alone_not_counted(self) -> None:
        recs = [_rec("S1", _flags(mhae=True)), _rec("S2", _flags(mhae=True))]
        t = analyse_bov_trend(recs)
        assert t.bacterial_coinfection_count == 0


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is BovAlertLevel.NEGATIVE

    def test_brsv_mhae_gives_critical_worst(self) -> None:
        recs = [_rec("S1"), _rec("S2", _flags(brsv=True, mhae=True)), _rec("S3")]
        t = analyse_bov_trend(recs)
        assert t.worst_alert_level is BovAlertLevel.CRITICAL

    def test_worst_preserved_after_clearance(self) -> None:
        recs = [
            _rec("S1", _flags(brsv=True, mhae=True)),
            _rec("S2"),
        ]
        t = analyse_bov_trend(recs)
        assert t.worst_alert_level is BovAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(brsv=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "ibr_newly_detected",
            "ibr_cleared",
            "bvdv_endemic",
            "bacterial_coinfection_count",
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
        t = _analyse(_rec("S1", _flags(brsv=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_ibr_newly_detected_is_bool(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ibr=True)))
        d = t.to_dict()
        assert isinstance(d["ibr_newly_detected"], bool)
        assert d["ibr_newly_detected"] is True


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
                    "sample_id": "HerdA-Jan",
                    "date": "2026-01-15",
                    "brsv": "1",
                    "bcov": "0",
                    "bvdv": "0",
                    "ibr": "0",
                    "mhae": "1",
                },
                {
                    "sample_id": "HerdA-Feb",
                    "date": "2026-02-15",
                    "brsv": "0",
                    "bcov": "0",
                    "bvdv": "0",
                    "ibr": "0",
                    "mhae": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "HerdA-Jan"
        assert recs[0].flags.brsv is True
        assert recs[0].flags.mhae is True
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.brsv is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "brsv": "true",
                    "bcov": "false",
                    "bvdv": "TRUE",
                    "ibr": "FALSE",
                    "mhae": "true",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.brsv is True
        assert recs[0].flags.bcov is False
        assert recs[0].flags.bvdv is True
        assert recs[0].flags.ibr is False
        assert recs[0].flags.mhae is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "brsv": "0", "bcov": "0", "bvdv": "0", "ibr": "0", "mhae": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "brsv": "0", "bcov": "0", "bvdv": "0", "ibr": "0"}],
        )
        with pytest.raises(ValueError, match="mhae"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "brsv": "1", "bcov": "0", "bvdv": "0", "ibr": "0", "mhae": "0"},
                {"sample_id": "B", "brsv": "1", "bcov": "0", "bvdv": "0", "ibr": "0", "mhae": "0"},
                {
                    "sample_id": "C",
                    "brsv": "1",
                    "bcov": "0",
                    "bvdv": "0",
                    "ibr": "1",
                    "mhae": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.EMERGING
        assert t.ibr_newly_detected is True


# ---------------------------------------------------------------------------
# records_from_bov_json
# ---------------------------------------------------------------------------


class TestRecordsFromBovJson:
    def _write_bov_json(
        self,
        path: Path,
        brsv: bool = False,
        bcov: bool = False,
        bvdv: bool = False,
        ibr: bool = False,
        mhae: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "brsv": brsv,
                "bcov": bcov,
                "bvdv": bvdv,
                "ibr": ibr,
                "mhae": mhae,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_bov_json(p1, brsv=True)
        self._write_bov_json(p2)
        recs = records_from_bov_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.brsv is True
        assert recs[1].flags.brsv is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_bov_json(p, date="2026-03-15")
        recs = records_from_bov_json([p])
        assert recs[0].date == "2026-03-15"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_bov_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, (brsv, ibr) in enumerate([(True, False), (True, False), (True, True)]):
            p = tmp_path / f"s{i}.json"
            self._write_bov_json(p, brsv=brsv, ibr=ibr)
            files.append(p)
        recs = records_from_bov_json(files)
        t = analyse_bov_trend(recs)
        assert t.direction is BovTrendDirection.EMERGING
        assert t.ibr_newly_detected is True


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(brsv=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(brsv=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"

    def test_ibr_newly_detected_in_json(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ibr=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["ibr_newly_detected"] is True

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
            "brsv",
            "bcov",
            "bvdv",
            "ibr",
            "mhae",
            "alert_level",
        ]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(brsv=True, mhae=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_brsv_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(brsv=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["brsv"] == "1"
        assert dr[1]["brsv"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestBovTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["bov-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        brsv: bool = False,
        ibr: bool = False,
        mhae: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "brsv": brsv,
                "bcov": False,
                "bvdv": False,
                "ibr": ibr,
                "mhae": mhae,
            },
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_bov_result_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, brsv=True)
        result = self._run(["--bov-result", str(p1), "--bov-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ibr=True)
        result = self._run(["--bov-result", str(p1), "--bov-result", str(p2)])
        assert "EMERGING" in result.output  # type: ignore[union-attr]

    def test_resolving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, brsv=True)
        self._write_json(p2)
        result = self._run(["--bov-result", str(p1), "--bov-result", str(p2)])
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2, p3 = tmp_path / "s1.json", tmp_path / "s2.json", tmp_path / "s3.json"
        self._write_json(p1)
        self._write_json(p2)
        self._write_json(p3)
        result = self._run(
            ["--bov-result", str(p1), "--bov-result", str(p2), "--bov-result", str(p3)]
        )
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, brsv=True)
        result = self._run(["--bov-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "brsv": "0", "bcov": "0", "bvdv": "0", "ibr": "0", "mhae": "0"},
                {"sample_id": "B", "brsv": "1", "bcov": "0", "bvdv": "0", "ibr": "0", "mhae": "0"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ibr=True)
        out = tmp_path / "trend.json"
        result = self._run(
            ["--bov-result", str(p1), "--bov-result", str(p2), "--out-json", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "EMERGING"
        assert data["ibr_newly_detected"] is True

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, brsv=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            ["--bov-result", str(p1), "--bov-result", str(p2), "--out-csv", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_bov_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "brsv": "0", "bcov": "0", "bvdv": "0", "ibr": "0", "mhae": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--bov-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "herd_a_jan.json"
        p2 = tmp_path / "herd_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, ibr=True)
        result = self._run(["--bov-result", str(p1), "--bov-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "herd_a_jan" in output
        assert "herd_a_feb" in output

    def test_ibr_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ibr=True)
        result = self._run(["--bov-result", str(p1), "--bov-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "IBR" in output
        assert "NEWLY DETECTED" in output or "newly" in output.lower()
