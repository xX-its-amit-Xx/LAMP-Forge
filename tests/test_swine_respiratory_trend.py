"""Tests for the swine respiratory / PRDC trend analysis module.

All tests are pure Python (no external binaries required).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import ClassVar

import pytest

from lamp_forge.swine_respiratory_risk import PrdcFlags
from lamp_forge.swine_respiratory_trend import (
    SwineRespRecord,
    SwineRespTrend,
    SwineRespTrendDirection,
    analyse_swine_resp_trend,
    records_from_csv,
    records_from_resp_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _rec(sample_id: str, **kwargs: bool) -> SwineRespRecord:
    defaults: dict[str, bool] = {"prrsv": False, "pcv2": False, "mhp": False}
    defaults.update(kwargs)
    return SwineRespRecord(sample_id=sample_id, flags=PrdcFlags(**defaults))


def _trend(*records: SwineRespRecord) -> SwineRespTrend:
    return analyse_swine_resp_trend(list(records))


# ---------------------------------------------------------------------------
# Empty input raises ValueError
# ---------------------------------------------------------------------------


class TestEmptyInput:
    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_swine_resp_trend([])


# ---------------------------------------------------------------------------
# Insufficient data (single sample)
# ---------------------------------------------------------------------------


class TestInsufficientData:
    def test_single_sample_direction(self) -> None:
        t = _trend(_rec("S1", prrsv=True))
        assert t.direction is SwineRespTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_event_flags_are_false(self) -> None:
        t = _trend(_rec("S1", prrsv=True))
        assert t.prrsv_newly_detected is False
        assert t.prrsv_cleared is False

    def test_single_negative_is_insufficient(self) -> None:
        t = _trend(_rec("S1"))
        assert t.direction is SwineRespTrendDirection.INSUFFICIENT_DATA


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_is_stable_clear(self) -> None:
        t = _trend(_rec("S1"), _rec("S2"), _rec("S3"))
        assert t.direction is SwineRespTrendDirection.STABLE_CLEAR

    def test_stable_clear_no_prrsv_events(self) -> None:
        t = _trend(_rec("S1"), _rec("S2"))
        assert t.prrsv_newly_detected is False
        assert t.prrsv_cleared is False


# ---------------------------------------------------------------------------
# PRRSV newly detected -> EMERGING
# ---------------------------------------------------------------------------


class TestPrrsvNewlyDetected:
    def test_prrsv_new_is_emerging(self) -> None:
        t = _trend(_rec("S1"), _rec("S2", prrsv=True))
        assert t.direction is SwineRespTrendDirection.EMERGING
        assert t.prrsv_newly_detected is True
        assert t.prrsv_cleared is False

    def test_prrsv_always_negative_then_new_is_emerging(self) -> None:
        t = _trend(
            _rec("S1", mhp=True),
            _rec("S2", mhp=True),
            _rec("S3", prrsv=True, mhp=True),
        )
        assert t.direction is SwineRespTrendDirection.EMERGING
        assert t.prrsv_newly_detected is True


# ---------------------------------------------------------------------------
# PRRSV cleared -> RESOLVING
# ---------------------------------------------------------------------------


class TestPrrsvCleared:
    def test_prrsv_cleared_is_resolving(self) -> None:
        t = _trend(_rec("S1", prrsv=True), _rec("S2"))
        assert t.direction is SwineRespTrendDirection.RESOLVING
        assert t.prrsv_cleared is True
        assert t.prrsv_newly_detected is False

    def test_prrsv_cleared_with_mhp_remaining(self) -> None:
        t = _trend(
            _rec("S1", prrsv=True, mhp=True),
            _rec("S2", mhp=True),
        )
        assert t.direction is SwineRespTrendDirection.RESOLVING
        assert t.prrsv_cleared is True


# ---------------------------------------------------------------------------
# Burden-based trend classification
# ---------------------------------------------------------------------------


class TestBurdenTrend:
    def test_increasing_burden_is_emerging(self) -> None:
        t = _trend(
            _rec("S1"),
            _rec("S2", mhp=True),
            _rec("S3", pcv2=True, mhp=True),
        )
        assert t.direction is SwineRespTrendDirection.EMERGING

    def test_decreasing_burden_is_resolving(self) -> None:
        t = _trend(
            _rec("S1", pcv2=True, mhp=True),
            _rec("S2", mhp=True),
            _rec("S3"),
        )
        assert t.direction is SwineRespTrendDirection.RESOLVING

    def test_stable_endemic_with_consistent_mhp(self) -> None:
        t = _trend(
            _rec("S1", mhp=True),
            _rec("S2", mhp=True),
            _rec("S3", mhp=True),
        )
        assert t.direction is SwineRespTrendDirection.STABLE_ENDEMIC


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_worst_level_reflects_highest_seen(self) -> None:
        from lamp_forge.swine_respiratory_risk import PrdcAlertLevel

        t = _trend(_rec("S1", prrsv=True, pcv2=True), _rec("S2"))
        assert t.worst_alert_level is PrdcAlertLevel.CRITICAL

    def test_all_negative_worst_is_negative(self) -> None:
        from lamp_forge.swine_respiratory_risk import PrdcAlertLevel

        t = _trend(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is PrdcAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# Text fields are non-empty for all combinations
# ---------------------------------------------------------------------------


class TestTextFields:
    _DIRECTIONS: ClassVar[list[SwineRespTrendDirection]] = [
        SwineRespTrendDirection.EMERGING,
        SwineRespTrendDirection.RESOLVING,
        SwineRespTrendDirection.STABLE_CLEAR,
        SwineRespTrendDirection.STABLE_ENDEMIC,
        SwineRespTrendDirection.INSUFFICIENT_DATA,
    ]

    def test_interpretation_nonempty(self) -> None:
        cases = [
            [_rec("S1")],
            [_rec("S1"), _rec("S2")],
            [_rec("S1", prrsv=True), _rec("S2", prrsv=True)],
            [_rec("S1", prrsv=True), _rec("S2")],
            [_rec("S1"), _rec("S2", prrsv=True)],
            [_rec("S1", mhp=True), _rec("S2", mhp=True)],
        ]
        for recs in cases:
            t = analyse_swine_resp_trend(recs)
            assert t.interpretation
            assert t.recommended_action


# ---------------------------------------------------------------------------
# records_from_csv
# ---------------------------------------------------------------------------


class TestRecordsFromCsv:
    def test_basic_csv_parsing(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.write_text(
            "sample_id,date,prrsv,pcv2,mhp\n"
            "BarnA-Jan,2026-01-15,1,0,0\n"
            "BarnA-Feb,2026-02-15,1,1,0\n"
        )
        recs = records_from_csv(csv_path)
        assert len(recs) == 2
        assert recs[0].sample_id == "BarnA-Jan"
        assert recs[0].flags.prrsv is True
        assert recs[0].flags.pcv2 is False
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.pcv2 is True

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "bad.csv"
        csv_path.write_text("sample_id,prrsv,pcv2\nS1,1,0\n")
        with pytest.raises(ValueError, match="mhp"):
            records_from_csv(csv_path)

    def test_truthy_values(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "truthy.csv"
        csv_path.write_text("sample_id,prrsv,pcv2,mhp\nS1,true,false,yes\nS2,0,1,0\n")
        recs = records_from_csv(csv_path)
        assert recs[0].flags.prrsv is True
        assert recs[0].flags.pcv2 is False
        assert recs[0].flags.mhp is True
        assert recs[1].flags.prrsv is False
        assert recs[1].flags.pcv2 is True

    def test_date_optional(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "nodate.csv"
        csv_path.write_text("sample_id,prrsv,pcv2,mhp\nS1,0,0,1\n")
        recs = records_from_csv(csv_path)
        assert recs[0].date is None


# ---------------------------------------------------------------------------
# records_from_resp_json
# ---------------------------------------------------------------------------


class TestRecordsFromRespJson:
    def test_load_single_json(self, tmp_path: Path) -> None:
        j = tmp_path / "2026-01.json"
        j.write_text(
            json.dumps(
                {
                    "alert_level": "HIGH",
                    "panel_flags": {"prrsv": True, "pcv2": False, "mhp": False},
                    "date": "2026-01-15",
                }
            )
        )
        recs = records_from_resp_json([j])
        assert len(recs) == 1
        assert recs[0].flags.prrsv is True
        assert recs[0].date == "2026-01-15"
        assert recs[0].sample_id == "2026-01"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        j = tmp_path / "bad.json"
        j.write_text(json.dumps({"alert_level": "NEGATIVE"}))
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_resp_json([j])


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_json_output_valid(self, tmp_path: Path) -> None:
        t = _trend(_rec("S1", prrsv=True), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"
        assert data["prrsv_cleared"] is True
        assert len(data["timeline"]) == 2


# ---------------------------------------------------------------------------
# write_trend_csv
# ---------------------------------------------------------------------------


class TestWriteTrendCsv:
    def test_csv_header_and_rows(self, tmp_path: Path) -> None:
        t = _trend(_rec("S1", mhp=True), _rec("S2", mhp=True))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["sample_id", "date", "prrsv", "pcv2", "mhp", "alert_level"]
        assert len(rows) == 3  # header + 2 data rows
        assert rows[1][0] == "S1"
        assert rows[1][4] == "1"  # mhp = True -> 1
        assert rows[1][5] == "LOW"
