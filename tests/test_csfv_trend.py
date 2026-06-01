"""Tests for csfv_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.csfv_risk import CsfvAlertLevel, CsfvContext, CsfvFlags
from lamp_forge.csfv_trend import (
    CsfvRecord,
    CsfvTrend,
    CsfvTrendDirection,
    analyse_csfv_trend,
    records_from_csfv_json,
    records_from_csv,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    ns5b: bool = False,
    context: CsfvContext = CsfvContext.ON_FARM,
) -> CsfvFlags:
    return CsfvFlags(ns5b_positive=ns5b, context=context)


def _rec(
    sample_id: str,
    flags: CsfvFlags | None = None,
    date: str | None = None,
) -> CsfvRecord:
    return CsfvRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: CsfvRecord) -> CsfvTrend:
    return analyse_csfv_trend(list(recs))


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_csfv_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)))
        assert t.direction is CsfvTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_event_flags_false(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)))
        assert t.newly_detected is False
        assert t.controlled is False
        assert t.sustained_outbreak is False

    def test_single_positive_worst_level_critical(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)))
        assert t.worst_alert_level is CsfvAlertLevel.CRITICAL

    def test_single_negative_worst_level_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=False)))
        assert t.worst_alert_level is CsfvAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is CsfvTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is CsfvTrendDirection.STABLE_CLEAR

    def test_stable_clear_text_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        assert "negative" in t.interpretation.lower() or "not detected" in t.interpretation.lower()

    def test_worst_alert_negative_when_all_clear(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is CsfvAlertLevel.NEGATIVE

    def test_no_event_flags_when_clear(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.newly_detected is False
        assert t.controlled is False
        assert t.sustained_outbreak is False


# ---------------------------------------------------------------------------
# NEWLY_DETECTED
# ---------------------------------------------------------------------------


class TestNewlyDetected:
    def test_first_positive_in_last_sample(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        assert t.direction is CsfvTrendDirection.NEWLY_DETECTED
        assert t.newly_detected is True

    def test_first_positive_after_three_negatives(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3"),
            _rec("S4", _flags(ns5b=True)),
        ]
        t = analyse_csfv_trend(recs)
        assert t.direction is CsfvTrendDirection.NEWLY_DETECTED
        assert t.newly_detected is True

    def test_newly_detected_action_mentions_quarantine(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        text = t.recommended_action.lower()
        assert "quarantine" in text or "notify" in text

    def test_newly_detected_interpretation_mentions_first_detection(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        text = t.interpretation.lower()
        assert "newly" in text or "first" in text or "new" in text

    def test_not_newly_detected_when_prior_positive(self) -> None:
        t = _analyse(
            _rec("S1", _flags(ns5b=True)),
            _rec("S2", _flags(ns5b=True)),
        )
        assert t.newly_detected is False

    def test_newly_detected_worst_level_critical(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        assert t.worst_alert_level is CsfvAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# ONGOING_OUTBREAK
# ---------------------------------------------------------------------------


class TestOngoingOutbreak:
    def test_two_consecutive_positives(self) -> None:
        t = _analyse(
            _rec("S1", _flags(ns5b=True)),
            _rec("S2", _flags(ns5b=True)),
        )
        assert t.direction is CsfvTrendDirection.ONGOING_OUTBREAK

    def test_positive_then_negative_then_positive(self) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True)),
            _rec("S2"),
            _rec("S3", _flags(ns5b=True)),
        ]
        t = analyse_csfv_trend(recs)
        assert t.direction is CsfvTrendDirection.ONGOING_OUTBREAK

    def test_sustained_outbreak_true_for_high_ratio(self) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True)),
            _rec("S2", _flags(ns5b=True)),
            _rec("S3", _flags(ns5b=True)),
        ]
        t = analyse_csfv_trend(recs)
        assert t.sustained_outbreak is True

    def test_sustained_outbreak_false_for_only_two_samples(self) -> None:
        t = _analyse(
            _rec("S1", _flags(ns5b=True)),
            _rec("S2", _flags(ns5b=True)),
        )
        assert t.sustained_outbreak is False

    def test_sustained_outbreak_false_for_low_ratio(self) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True)),
            _rec("S2"),
            _rec("S3"),
            _rec("S4", _flags(ns5b=True)),
        ]
        t = analyse_csfv_trend(recs)
        assert t.direction is CsfvTrendDirection.ONGOING_OUTBREAK
        assert t.sustained_outbreak is False

    def test_ongoing_action_mentions_quarantine(self) -> None:
        recs = [_rec("S1", _flags(ns5b=True)), _rec("S2", _flags(ns5b=True))]
        t = analyse_csfv_trend(recs)
        assert (
            "quarantine" in t.recommended_action.lower()
            or "movement" in t.recommended_action.lower()
        )

    def test_ongoing_interpretation_mentions_positive(self) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True)),
            _rec("S2", _flags(ns5b=True)),
            _rec("S3", _flags(ns5b=True)),
        ]
        t = analyse_csfv_trend(recs)
        assert "positive" in t.interpretation.lower() or "outbreak" in t.interpretation.lower()


# ---------------------------------------------------------------------------
# CONTROLLED
# ---------------------------------------------------------------------------


class TestControlled:
    def test_positive_then_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        assert t.direction is CsfvTrendDirection.CONTROLLED
        assert t.controlled is True

    def test_multiple_positives_then_negative(self) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True)),
            _rec("S2", _flags(ns5b=True)),
            _rec("S3"),
        ]
        t = analyse_csfv_trend(recs)
        assert t.direction is CsfvTrendDirection.CONTROLLED
        assert t.controlled is True

    def test_controlled_action_mentions_quarantine_or_consecutive(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "quarantine" in text or "consecutive" in text or "serology" in text

    def test_controlled_flag_set(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        assert t.controlled is True
        assert t.newly_detected is False

    def test_worst_level_preserved_after_clearance(self) -> None:
        t = _analyse(
            _rec("S1", _flags(ns5b=True)),
            _rec("S2"),
        )
        assert t.worst_alert_level is CsfvAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is CsfvAlertLevel.NEGATIVE

    def test_any_positive_worst_is_critical(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        assert t.worst_alert_level is CsfvAlertLevel.CRITICAL

    def test_worst_carried_through_controlled(self) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True, context=CsfvContext.BORDER_POST)),
            _rec("S2", _flags(ns5b=False, context=CsfvContext.BORDER_POST)),
        ]
        t = analyse_csfv_trend(recs)
        assert t.worst_alert_level is CsfvAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "newly_detected",
            "controlled",
            "sustained_outbreak",
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
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "CONTROLLED"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_newly_detected_is_bool(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        d = t.to_dict()
        assert isinstance(d["newly_detected"], bool)
        assert d["newly_detected"] is True

    def test_timeline_has_alert_level_column(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert "alert_level" in tl[0]
        assert tl[0]["alert_level"] == "CRITICAL"

    def test_timeline_ns5b_values(self) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert tl[0]["ns5b_positive"] is True
        assert tl[1]["ns5b_positive"] is False


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
                    "ns5b_positive": "1",
                    "context": "on_farm",
                },
                {
                    "sample_id": "FarmA-Feb",
                    "date": "2026-02-15",
                    "ns5b_positive": "0",
                    "context": "on_farm",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "FarmA-Jan"
        assert recs[0].flags.ns5b_positive is True
        assert recs[0].flags.context is CsfvContext.ON_FARM
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.ns5b_positive is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "ns5b_positive": "true",
                    "context": "border_post",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.ns5b_positive is True
        assert recs[0].flags.context is CsfvContext.BORDER_POST

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "ns5b_positive": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_optional_context_defaults_to_on_farm(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "ns5b_positive": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.context is CsfvContext.ON_FARM

    def test_market_context_parsed(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "ns5b_positive": "1", "context": "market"}],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.context is CsfvContext.MARKET

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1"}],
        )
        with pytest.raises(ValueError, match="ns5b_positive"):
            records_from_csv(p)

    def test_invalid_context_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "ns5b_positive": "1", "context": "village"}],
        )
        with pytest.raises(ValueError, match="context"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "ns5b_positive": "0"},
                {"sample_id": "B", "ns5b_positive": "0"},
                {"sample_id": "C", "ns5b_positive": "1"},
            ],
        )
        recs = records_from_csv(p)
        t = analyse_csfv_trend(recs)
        assert t.direction is CsfvTrendDirection.NEWLY_DETECTED
        assert t.newly_detected is True


# ---------------------------------------------------------------------------
# records_from_csfv_json
# ---------------------------------------------------------------------------


class TestRecordsFromCsfvJson:
    def _write_csfv_json(
        self,
        path: Path,
        ns5b_positive: bool = False,
        context: str = "on_farm",
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "CRITICAL" if ns5b_positive else "NEGATIVE",
            "alert_score": 100 if ns5b_positive else 0,
            "ns5b_result": "positive" if ns5b_positive else "negative",
            "context": context,
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_csfv_json(p1, ns5b_positive=True)
        self._write_csfv_json(p2)
        recs = records_from_csfv_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.ns5b_positive is True
        assert recs[1].flags.ns5b_positive is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_csfv_json(p, date="2026-03-15")
        recs = records_from_csfv_json([p])
        assert recs[0].date == "2026-03-15"

    def test_border_post_context_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "border.json"
        self._write_csfv_json(p, context="border_post")
        recs = records_from_csfv_json([p])
        assert recs[0].flags.context is CsfvContext.BORDER_POST

    def test_missing_ns5b_result_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "CRITICAL"}')
        with pytest.raises(ValueError, match="ns5b_result"):
            records_from_csfv_json([p])

    def test_unknown_context_defaults_to_on_farm(self, tmp_path: Path) -> None:
        p = tmp_path / "unk.json"
        data = {"ns5b_result": "negative", "context": "warehouse"}
        p.write_text(json.dumps(data))
        recs = records_from_csfv_json([p])
        assert recs[0].flags.context is CsfvContext.ON_FARM

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, pos in enumerate([False, True, True]):
            p = tmp_path / f"s{i}.json"
            self._write_csfv_json(p, ns5b_positive=pos)
            files.append(p)
        recs = records_from_csfv_json(files)
        t = analyse_csfv_trend(recs)
        assert t.direction is CsfvTrendDirection.ONGOING_OUTBREAK


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(ns5b=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "CONTROLLED"
        assert data["controlled"] is True

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
            "ns5b_positive",
            "context",
            "alert_level",
        ]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_ns5b_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(ns5b=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["ns5b_positive"] == "1"
        assert dr[1]["ns5b_positive"] == "0"

    def test_context_column_values(self, tmp_path: Path) -> None:
        recs = [
            _rec("S1", _flags(ns5b=True, context=CsfvContext.MARKET)),
            _rec("S2", _flags(ns5b=False, context=CsfvContext.MARKET)),
        ]
        t = analyse_csfv_trend(recs)
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["context"] == "market"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestCsfvTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["csfv-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        ns5b_positive: bool = False,
        context: str = "on_farm",
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "CRITICAL" if ns5b_positive else "NEGATIVE",
            "alert_score": 100 if ns5b_positive else 0,
            "ns5b_result": "positive" if ns5b_positive else "negative",
            "context": context,
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ns5b_positive=True)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ns5b_positive=True)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        assert "NEWLY_DETECTED" in result.output  # type: ignore[union-attr]

    def test_controlled_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, ns5b_positive=True)
        self._write_json(p2)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        assert "CONTROLLED" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, ns5b_positive=True)
        result = self._run(["--csfv-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "ns5b_positive": "0"},
                {"sample_id": "B", "ns5b_positive": "1"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ns5b_positive=True)
        out = tmp_path / "trend.json"
        result = self._run(
            [
                "--csfv-result",
                str(p1),
                "--csfv-result",
                str(p2),
                "--out-json",
                str(out),
            ]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "NEWLY_DETECTED"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, ns5b_positive=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            [
                "--csfv-result",
                str(p1),
                "--csfv-result",
                str(p2),
                "--out-csv",
                str(out),
            ]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_csfv_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "ns5b_positive": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--csfv-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "farm_a_jan.json"
        p2 = tmp_path / "farm_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, ns5b_positive=True)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "farm_a_jan" in output
        assert "farm_a_feb" in output

    def test_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, ns5b_positive=True)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "NEWLY DETECTED" in output or "newly" in output.lower()

    def test_ongoing_outbreak_for_consecutive_positives(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, ns5b_positive=True)
        self._write_json(p2, ns5b_positive=True)
        result = self._run(["--csfv-result", str(p1), "--csfv-result", str(p2)])
        assert "ONGOING_OUTBREAK" in result.output  # type: ignore[union-attr]

    def test_border_post_context_accepted(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "border.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "BP-Jan", "ns5b_positive": "0", "context": "border_post"},
                {"sample_id": "BP-Feb", "ns5b_positive": "1", "context": "border_post"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "NEWLY_DETECTED" in result.output  # type: ignore[union-attr]
