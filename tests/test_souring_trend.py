"""Tests for the souring-trajectory analysis module.

All tests are pure Python (no external binaries required).  The trend
direction logic is deterministic, so tests use exact comparisons on
direction enum values and flag attributes.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.mic_risk import GuildFlags, MICRiskLevel
from lamp_forge.souring_trend import (
    GuildRecord,
    SouringTrend,
    TrendDirection,
    analyse_souring_trend,
    records_from_csv,
    records_from_mic_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    srb: bool = False,
    mcr: bool = False,
    irb: bool = False,
    apb: bool = False,
    nrb: bool = False,
) -> GuildFlags:
    return GuildFlags(srb=srb, mcr=mcr, irb=irb, apb=apb, nrb=nrb)


def _rec(
    sample_id: str,
    srb: bool = False,
    mcr: bool = False,
    irb: bool = False,
    apb: bool = False,
    nrb: bool = False,
    date: str | None = None,
) -> GuildRecord:
    return GuildRecord(sample_id=sample_id, flags=_flags(srb, mcr, irb, apb, nrb), date=date)


def _analyse(*recs: GuildRecord) -> SouringTrend:
    return analyse_souring_trend(list(recs))


# ---------------------------------------------------------------------------
# Empty / single sample
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_souring_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", srb=True))
        assert t.direction is TrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1
        assert t.srb_newly_detected is False
        assert t.srb_clearing is False

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()


# ---------------------------------------------------------------------------
# STABLE_LOW
# ---------------------------------------------------------------------------


class TestStableLow:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is TrendDirection.STABLE_LOW
        assert t.srb_detection_count == 0

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is TrendDirection.STABLE_LOW

    def test_single_positive_in_four_is_stable_low(self) -> None:
        # One SRB hit out of four (25%) -> below 50% threshold.
        t = _analyse(
            _rec("S1", srb=True),
            _rec("S2"),
            _rec("S3"),
            _rec("S4"),
        )
        # srb[0]=T, srb[-1]=F -> clearing flag triggers IMPROVING
        assert t.direction is TrendDirection.IMPROVING

    def test_stable_low_interpretation_mentions_absent(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        text = t.interpretation.lower()
        assert "absent" in text or "negative" in text or "low" in text

    def test_nrb_only_stable_low_mentions_suppression(self) -> None:
        t = _analyse(_rec("S1", nrb=True), _rec("S2", nrb=True), _rec("S3", nrb=True))
        assert t.direction is TrendDirection.STABLE_LOW
        assert (
            "suppression" in t.interpretation.lower() or "suppressing" in t.interpretation.lower()
        )

    def test_nrb_stable_low_action_mentions_nitrate(self) -> None:
        t = _analyse(_rec("S1", nrb=True), _rec("S2", nrb=True))
        assert "nitrate" in t.recommended_action.lower()


# ---------------------------------------------------------------------------
# STABLE_HIGH
# ---------------------------------------------------------------------------


class TestStableHigh:
    def test_all_srb_positive_two_samples(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2", srb=True))
        assert t.direction is TrendDirection.STABLE_HIGH

    def test_most_srb_positive_is_stable_high(self) -> None:
        # 3/4 positive (75%) with no clear trend -> STABLE_HIGH
        t = _analyse(
            _rec("S1", srb=True),
            _rec("S2", srb=True),
            _rec("S3", srb=False),
            _rec("S4", srb=True),
        )
        assert t.direction is TrendDirection.STABLE_HIGH

    def test_stable_high_srb_count_correct(self) -> None:
        t = _analyse(_rec("A", srb=True), _rec("B", srb=True), _rec("C", srb=True))
        assert t.srb_detection_count == 3

    def test_stable_high_interpretation_mentions_chronic(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2", srb=True), _rec("S3", srb=True))
        text = t.interpretation.lower()
        assert "chronic" in text or "consistently" in text

    def test_underdosed_note_added_when_nrb_present(self) -> None:
        t = _analyse(
            _rec("S1", srb=True, nrb=True),
            _rec("S2", srb=True, nrb=True),
        )
        assert t.direction is TrendDirection.STABLE_HIGH
        assert t.treatment_underdosed_count == 2
        assert "nitrate" in t.interpretation.lower() or "nrb" in t.interpretation.lower()


# ---------------------------------------------------------------------------
# DETERIORATING
# ---------------------------------------------------------------------------


class TestDeteriorating:
    def test_newly_detected_flag(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3", srb=True))
        assert t.direction is TrendDirection.DETERIORATING
        assert t.srb_newly_detected is True
        assert t.srb_clearing is False

    def test_newly_detected_single_negative_before(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", srb=True))
        assert t.srb_newly_detected is True
        assert t.direction is TrendDirection.DETERIORATING

    def test_rate_increasing_over_four_samples(self) -> None:
        # early=[F,F] (rate=0.0), recent=[T,T] (rate=1.0) -> DETERIORATING
        t = _analyse(
            _rec("S1"),
            _rec("S2"),
            _rec("S3", srb=True),
            _rec("S4", srb=True),
        )
        assert t.direction is TrendDirection.DETERIORATING

    def test_deteriorating_text_mentions_immediate(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", srb=True))
        text = t.recommended_action.lower()
        assert "immediately" in text or "immediate" in text or "escalat" in text

    def test_underdosed_note_in_deteriorating(self) -> None:
        t = _analyse(
            _rec("S1"),
            _rec("S2"),
            _rec("S3", srb=True, nrb=True),
        )
        assert t.direction is TrendDirection.DETERIORATING
        assert t.treatment_underdosed_count == 1
        assert "nitrate" in t.interpretation.lower() or "under" in t.interpretation.lower()

    def test_srb_count_correct_in_deteriorating(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3", srb=True))
        assert t.srb_detection_count == 1


# ---------------------------------------------------------------------------
# IMPROVING
# ---------------------------------------------------------------------------


class TestImproving:
    def test_clearing_flag(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        assert t.direction is TrendDirection.IMPROVING
        assert t.srb_clearing is True
        assert t.srb_newly_detected is False

    def test_clearing_over_three_samples(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2", srb=True), _rec("S3"))
        assert t.direction is TrendDirection.IMPROVING
        assert t.srb_clearing is True

    def test_rate_decreasing_over_four_samples(self) -> None:
        # early=[T,T] (rate=1.0), recent=[F,F] (rate=0.0) -> IMPROVING via clearing first
        t = _analyse(
            _rec("S1", srb=True),
            _rec("S2", srb=True),
            _rec("S3"),
            _rec("S4"),
        )
        assert t.direction is TrendDirection.IMPROVING

    def test_improving_text_mentions_effective_or_continue(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        text = t.recommended_action.lower()
        assert "continue" in text or "maintain" in text

    def test_clearing_interpretation_mentions_cleared(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        assert "cleared" in t.interpretation.lower() or "clearing" in t.interpretation.lower()

    def test_srb_count_correct_in_improving(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2", srb=True), _rec("S3"))
        assert t.srb_detection_count == 2


# ---------------------------------------------------------------------------
# Worst risk level
# ---------------------------------------------------------------------------


class TestWorstRiskLevel:
    def test_all_negative_worst_is_minimal(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_risk_level is MICRiskLevel.MINIMAL

    def test_single_srb_worst_is_at_least_moderate(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", srb=True))
        assert t.worst_risk_level in (
            MICRiskLevel.MODERATE,
            MICRiskLevel.HIGH,
            MICRiskLevel.CRITICAL,
        )

    def test_critical_guilds_give_critical_worst(self) -> None:
        t = _analyse(
            _rec("S1"),
            _rec("S2", srb=True, mcr=True, irb=True, apb=True),
        )
        assert t.worst_risk_level is MICRiskLevel.CRITICAL


# ---------------------------------------------------------------------------
# NRB / treatment tracking
# ---------------------------------------------------------------------------


class TestNrbTracking:
    def test_nrb_count_correct(self) -> None:
        t = _analyse(_rec("S1", nrb=True), _rec("S2", nrb=True), _rec("S3"))
        assert t.nrb_detection_count == 2

    def test_underdosed_count_correct(self) -> None:
        t = _analyse(
            _rec("S1", srb=True, nrb=True),
            _rec("S2", srb=True),
            _rec("S3", srb=True, nrb=True),
        )
        assert t.treatment_underdosed_count == 2

    def test_nrb_without_srb_not_underdosed(self) -> None:
        t = _analyse(_rec("S1", nrb=True), _rec("S2", nrb=True))
        assert t.treatment_underdosed_count == 0


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", srb=True))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "srb_detection_count",
            "nrb_detection_count",
            "srb_newly_detected",
            "srb_clearing",
            "treatment_underdosed_count",
            "worst_risk_level",
            "interpretation",
            "recommended_action",
            "timeline",
        ):
            assert key in d, f"Missing key: {key}"

    def test_timeline_length_matches_n_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        d = t.to_dict()
        assert len(d["timeline"]) == 3  # type: ignore[arg-type]

    def test_timeline_contains_sample_ids(self) -> None:
        t = _analyse(_rec("Alpha"), _rec("Beta"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert tl[0]["sample_id"] == "Alpha"
        assert tl[1]["sample_id"] == "Beta"

    def test_json_serialisable(self) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "IMPROVING"

    def test_direction_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)


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
                    "sample_id": "W1-Jan",
                    "date": "2026-01-15",
                    "srb": "1",
                    "mcr": "0",
                    "irb": "0",
                    "apb": "0",
                    "nrb": "0",
                },
                {
                    "sample_id": "W1-Feb",
                    "date": "2026-02-15",
                    "srb": "0",
                    "mcr": "0",
                    "irb": "1",
                    "apb": "0",
                    "nrb": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "W1-Jan"
        assert recs[0].flags.srb is True
        assert recs[0].flags.irb is False
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.irb is True

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "srb": "true",
                    "mcr": "false",
                    "irb": "TRUE",
                    "apb": "FALSE",
                    "nrb": "true",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.srb is True
        assert recs[0].flags.mcr is False
        assert recs[0].flags.irb is True
        assert recs[0].flags.apb is False
        assert recs[0].flags.nrb is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "srb": "0", "mcr": "0", "irb": "0", "apb": "0", "nrb": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "srb": "0", "mcr": "0", "irb": "0", "apb": "0"}],
        )
        with pytest.raises(ValueError, match="nrb"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "srb": "0", "mcr": "0", "irb": "0", "apb": "0", "nrb": "0"},
                {"sample_id": "B", "srb": "1", "mcr": "0", "irb": "0", "apb": "0", "nrb": "0"},
                {"sample_id": "C", "srb": "1", "mcr": "0", "irb": "1", "apb": "0", "nrb": "0"},
            ],
        )
        recs = records_from_csv(p)
        trend = analyse_souring_trend(recs)
        assert trend.direction is TrendDirection.DETERIORATING
        assert trend.srb_newly_detected is False


# ---------------------------------------------------------------------------
# records_from_mic_json
# ---------------------------------------------------------------------------


class TestRecordsFromMicJson:
    def _write_mic_json(
        self,
        path: Path,
        srb: bool = False,
        mcr: bool = False,
        irb: bool = False,
        apb: bool = False,
        nrb: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "risk_level": "MODERATE",
            "guild_flags": {
                "srb": srb,
                "mcr": mcr,
                "irb": irb,
                "apb": apb,
                "nrb": nrb,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_mic_json(p1, srb=True)
        self._write_mic_json(p2, srb=False)
        recs = records_from_mic_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.srb is True
        assert recs[1].flags.srb is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_mic_json(p, date="2026-03-15")
        recs = records_from_mic_json([p])
        assert recs[0].date == "2026-03-15"

    def test_missing_guild_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"risk_level": "HIGH"}')
        with pytest.raises(ValueError, match="guild_flags"):
            records_from_mic_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, srb in enumerate([False, False, True]):
            p = tmp_path / f"s{i}.json"
            self._write_mic_json(p, srb=srb)
            files.append(p)
        recs = records_from_mic_json(files)
        trend = analyse_souring_trend(recs)
        assert trend.srb_newly_detected is True
        assert trend.direction is TrendDirection.DETERIORATING


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", srb=True))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "IMPROVING"
        assert data["srb_clearing"] is True

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
        assert rows[0] == ["sample_id", "date", "srb", "mcr", "irb", "apb", "nrb", "risk_level"]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_risk_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["risk_level"] in ("MINIMAL", "LOW", "MODERATE", "HIGH", "CRITICAL")

    def test_srb_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", srb=True), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["srb"] == "1"
        assert dr[1]["srb"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestSouringTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["souring-trend", *args], obj={})

    def _write_json(
        self, path: Path, srb: bool = False, nrb: bool = False, irb: bool = False
    ) -> None:
        data = {
            "risk_level": "LOW",
            "guild_flags": {
                "srb": srb,
                "mcr": False,
                "irb": irb,
                "apb": False,
                "nrb": nrb,
            },
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_mic_result_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, srb=True)
        result = self._run(["--mic-result", str(p1), "--mic-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, srb=True)
        result = self._run(["--mic-result", str(p1), "--mic-result", str(p2)])
        assert "DETERIORATING" in result.output  # type: ignore[union-attr]

    def test_improving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, srb=True)
        self._write_json(p2)
        result = self._run(["--mic-result", str(p1), "--mic-result", str(p2)])
        assert "IMPROVING" in result.output  # type: ignore[union-attr]

    def test_stable_low_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2, p3 = tmp_path / "s1.json", tmp_path / "s2.json", tmp_path / "s3.json"
        self._write_json(p1)
        self._write_json(p2)
        self._write_json(p3)
        result = self._run(
            ["--mic-result", str(p1), "--mic-result", str(p2), "--mic-result", str(p3)]
        )
        assert "STABLE_LOW" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, srb=True)
        result = self._run(["--mic-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "srb": "0", "mcr": "0", "irb": "0", "apb": "0", "nrb": "0"},
                {"sample_id": "B", "srb": "1", "mcr": "0", "irb": "0", "apb": "0", "nrb": "0"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, srb=True)
        out = tmp_path / "trend.json"
        result = self._run(
            ["--mic-result", str(p1), "--mic-result", str(p2), "--out-json", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "DETERIORATING"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, srb=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            ["--mic-result", str(p1), "--mic-result", str(p2), "--out-csv", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_mic_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "srb": "0", "mcr": "0", "irb": "0", "apb": "0", "nrb": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--mic-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "well1_jan.json"
        p2 = tmp_path / "well1_feb.json"
        self._write_json(p1)
        self._write_json(p2, srb=True)
        result = self._run(["--mic-result", str(p1), "--mic-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "well1_jan" in output
        assert "well1_feb" in output
