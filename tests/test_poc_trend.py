"""Tests for poc_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.poc_risk import PocAlertLevel, PocPanelFlags
from lamp_forge.poc_trend import (
    PocRecord,
    PocTrend,
    PocTrendDirection,
    analyse_poc_trend,
    records_from_csv,
    records_from_poc_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    mtb: bool = False,
    cdiff: bool = False,
    sars2: bool = False,
    flu_a: bool = False,
    gas: bool = False,
) -> PocPanelFlags:
    return PocPanelFlags(mtb=mtb, cdiff=cdiff, sars2=sars2, flu_a=flu_a, gas=gas)


def _rec(
    sample_id: str,
    flags: PocPanelFlags | None = None,
    date: str | None = None,
) -> PocRecord:
    return PocRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: PocRecord) -> PocTrend:
    return analyse_poc_trend(list(recs))


# ---------------------------------------------------------------------------
# Edge cases / insufficient data
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_poc_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(flu_a=True)))
        assert t.direction is PocTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_no_tb_flags(self) -> None:
        t = _analyse(_rec("S1", _flags(mtb=True)))
        assert t.tb_newly_detected is False
        assert t.tb_clearing is False


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is PocTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is PocTrendDirection.STABLE_CLEAR

    def test_stable_clear_interpretation_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        text = t.interpretation.lower()
        assert "negative" in text or "no active" in text

    def test_worst_alert_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is PocAlertLevel.NEGATIVE

    def test_respiratory_wave_false_when_all_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.respiratory_wave_active is False


# ---------------------------------------------------------------------------
# EMERGING
# ---------------------------------------------------------------------------


class TestEmerging:
    def test_tb_newly_detected_triggers_emerging(self) -> None:
        recs = [
            _rec("S1", _flags(flu_a=True)),
            _rec("S2", _flags(flu_a=True)),
            _rec("S3", _flags(flu_a=True, mtb=True)),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.EMERGING
        assert t.tb_newly_detected is True

    def test_tb_newly_detected_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(mtb=True)))
        assert t.direction is PocTrendDirection.EMERGING
        assert t.tb_newly_detected is True

    def test_tb_newly_detected_action_mentions_isolation(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(mtb=True)))
        text = t.recommended_action.lower()
        assert "isolation" in text or "notify" in text or "airborne" in text

    def test_burden_increasing_triggers_emerging(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3", _flags(flu_a=True)),
            _rec("S4", _flags(flu_a=True, sars2=True)),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.EMERGING

    def test_emerging_action_mentions_outbreak_or_notify(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(flu_a=True, sars2=True, cdiff=True)))
        text = t.recommended_action.lower()
        assert "outbreak" in text or "notify" in text or "isolation" in text

    def test_tb_onset_overrides_resolving_signal(self) -> None:
        # Even if earlier sample had more pathogens, TB onset in last = EMERGING
        recs = [
            _rec("S1", _flags(flu_a=True, sars2=True)),
            _rec("S2", _flags(mtb=True)),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.EMERGING
        assert t.tb_newly_detected is True


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_last_sample_negative_after_positives(self) -> None:
        t = _analyse(_rec("S1", _flags(flu_a=True)), _rec("S2"))
        assert t.direction is PocTrendDirection.RESOLVING

    def test_burden_decreasing_triggers_resolving(self) -> None:
        recs = [
            _rec("S1", _flags(flu_a=True, sars2=True)),
            _rec("S2", _flags(flu_a=True, sars2=True)),
            _rec("S3"),
            _rec("S4"),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.RESOLVING

    def test_tb_clearing_flag_set(self) -> None:
        recs = [
            _rec("S1", _flags(mtb=True)),
            _rec("S2", _flags(mtb=True)),
            _rec("S3"),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.RESOLVING
        assert t.tb_clearing is True

    def test_resolving_action_mentions_maintain_or_surveillance(self) -> None:
        t = _analyse(_rec("S1", _flags(flu_a=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "maintain" in text or "continue" in text or "surveillance" in text


# ---------------------------------------------------------------------------
# STABLE_ENDEMIC
# ---------------------------------------------------------------------------


class TestStableEndemic:
    def test_consistently_flu_positive(self) -> None:
        recs = [
            _rec("S1", _flags(flu_a=True)),
            _rec("S2", _flags(flu_a=True)),
            _rec("S3", _flags(flu_a=True)),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.STABLE_ENDEMIC

    def test_stable_endemic_interpretation_non_empty(self) -> None:
        recs = [_rec("S1", _flags(gas=True)), _rec("S2", _flags(gas=True))]
        t = analyse_poc_trend(recs)
        assert t.interpretation.strip()

    def test_respiratory_wave_active_when_recent_half_dominant(self) -> None:
        recs = [
            _rec("S1", _flags(flu_a=True)),
            _rec("S2", _flags(flu_a=True)),
            _rec("S3", _flags(flu_a=True)),
            _rec("S4", _flags(flu_a=True)),
        ]
        t = analyse_poc_trend(recs)
        assert t.respiratory_wave_active is True

    def test_respiratory_wave_note_in_interpretation(self) -> None:
        recs = [
            _rec("S1", _flags(flu_a=True, sars2=True)),
            _rec("S2", _flags(flu_a=True, sars2=True)),
        ]
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.STABLE_ENDEMIC
        text = t.interpretation.lower()
        assert "respiratory" in text or "wave" in text or "seasonal" in text


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is PocAlertLevel.NEGATIVE

    def test_mtb_detection_gives_critical_worst(self) -> None:
        recs = [_rec("S1"), _rec("S2", _flags(mtb=True)), _rec("S3")]
        t = analyse_poc_trend(recs)
        assert t.worst_alert_level is PocAlertLevel.CRITICAL

    def test_gas_only_gives_low_worst(self) -> None:
        recs = [_rec("S1", _flags(gas=True)), _rec("S2", _flags(gas=True))]
        t = analyse_poc_trend(recs)
        assert t.worst_alert_level is PocAlertLevel.LOW

    def test_worst_preserved_after_clearance(self) -> None:
        recs = [
            _rec("S1", _flags(cdiff=True)),
            _rec("S2"),
        ]
        t = analyse_poc_trend(recs)
        assert t.worst_alert_level is PocAlertLevel.HIGH


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(flu_a=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "tb_newly_detected",
            "tb_clearing",
            "respiratory_wave_active",
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
        t = _analyse(_rec("S1", _flags(flu_a=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
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
                    "sample_id": "Clinic-Jan",
                    "date": "2026-01-15",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "1",
                    "flu_a": "1",
                    "gas": "0",
                },
                {
                    "sample_id": "Clinic-Feb",
                    "date": "2026-02-15",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "0",
                    "gas": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "Clinic-Jan"
        assert recs[0].flags.sars2 is True
        assert recs[0].flags.flu_a is True
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.sars2 is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "mtb": "true",
                    "cdiff": "false",
                    "sars2": "TRUE",
                    "flu_a": "FALSE",
                    "gas": "1",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.mtb is True
        assert recs[0].flags.cdiff is False
        assert recs[0].flags.sars2 is True
        assert recs[0].flags.flu_a is False
        assert recs[0].flags.gas is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "0",
                    "gas": "0",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "mtb": "0", "cdiff": "0", "sars2": "0", "flu_a": "0"}],
        )
        with pytest.raises(ValueError, match="gas"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "A",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "1",
                    "gas": "0",
                },
                {
                    "sample_id": "B",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "1",
                    "gas": "0",
                },
                {
                    "sample_id": "C",
                    "mtb": "1",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "1",
                    "gas": "0",
                },
            ],
        )
        recs = records_from_csv(p)
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.EMERGING
        assert t.tb_newly_detected is True


# ---------------------------------------------------------------------------
# records_from_poc_json
# ---------------------------------------------------------------------------


class TestRecordsFromPocJson:
    def _write_poc_json(
        self,
        path: Path,
        mtb: bool = False,
        cdiff: bool = False,
        sars2: bool = False,
        flu_a: bool = False,
        gas: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "NEGATIVE",
            "panel_flags": {
                "mtb": mtb,
                "cdiff": cdiff,
                "sars2": sars2,
                "flu_a": flu_a,
                "gas": gas,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_poc_json(p1, flu_a=True)
        self._write_poc_json(p2)
        recs = records_from_poc_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.flu_a is True
        assert recs[1].flags.flu_a is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_poc_json(p, date="2026-03-15")
        recs = records_from_poc_json([p])
        assert recs[0].date == "2026-03-15"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_poc_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, (flu_a, mtb) in enumerate([(True, False), (True, False), (False, True)]):
            p = tmp_path / f"s{i}.json"
            self._write_poc_json(p, flu_a=flu_a, mtb=mtb)
            files.append(p)
        recs = records_from_poc_json(files)
        t = analyse_poc_trend(recs)
        assert t.direction is PocTrendDirection.EMERGING
        assert t.tb_newly_detected is True


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(flu_a=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(flu_a=True)), _rec("S2"))
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
        assert rows[0] == [
            "sample_id",
            "date",
            "mtb",
            "cdiff",
            "sars2",
            "flu_a",
            "gas",
            "alert_level",
        ]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(mtb=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_flu_a_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(flu_a=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["flu_a"] == "1"
        assert dr[1]["flu_a"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestPocTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["poc-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        flu_a: bool = False,
        sars2: bool = False,
        mtb: bool = False,
        gas: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "NEGATIVE",
            "panel_flags": {
                "mtb": mtb,
                "cdiff": False,
                "sars2": sars2,
                "flu_a": flu_a,
                "gas": gas,
            },
        }
        with path.open("w") as fh:
            json.dump(data, fh)

    def _write_csv(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_poc_result_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, flu_a=True)
        result = self._run(["--poc-result", str(p1), "--poc-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, mtb=True)
        result = self._run(["--poc-result", str(p1), "--poc-result", str(p2)])
        assert "EMERGING" in result.output  # type: ignore[union-attr]

    def test_tb_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, mtb=True)
        result = self._run(["--poc-result", str(p1), "--poc-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "TB NEWLY DETECTED" in output or "tb" in output.lower()

    def test_resolving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, flu_a=True)
        self._write_json(p2)
        result = self._run(["--poc-result", str(p1), "--poc-result", str(p2)])
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2, p3 = tmp_path / "s1.json", tmp_path / "s2.json", tmp_path / "s3.json"
        self._write_json(p1)
        self._write_json(p2)
        self._write_json(p3)
        result = self._run(
            [
                "--poc-result",
                str(p1),
                "--poc-result",
                str(p2),
                "--poc-result",
                str(p3),
            ]
        )
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, flu_a=True)
        result = self._run(["--poc-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {
                    "sample_id": "A",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "0",
                    "gas": "0",
                },
                {
                    "sample_id": "B",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "1",
                    "flu_a": "1",
                    "gas": "0",
                },
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, mtb=True)
        out = tmp_path / "trend.json"
        result = self._run(
            ["--poc-result", str(p1), "--poc-result", str(p2), "--out-json", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "EMERGING"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, flu_a=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            ["--poc-result", str(p1), "--poc-result", str(p2), "--out-csv", str(out)]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_no_input_exits_nonzero(self) -> None:
        result = self._run([])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_csv_and_poc_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [
                {
                    "sample_id": "A",
                    "mtb": "0",
                    "cdiff": "0",
                    "sars2": "0",
                    "flu_a": "0",
                    "gas": "0",
                }
            ],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--poc-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "clinic_a_jan.json"
        p2 = tmp_path / "clinic_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, flu_a=True)
        result = self._run(["--poc-result", str(p1), "--poc-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "clinic_a_jan" in output
        assert "clinic_a_feb" in output
