"""Tests for bov_enteric_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.bov_enteric_risk import BovEntericAlertLevel, BovEntericFlags
from lamp_forge.bov_enteric_trend import (
    BovEntericRecord,
    BovEntericTrend,
    BovEntericTrendDirection,
    analyse_bov_enteric_trend,
    records_from_csv,
    records_from_enteric_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    etec_k99: bool = False,
    salmonella: bool = False,
    crypto: bool = False,
    bcov: bool = False,
    rota_a: bool = False,
) -> BovEntericFlags:
    return BovEntericFlags(
        etec_k99=etec_k99,
        salmonella=salmonella,
        crypto=crypto,
        bcov=bcov,
        rota_a=rota_a,
    )


def _rec(
    sample_id: str,
    flags: BovEntericFlags | None = None,
    date: str | None = None,
) -> BovEntericRecord:
    return BovEntericRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: BovEntericRecord) -> BovEntericTrend:
    return analyse_bov_enteric_trend(list(recs))


# ---------------------------------------------------------------------------
# Edge cases / insufficient data
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_bov_enteric_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(rota_a=True)))
        assert t.direction is BovEntericTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_etec_flags_empty(self) -> None:
        t = _analyse(_rec("S1", _flags(etec_k99=True)))
        assert t.etec_newly_detected is False
        assert t.etec_cleared is False


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is BovEntericTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is BovEntericTrendDirection.STABLE_CLEAR

    def test_stable_clear_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        text = t.interpretation.lower()
        assert "negative" in text or "no active" in text

    def test_worst_alert_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is BovEntericAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# EMERGING
# ---------------------------------------------------------------------------


class TestEmerging:
    def test_etec_newly_detected(self) -> None:
        recs = [
            _rec("S1", _flags(rota_a=True)),
            _rec("S2", _flags(rota_a=True)),
            _rec("S3", _flags(etec_k99=True, rota_a=True)),
        ]
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.EMERGING
        assert t.etec_newly_detected is True

    def test_etec_newly_detected_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(etec_k99=True)))
        assert t.direction is BovEntericTrendDirection.EMERGING
        assert t.etec_newly_detected is True

    def test_burden_increasing_triggers_emerging(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3", _flags(rota_a=True)),
            _rec("S4", _flags(rota_a=True, bcov=True)),
        ]
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.EMERGING

    def test_emerging_action_mentions_antibiotic_or_monitoring(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(etec_k99=True)))
        text = t.recommended_action.lower()
        assert "antibiotic" in text or "monitoring" in text or "isolate" in text

    def test_etec_emerging_action_mentions_iv(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(etec_k99=True)))
        text = t.recommended_action.lower()
        assert "iv" in text or "intravenous" in text or "antibiotic" in text

    def test_etec_emerging_action_mentions_2h(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(etec_k99=True)))
        assert "2 h" in t.recommended_action or "2h" in t.recommended_action


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_first_positive_last_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        assert t.direction is BovEntericTrendDirection.RESOLVING

    def test_etec_cleared(self) -> None:
        recs = [
            _rec("S1", _flags(etec_k99=True)),
            _rec("S2", _flags(etec_k99=True)),
            _rec("S3"),
        ]
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.RESOLVING
        assert t.etec_cleared is True

    def test_etec_cleared_action_mentions_hygiene(self) -> None:
        recs = [_rec("S1", _flags(etec_k99=True)), _rec("S2")]
        t = analyse_bov_enteric_trend(recs)
        text = t.recommended_action.lower()
        assert "hygiene" in text or "disinfect" in text or "sanit" in text

    def test_burden_decreasing_triggers_resolving(self) -> None:
        recs = [
            _rec("S1", _flags(rota_a=True, bcov=True)),
            _rec("S2", _flags(rota_a=True, bcov=True)),
            _rec("S3"),
            _rec("S4"),
        ]
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.RESOLVING

    def test_resolving_action_mentions_maintain_or_continue(self) -> None:
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "maintain" in text or "continue" in text or "surveillance" in text


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
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.STABLE_ENDEMIC

    def test_stable_endemic_interpretation_non_empty(self) -> None:
        recs = [_rec("S1", _flags(rota_a=True)), _rec("S2", _flags(rota_a=True))]
        t = analyse_bov_enteric_trend(recs)
        assert t.interpretation.strip()


# ---------------------------------------------------------------------------
# Salmonella persistence
# ---------------------------------------------------------------------------


class TestSalmonellaPersistent:
    def test_persistent_when_all_positive(self) -> None:
        recs = [_rec("S1", _flags(salmonella=True)), _rec("S2", _flags(salmonella=True))]
        t = analyse_bov_enteric_trend(recs)
        assert t.salmonella_persistent is True

    def test_not_persistent_when_some_negative(self) -> None:
        recs = [
            _rec("S1", _flags(salmonella=True)),
            _rec("S2"),
            _rec("S3", _flags(salmonella=True)),
        ]
        t = analyse_bov_enteric_trend(recs)
        assert t.salmonella_persistent is False

    def test_persistent_note_in_interpretation(self) -> None:
        recs = [_rec("S1", _flags(salmonella=True)), _rec("S2", _flags(salmonella=True))]
        t = analyse_bov_enteric_trend(recs)
        text = t.interpretation.lower()
        assert "salmonella" in text or "persistent" in text or "herd" in text


# ---------------------------------------------------------------------------
# Zoonotic risk count
# ---------------------------------------------------------------------------


class TestZoonoticRiskCount:
    def test_zero_when_no_zoonotic(self) -> None:
        recs = [_rec("S1", _flags(rota_a=True)), _rec("S2", _flags(bcov=True))]
        t = analyse_bov_enteric_trend(recs)
        assert t.zoonotic_risk_count == 0

    def test_counts_salmonella_intervals(self) -> None:
        recs = [
            _rec("S1", _flags(salmonella=True)),
            _rec("S2"),
            _rec("S3", _flags(salmonella=True)),
        ]
        t = analyse_bov_enteric_trend(recs)
        assert t.zoonotic_risk_count == 2

    def test_counts_crypto_intervals(self) -> None:
        recs = [_rec("S1", _flags(crypto=True)), _rec("S2")]
        t = analyse_bov_enteric_trend(recs)
        assert t.zoonotic_risk_count == 1

    def test_counts_salmonella_and_crypto_same_interval_once(self) -> None:
        recs = [_rec("S1", _flags(salmonella=True, crypto=True)), _rec("S2")]
        t = analyse_bov_enteric_trend(recs)
        assert t.zoonotic_risk_count == 1


# ---------------------------------------------------------------------------
# Antibiotic required count
# ---------------------------------------------------------------------------


class TestAntibioticRequiredCount:
    def test_zero_when_no_bacterial(self) -> None:
        recs = [_rec("S1", _flags(rota_a=True)), _rec("S2", _flags(bcov=True))]
        t = analyse_bov_enteric_trend(recs)
        assert t.antibiotic_required_count == 0

    def test_counts_etec_intervals(self) -> None:
        recs = [_rec("S1", _flags(etec_k99=True)), _rec("S2"), _rec("S3", _flags(etec_k99=True))]
        t = analyse_bov_enteric_trend(recs)
        assert t.antibiotic_required_count == 2

    def test_counts_salmonella_intervals(self) -> None:
        recs = [_rec("S1", _flags(salmonella=True)), _rec("S2")]
        t = analyse_bov_enteric_trend(recs)
        assert t.antibiotic_required_count == 1

    def test_etec_salmonella_same_interval_counted_once(self) -> None:
        recs = [_rec("S1", _flags(etec_k99=True, salmonella=True)), _rec("S2")]
        t = analyse_bov_enteric_trend(recs)
        assert t.antibiotic_required_count == 1


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is BovEntericAlertLevel.NEGATIVE

    def test_etec_gives_critical_worst(self) -> None:
        recs = [_rec("S1"), _rec("S2", _flags(etec_k99=True)), _rec("S3")]
        t = analyse_bov_enteric_trend(recs)
        assert t.worst_alert_level is BovEntericAlertLevel.CRITICAL

    def test_worst_preserved_after_clearance(self) -> None:
        recs = [_rec("S1", _flags(etec_k99=True)), _rec("S2")]
        t = analyse_bov_enteric_trend(recs)
        assert t.worst_alert_level is BovEntericAlertLevel.CRITICAL

    def test_rota_a_only_worst_is_low(self) -> None:
        recs = [_rec("S1", _flags(rota_a=True)), _rec("S2", _flags(rota_a=True))]
        t = analyse_bov_enteric_trend(recs)
        assert t.worst_alert_level is BovEntericAlertLevel.LOW


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
            "etec_newly_detected",
            "etec_cleared",
            "salmonella_persistent",
            "zoonotic_risk_count",
            "antibiotic_required_count",
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
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_etec_newly_detected_is_bool(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(etec_k99=True)))
        d = t.to_dict()
        assert isinstance(d["etec_newly_detected"], bool)
        assert d["etec_newly_detected"] is True

    def test_timeline_flags_match_input(self) -> None:
        recs = [
            _rec("W1", _flags(rota_a=True)),
            _rec("W2", _flags(etec_k99=True, rota_a=True)),
        ]
        t = analyse_bov_enteric_trend(recs)
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert tl[0]["rota_a"] is True
        assert tl[0]["etec_k99"] is False
        assert tl[1]["etec_k99"] is True


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
                    "sample_id": "CG1-wk1",
                    "date": "2026-01-07",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "1",
                },
                {
                    "sample_id": "CG1-wk2",
                    "date": "2026-01-14",
                    "etec_k99": "1",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "1",
                },
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "CG1-wk1"
        assert recs[0].flags.rota_a is True
        assert recs[0].flags.etec_k99 is False
        assert recs[0].date == "2026-01-07"
        assert recs[1].flags.etec_k99 is True

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "etec_k99": "true",
                    "salmonella": "false",
                    "crypto": "TRUE",
                    "bcov": "FALSE",
                    "rota_a": "1",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.etec_k99 is True
        assert recs[0].flags.salmonella is False
        assert recs[0].flags.crypto is True
        assert recs[0].flags.bcov is False
        assert recs[0].flags.rota_a is True

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "0",
                }
            ],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "S1",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                }
            ],
        )
        with pytest.raises(ValueError, match="rota_a"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {
                    "sample_id": "A",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "1",
                },
                {
                    "sample_id": "B",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "1",
                },
                {
                    "sample_id": "C",
                    "etec_k99": "1",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "1",
                },
            ],
        )
        recs = records_from_csv(p)
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.EMERGING
        assert t.etec_newly_detected is True


# ---------------------------------------------------------------------------
# records_from_enteric_json
# ---------------------------------------------------------------------------


class TestRecordsFromEntericJson:
    def _write_enteric_json(
        self,
        path: Path,
        etec_k99: bool = False,
        salmonella: bool = False,
        crypto: bool = False,
        bcov: bool = False,
        rota_a: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "etec_k99": etec_k99,
                "salmonella": salmonella,
                "crypto": crypto,
                "bcov": bcov,
                "rota_a": rota_a,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "week-01.json"
        p2 = tmp_path / "week-02.json"
        self._write_enteric_json(p1, rota_a=True)
        self._write_enteric_json(p2)
        recs = records_from_enteric_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "week-01"
        assert recs[0].flags.rota_a is True
        assert recs[1].flags.rota_a is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_enteric_json(p, date="2026-01-14")
        recs = records_from_enteric_json([p])
        assert recs[0].date == "2026-01-14"

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_enteric_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, (rota, etec) in enumerate([(True, False), (True, False), (True, True)]):
            p = tmp_path / f"s{i}.json"
            self._write_enteric_json(p, rota_a=rota, etec_k99=etec)
            files.append(p)
        recs = records_from_enteric_json(files)
        t = analyse_bov_enteric_trend(recs)
        assert t.direction is BovEntericTrendDirection.EMERGING
        assert t.etec_newly_detected is True


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
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"

    def test_etec_newly_detected_in_json(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(etec_k99=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["etec_newly_detected"] is True

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "trend.json"
        write_trend_json(_analyse(_rec("S1"), _rec("S2")), out)
        assert out.exists()

    def test_ends_with_newline(self, tmp_path: Path) -> None:
        out = tmp_path / "trend.json"
        write_trend_json(_analyse(_rec("S1"), _rec("S2")), out)
        assert out.read_bytes().endswith(b"\n")


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
            "etec_k99",
            "salmonella",
            "crypto",
            "bcov",
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
        t = _analyse(_rec("S1", _flags(etec_k99=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_rota_a_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(rota_a=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["rota_a"] == "1"
        assert dr[1]["rota_a"] == "0"


# ---------------------------------------------------------------------------
# Immutability
# ---------------------------------------------------------------------------


class TestImmutability:
    def test_frozen_trend(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        with pytest.raises((AttributeError, TypeError)):
            t.n_samples = 99  # type: ignore[misc]

    def test_frozen_record(self) -> None:
        r = _rec("S1")
        with pytest.raises((AttributeError, TypeError)):
            r.sample_id = "X"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestBovEntericTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["bov-enteric-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        etec_k99: bool = False,
        salmonella: bool = False,
        rota_a: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "LOW",
            "panel_flags": {
                "etec_k99": etec_k99,
                "salmonella": salmonella,
                "crypto": False,
                "bcov": False,
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

    def test_two_files_exits_0(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, rota_a=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, etec_k99=True)
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
                {
                    "sample_id": "A",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "0",
                },
                {
                    "sample_id": "B",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "1",
                },
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, etec_k99=True)
        out = tmp_path / "trend.json"
        result = self._run(
            [
                "--enteric-result",
                str(p1),
                "--enteric-result",
                str(p2),
                "--out-json",
                str(out),
            ]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["direction"] == "EMERGING"
        assert data["etec_newly_detected"] is True

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, rota_a=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            [
                "--enteric-result",
                str(p1),
                "--enteric-result",
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

    def test_csv_and_enteric_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [
                {
                    "sample_id": "A",
                    "etec_k99": "0",
                    "salmonella": "0",
                    "crypto": "0",
                    "bcov": "0",
                    "rota_a": "0",
                }
            ],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--enteric-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "calf_grp1_jan.json"
        p2 = tmp_path / "calf_grp1_feb.json"
        self._write_json(p1)
        self._write_json(p2, etec_k99=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "calf_grp1_jan" in output
        assert "calf_grp1_feb" in output

    def test_etec_newly_detected_highlighted(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, etec_k99=True)
        result = self._run(["--enteric-result", str(p1), "--enteric-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "ETEC" in output
        assert "NEWLY DETECTED" in output or "newly" in output.lower()
