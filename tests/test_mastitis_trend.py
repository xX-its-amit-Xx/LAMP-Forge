"""Tests for mastitis_trend -- no external binaries required.

Direction logic is deterministic, so tests use exact comparisons on enum values.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lamp_forge.mastitis_risk import MastitisAlertLevel, MastitisPanelFlags
from lamp_forge.mastitis_trend import (
    MastitisRecord,
    MastitisTrend,
    MastitisTrendDirection,
    analyse_mastitis_trend,
    records_from_csv,
    records_from_mastitis_json,
    write_trend_csv,
    write_trend_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    saur: bool = False,
    saga: bool = False,
    sube: bool = False,
    ecoli: bool = False,
) -> MastitisPanelFlags:
    return MastitisPanelFlags(saur=saur, saga=saga, sube=sube, ecoli=ecoli)


def _rec(
    sample_id: str,
    flags: MastitisPanelFlags | None = None,
    date: str | None = None,
) -> MastitisRecord:
    return MastitisRecord(sample_id=sample_id, flags=flags or _flags(), date=date)


def _analyse(*recs: MastitisRecord) -> MastitisTrend:
    return analyse_mastitis_trend(list(recs))


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            analyse_mastitis_trend([])

    def test_single_sample_insufficient_data(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)))
        assert t.direction is MastitisTrendDirection.INSUFFICIENT_DATA
        assert t.n_samples == 1

    def test_single_sample_text_non_empty(self) -> None:
        t = _analyse(_rec("S1"))
        assert t.interpretation.strip()
        assert t.recommended_action.strip()

    def test_single_sample_event_flags_false(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)))
        assert t.contagious_newly_detected is False
        assert t.contagious_cleared is False
        assert t.persistent_contagious_likely is False

    def test_single_positive_worst_level_high(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)))
        assert t.worst_alert_level is MastitisAlertLevel.HIGH

    def test_single_critical_worst_level(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True, saga=True)))
        assert t.worst_alert_level is MastitisAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# STABLE_CLEAR
# ---------------------------------------------------------------------------


class TestStableClear:
    def test_all_negative_two_samples(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.direction is MastitisTrendDirection.STABLE_CLEAR

    def test_all_negative_four_samples(self) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"), _rec("D"))
        assert t.direction is MastitisTrendDirection.STABLE_CLEAR

    def test_stable_clear_text_mentions_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"), _rec("S3"))
        assert "negative" in t.interpretation.lower() or "no mastitis" in t.interpretation.lower()

    def test_worst_alert_negative_when_all_clear(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is MastitisAlertLevel.NEGATIVE

    def test_no_event_flags_when_clear(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.contagious_newly_detected is False
        assert t.contagious_cleared is False
        assert t.persistent_contagious_likely is False


# ---------------------------------------------------------------------------
# STABLE_ENDEMIC
# ---------------------------------------------------------------------------


class TestStableEndemic:
    def test_two_environmental_only_samples(self) -> None:
        t = _analyse(
            _rec("S1", _flags(sube=True)),
            _rec("S2", _flags(sube=True)),
        )
        assert t.direction is MastitisTrendDirection.STABLE_ENDEMIC

    def test_ecoli_recurring(self) -> None:
        t = _analyse(
            _rec("S1", _flags(ecoli=True)),
            _rec("S2", _flags(ecoli=True)),
        )
        assert t.direction is MastitisTrendDirection.STABLE_ENDEMIC

    def test_mixed_env_recurring(self) -> None:
        t = _analyse(
            _rec("S1", _flags(sube=True)),
            _rec("S2", _flags(ecoli=True)),
        )
        assert t.direction is MastitisTrendDirection.STABLE_ENDEMIC

    def test_single_env_not_endemic(self) -> None:
        t = _analyse(
            _rec("S1"),
            _rec("S2", _flags(sube=True)),
        )
        assert t.direction is MastitisTrendDirection.STABLE_CLEAR

    def test_endemic_action_mentions_bedding(self) -> None:
        t = _analyse(
            _rec("S1", _flags(sube=True)),
            _rec("S2", _flags(sube=True)),
        )
        text = t.recommended_action.lower()
        assert "bedding" in text or "cubicle" in text or "husbandry" in text

    def test_no_contagious_flags_in_endemic(self) -> None:
        t = _analyse(
            _rec("S1", _flags(sube=True)),
            _rec("S2", _flags(sube=True)),
        )
        assert t.contagious_newly_detected is False
        assert t.contagious_cleared is False


# ---------------------------------------------------------------------------
# NEWLY_DETECTED
# ---------------------------------------------------------------------------


class TestNewlyDetected:
    def test_saur_first_positive_in_last_sample(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        assert t.direction is MastitisTrendDirection.NEWLY_DETECTED
        assert t.contagious_newly_detected is True

    def test_saga_first_positive(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saga=True)))
        assert t.direction is MastitisTrendDirection.NEWLY_DETECTED

    def test_first_positive_after_three_negatives(self) -> None:
        recs = [
            _rec("S1"),
            _rec("S2"),
            _rec("S3"),
            _rec("S4", _flags(saur=True)),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.direction is MastitisTrendDirection.NEWLY_DETECTED
        assert t.contagious_newly_detected is True

    def test_newly_detected_action_mentions_segregate(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        text = t.recommended_action.lower()
        assert "segregate" in text or "screen" in text

    def test_newly_detected_interpretation_mentions_first(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        text = t.interpretation.lower()
        assert "newly" in text or "first" in text or "absent" in text

    def test_saur_name_in_interpretation(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        assert "s. aureus" in t.interpretation.lower() or "aureus" in t.interpretation.lower()

    def test_not_newly_detected_when_prior_positive(self) -> None:
        t = _analyse(
            _rec("S1", _flags(saur=True)),
            _rec("S2", _flags(saur=True)),
        )
        assert t.contagious_newly_detected is False

    def test_env_only_then_contagious_is_newly_detected(self) -> None:
        t = _analyse(
            _rec("S1", _flags(sube=True)),
            _rec("S2", _flags(saur=True)),
        )
        assert t.direction is MastitisTrendDirection.NEWLY_DETECTED


# ---------------------------------------------------------------------------
# PERSISTENT
# ---------------------------------------------------------------------------


class TestPersistent:
    def test_two_consecutive_saur_positives(self) -> None:
        t = _analyse(
            _rec("S1", _flags(saur=True)),
            _rec("S2", _flags(saur=True)),
        )
        assert t.direction is MastitisTrendDirection.PERSISTENT

    def test_positive_then_negative_then_positive(self) -> None:
        recs = [
            _rec("S1", _flags(saur=True)),
            _rec("S2"),
            _rec("S3", _flags(saur=True)),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.direction is MastitisTrendDirection.PERSISTENT

    def test_persistent_contagious_likely_true_for_high_ratio(self) -> None:
        recs = [
            _rec("S1", _flags(saur=True)),
            _rec("S2", _flags(saur=True)),
            _rec("S3", _flags(saur=True)),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.persistent_contagious_likely is True

    def test_persistent_contagious_likely_false_for_only_two_samples(self) -> None:
        t = _analyse(
            _rec("S1", _flags(saur=True)),
            _rec("S2", _flags(saur=True)),
        )
        assert t.persistent_contagious_likely is False

    def test_persistent_contagious_likely_false_for_low_ratio(self) -> None:
        recs = [
            _rec("S1", _flags(saur=True)),
            _rec("S2"),
            _rec("S3"),
            _rec("S4", _flags(saur=True)),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.direction is MastitisTrendDirection.PERSISTENT
        assert t.persistent_contagious_likely is False

    def test_persistent_action_mentions_screening(self) -> None:
        recs = [_rec("S1", _flags(saur=True)), _rec("S2", _flags(saur=True))]
        t = analyse_mastitis_trend(recs)
        assert (
            "screen" in t.recommended_action.lower() or "segregat" in t.recommended_action.lower()
        )

    def test_saga_persistent(self) -> None:
        t = _analyse(
            _rec("S1", _flags(saga=True)),
            _rec("S2", _flags(saga=True)),
        )
        assert t.direction is MastitisTrendDirection.PERSISTENT


# ---------------------------------------------------------------------------
# RESOLVING
# ---------------------------------------------------------------------------


class TestResolving:
    def test_saur_then_negative(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        assert t.direction is MastitisTrendDirection.RESOLVING
        assert t.contagious_cleared is True

    def test_multiple_positives_then_negative(self) -> None:
        recs = [
            _rec("S1", _flags(saur=True)),
            _rec("S2", _flags(saur=True)),
            _rec("S3"),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.direction is MastitisTrendDirection.RESOLVING
        assert t.contagious_cleared is True

    def test_resolving_action_mentions_consecutive(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        text = t.recommended_action.lower()
        assert "consecutive" in text or "confirm" in text

    def test_cleared_flag_set_when_resolving(self) -> None:
        t = _analyse(_rec("S1", _flags(saga=True)), _rec("S2"))
        assert t.contagious_cleared is True
        assert t.contagious_newly_detected is False

    def test_worst_level_preserved_after_clearance(self) -> None:
        t = _analyse(
            _rec("S1", _flags(saur=True)),
            _rec("S2"),
        )
        assert t.worst_alert_level is MastitisAlertLevel.HIGH

    def test_env_after_contagious_is_resolving(self) -> None:
        t = _analyse(
            _rec("S1", _flags(saur=True)),
            _rec("S2", _flags(sube=True)),
        )
        assert t.direction is MastitisTrendDirection.RESOLVING


# ---------------------------------------------------------------------------
# Worst alert level
# ---------------------------------------------------------------------------


class TestWorstAlertLevel:
    def test_all_negative_worst_is_negative(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert t.worst_alert_level is MastitisAlertLevel.NEGATIVE

    def test_saur_gives_high(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        assert t.worst_alert_level is MastitisAlertLevel.HIGH

    def test_ecoli_alone_gives_low(self) -> None:
        t = _analyse(_rec("S1", _flags(ecoli=True)), _rec("S2"))
        assert t.worst_alert_level is MastitisAlertLevel.LOW

    def test_saur_and_saga_gives_critical(self) -> None:
        recs = [
            _rec("S1", _flags(saur=True, saga=True)),
            _rec("S2"),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.worst_alert_level is MastitisAlertLevel.CRITICAL

    def test_worst_higher_than_most_recent(self) -> None:
        recs = [
            _rec("S1", _flags(saur=True, saga=True)),
            _rec("S2", _flags(ecoli=True)),
        ]
        t = analyse_mastitis_trend(recs)
        assert t.worst_alert_level is MastitisAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        d = t.to_dict()
        for key in (
            "n_samples",
            "direction",
            "contagious_newly_detected",
            "contagious_cleared",
            "persistent_contagious_likely",
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
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        serialised = json.dumps(t.to_dict())
        recovered = json.loads(serialised)
        assert recovered["direction"] == "RESOLVING"

    def test_direction_value_is_string(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2"))
        assert isinstance(t.to_dict()["direction"], str)

    def test_newly_detected_is_bool(self) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        d = t.to_dict()
        assert isinstance(d["contagious_newly_detected"], bool)
        assert d["contagious_newly_detected"] is True

    def test_timeline_has_alert_level_column(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        assert "alert_level" in tl[0]
        assert tl[0]["alert_level"] == "HIGH"

    def test_timeline_has_all_flag_columns(self) -> None:
        t = _analyse(_rec("S1", _flags(saur=True, sube=True)), _rec("S2"))
        tl = t.to_dict()["timeline"]
        assert isinstance(tl, list)
        row = tl[0]
        for key in ("saur", "saga", "sube", "ecoli"):
            assert key in row


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
                {"sample_id": "HerdA-Jan", "date": "2026-01-15", "saur": "1", "saga": "0"},
                {"sample_id": "HerdA-Feb", "date": "2026-02-15", "saur": "0", "saga": "0"},
            ],
        )
        recs = records_from_csv(p)
        assert len(recs) == 2
        assert recs[0].sample_id == "HerdA-Jan"
        assert recs[0].flags.saur is True
        assert recs[0].flags.saga is False
        assert recs[0].date == "2026-01-15"
        assert recs[1].flags.saur is False

    def test_true_false_strings(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "saur": "true", "saga": "false"}],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.saur is True
        assert recs[0].flags.saga is False

    def test_optional_sube_ecoli_default_false(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "saur": "0", "saga": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].flags.sube is False
        assert recs[0].flags.ecoli is False

    def test_optional_date_column_absent(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "saur": "0", "saga": "0"}],
        )
        recs = records_from_csv(p)
        assert recs[0].date is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"sample_id": "S1", "saur": "0"}],
        )
        with pytest.raises(ValueError, match="saga"):
            records_from_csv(p)

    def test_missing_sample_id_raises(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [{"saur": "0", "saga": "0"}],
        )
        with pytest.raises(ValueError, match="sample_id"):
            records_from_csv(p)

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        p = self._write_csv(
            tmp_path,
            [
                {"sample_id": "A", "saur": "0", "saga": "0"},
                {"sample_id": "B", "saur": "0", "saga": "0"},
                {"sample_id": "C", "saur": "1", "saga": "0"},
            ],
        )
        recs = records_from_csv(p)
        t = analyse_mastitis_trend(recs)
        assert t.direction is MastitisTrendDirection.NEWLY_DETECTED
        assert t.contagious_newly_detected is True


# ---------------------------------------------------------------------------
# records_from_mastitis_json
# ---------------------------------------------------------------------------


class TestRecordsFromMastitisJson:
    def _write_mastitis_json(
        self,
        path: Path,
        saur: bool = False,
        saga: bool = False,
        sube: bool = False,
        ecoli: bool = False,
        date: str | None = None,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "NEGATIVE",
            "panel_flags": {
                "saur": saur,
                "saga": saga,
                "sube": sube,
                "ecoli": ecoli,
            },
        }
        if date is not None:
            data["date"] = date
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def test_basic_parse(self, tmp_path: Path) -> None:
        p1 = tmp_path / "2026-01.json"
        p2 = tmp_path / "2026-02.json"
        self._write_mastitis_json(p1, saur=True)
        self._write_mastitis_json(p2)
        recs = records_from_mastitis_json([p1, p2])
        assert len(recs) == 2
        assert recs[0].sample_id == "2026-01"
        assert recs[0].flags.saur is True
        assert recs[1].flags.saur is False

    def test_date_field_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "sample.json"
        self._write_mastitis_json(p, date="2026-03-15")
        recs = records_from_mastitis_json([p])
        assert recs[0].date == "2026-03-15"

    def test_all_flags_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "full.json"
        self._write_mastitis_json(p, saur=True, saga=False, sube=True, ecoli=False)
        recs = records_from_mastitis_json([p])
        assert recs[0].flags.saur is True
        assert recs[0].flags.saga is False
        assert recs[0].flags.sube is True
        assert recs[0].flags.ecoli is False

    def test_missing_panel_flags_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"alert_level": "HIGH"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_mastitis_json([p])

    def test_panel_flags_not_dict_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.json"
        p.write_text('{"panel_flags": "bad"}')
        with pytest.raises(ValueError, match="panel_flags"):
            records_from_mastitis_json([p])

    def test_round_trip_through_analyse(self, tmp_path: Path) -> None:
        files = []
        for i, saur in enumerate([False, True, True]):
            p = tmp_path / f"s{i}.json"
            self._write_mastitis_json(p, saur=saur)
            files.append(p)
        recs = records_from_mastitis_json(files)
        t = analyse_mastitis_trend(recs)
        assert t.direction is MastitisTrendDirection.PERSISTENT


# ---------------------------------------------------------------------------
# write_trend_json
# ---------------------------------------------------------------------------


class TestWriteTrendJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1"), _rec("S2", _flags(saur=True)))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        out = tmp_path / "trend.json"
        write_trend_json(t, out)
        data = json.loads(out.read_text())
        assert data["direction"] == "RESOLVING"
        assert data["contagious_cleared"] is True

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
            "saur",
            "saga",
            "sube",
            "ecoli",
            "alert_level",
        ]

    def test_row_count_equals_n_samples(self, tmp_path: Path) -> None:
        t = _analyse(_rec("A"), _rec("B"), _rec("C"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 4  # header + 3 data rows

    def test_alert_level_column_filled(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        valid_levels = {"NEGATIVE", "LOW", "MODERATE", "HIGH", "CRITICAL"}
        assert dr[0]["alert_level"] in valid_levels

    def test_saur_column_values(self, tmp_path: Path) -> None:
        t = _analyse(_rec("S1", _flags(saur=True)), _rec("S2"))
        out = tmp_path / "trend.csv"
        write_trend_csv(t, out)
        dr = list(csv.DictReader(out.open()))
        assert dr[0]["saur"] == "1"
        assert dr[1]["saur"] == "0"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestMastitisTrendCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["mastitis-trend", *args], obj={})

    def _write_json(
        self,
        path: Path,
        saur: bool = False,
        saga: bool = False,
    ) -> None:
        data: dict[str, object] = {
            "alert_level": "NEGATIVE",
            "panel_flags": {
                "saur": saur,
                "saga": saga,
                "sube": False,
                "ecoli": False,
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
        self._write_json(p2, saur=True)
        result = self._run(["--mastitis-result", str(p1), "--mastitis-result", str(p2)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_direction(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, saur=True)
        result = self._run(["--mastitis-result", str(p1), "--mastitis-result", str(p2)])
        assert "NEWLY_DETECTED" in result.output  # type: ignore[union-attr]

    def test_resolving_shown_for_clearing(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, saur=True)
        self._write_json(p2)
        result = self._run(["--mastitis-result", str(p1), "--mastitis-result", str(p2)])
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_stable_clear_for_all_negative(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2)
        result = self._run(["--mastitis-result", str(p1), "--mastitis-result", str(p2)])
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_insufficient_data_for_one_file(self, tmp_path: Path) -> None:
        p = tmp_path / "s1.json"
        self._write_json(p, saur=True)
        result = self._run(["--mastitis-result", str(p)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT" in result.output  # type: ignore[union-attr]

    def test_csv_input_exits_0(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "monitoring.csv"
        self._write_csv(
            csv_path,
            [
                {"sample_id": "A", "saur": "0", "saga": "0"},
                {"sample_id": "B", "saur": "1", "saga": "0"},
            ],
        )
        result = self._run(["--csv", str(csv_path)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1)
        self._write_json(p2, saur=True)
        out = tmp_path / "trend.json"
        result = self._run(
            [
                "--mastitis-result",
                str(p1),
                "--mastitis-result",
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
        self._write_json(p1, saur=True)
        self._write_json(p2)
        out = tmp_path / "trend.csv"
        result = self._run(
            [
                "--mastitis-result",
                str(p1),
                "--mastitis-result",
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

    def test_csv_and_mastitis_result_mutual_exclusion(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "m.csv"
        self._write_csv(
            csv_path,
            [{"sample_id": "A", "saur": "0", "saga": "0"}],
        )
        p = tmp_path / "s.json"
        self._write_json(p)
        result = self._run(["--csv", str(csv_path), "--mastitis-result", str(p)])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_output_table_shows_sample_ids(self, tmp_path: Path) -> None:
        p1 = tmp_path / "herd_a_jan.json"
        p2 = tmp_path / "herd_a_feb.json"
        self._write_json(p1)
        self._write_json(p2, saur=True)
        result = self._run(["--mastitis-result", str(p1), "--mastitis-result", str(p2)])
        output = result.output  # type: ignore[union-attr]
        assert "herd_a_jan" in output
        assert "herd_a_feb" in output

    def test_persistent_shown_for_consecutive_positives(self, tmp_path: Path) -> None:
        p1, p2 = tmp_path / "s1.json", tmp_path / "s2.json"
        self._write_json(p1, saur=True)
        self._write_json(p2, saur=True)
        result = self._run(["--mastitis-result", str(p1), "--mastitis-result", str(p2)])
        assert "PERSISTENT" in result.output  # type: ignore[union-attr]
