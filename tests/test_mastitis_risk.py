"""Tests for the bovine mastitis LAMP panel risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical alert levels.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

from lamp_forge.mastitis_risk import (
    MastitisAlertLevel,
    MastitisPanelFlags,
    MastitisRiskAssessment,
    assess_mastitis_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> MastitisPanelFlags:
    defaults: dict[str, bool] = {
        "saur": False,
        "saga": False,
        "sube": False,
        "ecoli": False,
    }
    defaults.update(kwargs)
    return MastitisPanelFlags(**defaults)


def _assess(**kwargs: bool) -> MastitisRiskAssessment:
    return assess_mastitis_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is MastitisAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_saur_alone_is_high(self) -> None:
        a = _assess(saur=True)
        assert a.alert_level is MastitisAlertLevel.HIGH
        assert a.alert_score == 40

    def test_saga_alone_is_high(self) -> None:
        a = _assess(saga=True)
        assert a.alert_level is MastitisAlertLevel.HIGH
        assert a.alert_score == 35

    def test_sube_alone_is_moderate(self) -> None:
        a = _assess(sube=True)
        assert a.alert_level is MastitisAlertLevel.MODERATE
        assert a.alert_score == 20

    def test_ecoli_alone_is_low(self) -> None:
        a = _assess(ecoli=True)
        assert a.alert_level is MastitisAlertLevel.LOW
        assert a.alert_score == 15

    def test_score_increases_with_targets(self) -> None:
        s1 = _assess(ecoli=True).alert_score
        s2 = _assess(ecoli=True, sube=True).alert_score
        s3 = _assess(ecoli=True, sube=True, saga=True).alert_score
        assert s1 < s2 < s3

    def test_score_capped_at_100(self) -> None:
        a = _assess(saur=True, saga=True, sube=True, ecoli=True)
        assert a.alert_score <= 100

    def test_score_non_negative(self) -> None:
        for kwargs in [{}, {"saur": True}, {"ecoli": True}]:
            assert _assess(**kwargs).alert_score >= 0


# ---------------------------------------------------------------------------
# Multi-target alert levels
# ---------------------------------------------------------------------------


class TestMultiTargetAlertLevels:
    def test_saur_saga_is_critical(self) -> None:
        # 40 + 35 = 75 -> CRITICAL (>= 70)
        a = _assess(saur=True, saga=True)
        assert a.alert_level is MastitisAlertLevel.CRITICAL
        assert a.alert_score == 75

    def test_saur_sube_is_high(self) -> None:
        # 40 + 20 = 60 -> HIGH (>= 35 but < 70)
        a = _assess(saur=True, sube=True)
        assert a.alert_level is MastitisAlertLevel.HIGH
        assert a.alert_score == 60

    def test_saur_ecoli_is_high(self) -> None:
        # 40 + 15 = 55 -> HIGH
        a = _assess(saur=True, ecoli=True)
        assert a.alert_level is MastitisAlertLevel.HIGH
        assert a.alert_score == 55

    def test_sube_ecoli_is_high(self) -> None:
        # 20 + 15 = 35 -> HIGH (boundary)
        a = _assess(sube=True, ecoli=True)
        assert a.alert_level is MastitisAlertLevel.HIGH
        assert a.alert_score == 35

    def test_all_four_is_critical(self) -> None:
        a = _assess(saur=True, saga=True, sube=True, ecoli=True)
        assert a.alert_level is MastitisAlertLevel.CRITICAL

    def test_saga_sube_is_high(self) -> None:
        # 35 + 20 = 55 -> HIGH
        a = _assess(saga=True, sube=True)
        assert a.alert_level is MastitisAlertLevel.HIGH
        assert a.alert_score == 55


# ---------------------------------------------------------------------------
# Contagious pathogen flags
# ---------------------------------------------------------------------------


class TestContagiousFlags:
    def test_saur_sets_contagious_detected(self) -> None:
        a = _assess(saur=True)
        assert a.contagious_pathogen_detected is True

    def test_saga_sets_contagious_detected(self) -> None:
        a = _assess(saga=True)
        assert a.contagious_pathogen_detected is True

    def test_sube_alone_not_contagious(self) -> None:
        a = _assess(sube=True)
        assert a.contagious_pathogen_detected is False

    def test_ecoli_alone_not_contagious(self) -> None:
        a = _assess(ecoli=True)
        assert a.contagious_pathogen_detected is False

    def test_all_negative_not_contagious(self) -> None:
        a = _assess()
        assert a.contagious_pathogen_detected is False

    def test_saur_saga_is_contagious_coinfection(self) -> None:
        a = _assess(saur=True, saga=True)
        assert a.contagious_coinfection is True

    def test_saur_alone_not_coinfection(self) -> None:
        a = _assess(saur=True)
        assert a.contagious_coinfection is False

    def test_saga_alone_not_coinfection(self) -> None:
        a = _assess(saga=True)
        assert a.contagious_coinfection is False

    def test_sube_ecoli_not_coinfection(self) -> None:
        a = _assess(sube=True, ecoli=True)
        assert a.contagious_coinfection is False


# ---------------------------------------------------------------------------
# S. aureus biofilm risk flag
# ---------------------------------------------------------------------------


class TestSaureusFlag:
    def test_saur_positive_sets_biofilm_risk(self) -> None:
        a = _assess(saur=True)
        assert a.s_aureus_biofilm_risk is True

    def test_saga_alone_no_biofilm_risk(self) -> None:
        a = _assess(saga=True)
        assert a.s_aureus_biofilm_risk is False

    def test_all_negative_no_biofilm_risk(self) -> None:
        a = _assess()
        assert a.s_aureus_biofilm_risk is False

    def test_saur_saga_still_sets_biofilm_risk(self) -> None:
        a = _assess(saur=True, saga=True)
        assert a.s_aureus_biofilm_risk is True


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    def test_empty_when_all_negative(self) -> None:
        a = _assess()
        assert a.active_targets == ()

    def test_canonical_order_saur_saga_sube_ecoli(self) -> None:
        a = _assess(ecoli=True, saur=True, sube=True)
        assert a.active_targets == ("SAUR", "SUBE", "ECOLI")

    def test_all_four_listed_when_all_positive(self) -> None:
        a = _assess(saur=True, saga=True, sube=True, ecoli=True)
        assert len(a.active_targets) == 4
        assert a.active_targets == ("SAUR", "SAGA", "SUBE", "ECOLI")

    def test_only_positive_targets_listed(self) -> None:
        a = _assess(saga=True, ecoli=True)
        assert "SAGA" in a.active_targets
        assert "ECOLI" in a.active_targets
        assert "SAUR" not in a.active_targets
        assert "SUBE" not in a.active_targets


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestTextContent:
    def test_all_negative_mentions_negative(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_saur_mentions_biofilm_or_treatment(self) -> None:
        a = _assess(saur=True)
        text = (a.interpretation + a.recommended_action).lower()
        assert "biofilm" in text or "treatment" in text or "culling" in text

    def test_saur_saga_mentions_herd_screen(self) -> None:
        a = _assess(saur=True, saga=True)
        text = (a.interpretation + a.recommended_action).lower()
        assert "herd" in text

    def test_saga_mentions_penicillin(self) -> None:
        a = _assess(saga=True)
        text = (a.interpretation + a.recommended_action).lower()
        assert "penicillin" in text

    def test_sube_mentions_bedding_or_environment(self) -> None:
        a = _assess(sube=True)
        text = (a.interpretation + a.recommended_action).lower()
        assert "bedding" in text or "environment" in text or "hygiene" in text

    def test_ecoli_mentions_nsaids_or_supportive(self) -> None:
        a = _assess(ecoli=True)
        text = (a.interpretation + a.recommended_action).lower()
        assert "nsaid" in text or "supportive" in text or "endotox" in text

    def test_all_combos_have_non_empty_text(self) -> None:
        keys = ["saur", "saga", "sube", "ecoli"]
        for bits in product([False, True], repeat=4):
            flags = MastitisPanelFlags(**dict(zip(keys, bits, strict=True)))
            a = assess_mastitis_risk(flags)
            assert a.interpretation.strip(), f"Empty interpretation for {flags}"
            assert a.recommended_action.strip(), f"Empty action for {flags}"

    def test_recommended_action_non_empty(self) -> None:
        for a in [_assess(), _assess(saur=True), _assess(ecoli=True), _assess(saga=True)]:
            assert len(a.recommended_action) > 20

    def test_ecoli_action_does_not_mandate_antibiotics(self) -> None:
        a = _assess(ecoli=True)
        text = a.recommended_action.lower()
        # E. coli mastitis is often self-limiting; the action should not
        # unconditionally prescribe antimicrobials
        assert "strip" in text or "nsaid" in text or "supportive" in text


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        d = _assess(saur=True).to_dict()
        for key in (
            "alert_level",
            "alert_score",
            "active_targets",
            "contagious_pathogen_detected",
            "contagious_coinfection",
            "s_aureus_biofilm_risk",
            "interpretation",
            "recommended_action",
            "panel_flags",
        ):
            assert key in d, f"Missing key: {key}"

    def test_alert_level_is_string(self) -> None:
        d = _assess().to_dict()
        assert isinstance(d["alert_level"], str)

    def test_panel_flags_inlined(self) -> None:
        d = _assess(saur=True, ecoli=True).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf["saur"] is True  # type: ignore[index]
        assert pf["ecoli"] is True  # type: ignore[index]
        assert pf["saga"] is False  # type: ignore[index]
        assert pf["sube"] is False  # type: ignore[index]

    def test_json_serialisable_critical(self) -> None:
        d = _assess(saur=True, saga=True).to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "CRITICAL"
        assert recovered["contagious_coinfection"] is True

    def test_json_serialisable_negative(self) -> None:
        d = _assess().to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "NEGATIVE"
        assert recovered["contagious_pathogen_detected"] is False


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_true(self) -> None:
        flags = flags_from_dict({"saur": True, "saga": True, "sube": True, "ecoli": True})
        assert flags == MastitisPanelFlags(saur=True, saga=True, sube=True, ecoli=True)

    def test_missing_keys_default_false(self) -> None:
        flags = flags_from_dict({"saur": True})
        assert flags.saur is True
        assert flags.saga is False
        assert flags.sube is False
        assert flags.ecoli is False

    def test_empty_dict_all_false(self) -> None:
        flags = flags_from_dict({})
        assert flags == MastitisPanelFlags(saur=False, saga=False, sube=False, ecoli=False)

    def test_unknown_keys_ignored(self) -> None:
        flags = flags_from_dict({"saur": True, "extra_key": "ignored"})
        assert flags.saur is True
        assert flags.saga is False


# ---------------------------------------------------------------------------
# write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "mastitis_assessment.json"
        write_assessment_json(_assess(saur=True, saga=True), out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        out = tmp_path / "mastitis_assessment.json"
        write_assessment_json(_assess(saur=True, saga=True), out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert "SAUR" in data["active_targets"]
        assert data["contagious_coinfection"] is True

    def test_negative_result_serialised(self, tmp_path: Path) -> None:
        out = tmp_path / "neg.json"
        write_assessment_json(_assess(), out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "NEGATIVE"
        assert data["active_targets"] == []

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "result.json"
        write_assessment_json(_assess(), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "mastitis_assessment.csv"
        write_assessment_csv(_assess(saur=True), out)
        assert out.exists()

    def test_has_header_row(self, tmp_path: Path) -> None:
        out = tmp_path / "mastitis_assessment.csv"
        write_assessment_csv(_assess(), out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["field", "value"]

    def test_alert_level_row_present(self, tmp_path: Path) -> None:
        out = tmp_path / "mastitis_assessment.csv"
        write_assessment_csv(_assess(saur=True, saga=True), out)
        rows = list(csv.DictReader(out.open()))
        level_rows = [r for r in rows if r["field"] == "alert_level"]
        assert len(level_rows) == 1
        assert level_rows[0]["value"] == "CRITICAL"

    def test_panel_flags_present_as_separate_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "mastitis_assessment.csv"
        write_assessment_csv(_assess(saur=True, ecoli=True), out)
        rows = list(csv.DictReader(out.open()))
        fields = {r["field"] for r in rows}
        assert "target_saur" in fields
        assert "target_ecoli" in fields
        assert "target_saga" in fields
        assert "target_sube" in fields


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestMastitisRiskCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["mastitis-risk", *args], obj={})

    def test_all_negative_exits_0(self) -> None:
        result = self._run([])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_target_table(self) -> None:
        result = self._run(["--saur", "--saga"])
        output = result.output  # type: ignore[union-attr]
        assert "SAUR" in output
        assert "SAGA" in output
        assert "[+]" in output

    def test_saur_saga_shows_critical(self) -> None:
        result = self._run(["--saur", "--saga"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_saur_shows_high(self) -> None:
        result = self._run(["--saur"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_saga_shows_high(self) -> None:
        result = self._run(["--saga"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_sube_shows_moderate(self) -> None:
        result = self._run(["--sube"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]

    def test_ecoli_shows_low(self) -> None:
        result = self._run(["--ecoli"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "LOW" in result.output  # type: ignore[union-attr]

    def test_all_negative_shows_negative(self) -> None:
        result = self._run([])
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_saur_warns_biofilm(self) -> None:
        result = self._run(["--saur"])
        output = result.output.lower()  # type: ignore[union-attr]
        assert "biofilm" in output

    def test_saur_saga_warns_contagious_coinfection(self) -> None:
        result = self._run(["--saur", "--saga"])
        output = result.output.lower()  # type: ignore[union-attr]
        assert "contagious" in output

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        result = self._run(["--saur", "--saga", "--out-json", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert data["contagious_coinfection"] is True

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        result = self._run(["--sube", "--out-csv", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_flag(self, tmp_path: Path) -> None:
        flags_file = tmp_path / "flags.json"
        flags_file.write_text('{"saur": false, "saga": true, "sube": false, "ecoli": false}')
        result = self._run(["--input-json", str(flags_file)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_saga_interpretation_mentions_penicillin(self) -> None:
        result = self._run(["--saga"])
        output = result.output.lower()  # type: ignore[union-attr]
        assert "penicillin" in output

    def test_ecoli_action_does_not_only_mention_antibiotics(self) -> None:
        result = self._run(["--ecoli"])
        output = result.output.lower()  # type: ignore[union-attr]
        assert "strip" in output or "nsaid" in output or "supportive" in output
