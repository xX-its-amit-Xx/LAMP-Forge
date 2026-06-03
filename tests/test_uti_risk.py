"""Tests for the POC UTI LAMP panel risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical alert levels.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

import pytest

from lamp_forge.uti_risk import (
    UtiAlertLevel,
    UtiPanelFlags,
    UtiRiskAssessment,
    assess_uti_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> UtiPanelFlags:
    defaults: dict[str, bool] = {
        "ecoli": False,
        "kpneu": False,
        "efaec": False,
        "ssaph": False,
    }
    defaults.update(kwargs)
    return UtiPanelFlags(**defaults)


def _assess(**kwargs: bool) -> UtiRiskAssessment:
    return assess_uti_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """Each single-pathogen detection must map to the correct alert level and score."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is UtiAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_ecoli_alone_is_high(self) -> None:
        a = _assess(ecoli=True)
        assert a.alert_level is UtiAlertLevel.HIGH
        assert a.alert_score == 40

    def test_kpneu_alone_is_high(self) -> None:
        a = _assess(kpneu=True)
        assert a.alert_level is UtiAlertLevel.HIGH
        assert a.alert_score == 35

    def test_efaec_alone_is_moderate(self) -> None:
        a = _assess(efaec=True)
        assert a.alert_level is UtiAlertLevel.MODERATE
        assert a.alert_score == 25

    def test_ssaph_alone_is_low(self) -> None:
        a = _assess(ssaph=True)
        assert a.alert_level is UtiAlertLevel.LOW
        assert a.alert_score == 15


# ---------------------------------------------------------------------------
# Co-detection scoring
# ---------------------------------------------------------------------------


class TestCoDetectionScores:
    """Combined flags must sum weights correctly (capped at 100)."""

    def test_ecoli_kpneu_score(self) -> None:
        a = _assess(ecoli=True, kpneu=True)
        assert a.alert_score == 75
        assert a.alert_level is UtiAlertLevel.HIGH

    def test_ecoli_efaec_score(self) -> None:
        a = _assess(ecoli=True, efaec=True)
        assert a.alert_score == 65
        assert a.alert_level is UtiAlertLevel.HIGH

    def test_ecoli_ssaph_score(self) -> None:
        a = _assess(ecoli=True, ssaph=True)
        assert a.alert_score == 55
        assert a.alert_level is UtiAlertLevel.HIGH

    def test_kpneu_efaec_score(self) -> None:
        a = _assess(kpneu=True, efaec=True)
        assert a.alert_score == 60
        assert a.alert_level is UtiAlertLevel.HIGH

    def test_efaec_ssaph_score(self) -> None:
        a = _assess(efaec=True, ssaph=True)
        assert a.alert_score == 40
        assert a.alert_level is UtiAlertLevel.HIGH

    def test_all_positive_score_capped(self) -> None:
        a = _assess(ecoli=True, kpneu=True, efaec=True, ssaph=True)
        assert a.alert_score == 100
        assert a.alert_level is UtiAlertLevel.HIGH

    def test_score_never_exceeds_100(self) -> None:
        for ecoli, kpneu, efaec, ssaph in product([False, True], repeat=4):
            a = _assess(ecoli=ecoli, kpneu=kpneu, efaec=efaec, ssaph=ssaph)
            assert 0 <= a.alert_score <= 100

    def test_score_never_below_zero(self) -> None:
        a = _assess()
        assert a.alert_score >= 0


# ---------------------------------------------------------------------------
# Clinical decision flags
# ---------------------------------------------------------------------------


class TestClinicalFlags:
    """gram_negative_detected, extended_spectrum_resistance_risk, antibiotic_guidance_required."""

    def test_all_negative_clears_all_flags(self) -> None:
        a = _assess()
        assert a.gram_negative_detected is False
        assert a.extended_spectrum_resistance_risk is False
        assert a.antibiotic_guidance_required is False

    def test_ecoli_sets_gram_negative(self) -> None:
        a = _assess(ecoli=True)
        assert a.gram_negative_detected is True
        assert a.extended_spectrum_resistance_risk is False
        assert a.antibiotic_guidance_required is True

    def test_kpneu_sets_both_gram_neg_flags(self) -> None:
        a = _assess(kpneu=True)
        assert a.gram_negative_detected is True
        assert a.extended_spectrum_resistance_risk is True
        assert a.antibiotic_guidance_required is True

    def test_efaec_no_gram_neg_flag(self) -> None:
        a = _assess(efaec=True)
        assert a.gram_negative_detected is False
        assert a.extended_spectrum_resistance_risk is False
        assert a.antibiotic_guidance_required is True

    def test_ssaph_no_gram_neg_flag(self) -> None:
        a = _assess(ssaph=True)
        assert a.gram_negative_detected is False
        assert a.extended_spectrum_resistance_risk is False
        assert a.antibiotic_guidance_required is True

    def test_ecoli_kpneu_co_detection_sets_resistance_risk(self) -> None:
        a = _assess(ecoli=True, kpneu=True)
        assert a.gram_negative_detected is True
        assert a.extended_spectrum_resistance_risk is True

    def test_ssaph_efaec_no_gram_neg(self) -> None:
        a = _assess(ssaph=True, efaec=True)
        assert a.gram_negative_detected is False
        assert a.extended_spectrum_resistance_risk is False

    def test_all_positive_sets_all_flags(self) -> None:
        a = _assess(ecoli=True, kpneu=True, efaec=True, ssaph=True)
        assert a.gram_negative_detected is True
        assert a.extended_spectrum_resistance_risk is True
        assert a.antibiotic_guidance_required is True


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    """active_targets must list positive detections in canonical order."""

    def test_all_negative_empty(self) -> None:
        assert _assess().active_targets == ()

    def test_ecoli_only(self) -> None:
        assert _assess(ecoli=True).active_targets == ("E.coli",)

    def test_kpneu_only(self) -> None:
        assert _assess(kpneu=True).active_targets == ("K.pneu",)

    def test_efaec_only(self) -> None:
        assert _assess(efaec=True).active_targets == ("E.faec",)

    def test_ssaph_only(self) -> None:
        assert _assess(ssaph=True).active_targets == ("S.saph",)

    def test_all_positive_canonical_order(self) -> None:
        a = _assess(ecoli=True, kpneu=True, efaec=True, ssaph=True)
        assert a.active_targets == ("E.coli", "K.pneu", "E.faec", "S.saph")

    def test_kpneu_efaec_order(self) -> None:
        a = _assess(kpneu=True, efaec=True)
        assert a.active_targets == ("K.pneu", "E.faec")

    def test_ecoli_ssaph_order(self) -> None:
        a = _assess(ecoli=True, ssaph=True)
        assert a.active_targets == ("E.coli", "S.saph")


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestInterpretationText:
    """Key phrases must appear in interpretation/action strings."""

    def test_negative_interpretation_mentions_negative(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_ecoli_interpretation_mentions_ecoli(self) -> None:
        a = _assess(ecoli=True)
        assert "coli" in a.interpretation.lower()

    def test_ecoli_action_mentions_nitrofurantoin(self) -> None:
        a = _assess(ecoli=True)
        assert "nitrofurantoin" in a.recommended_action.lower()

    def test_kpneu_interpretation_mentions_klebsiella(self) -> None:
        a = _assess(kpneu=True)
        assert "klebsiella" in a.interpretation.lower()

    def test_kpneu_action_mentions_culture(self) -> None:
        a = _assess(kpneu=True)
        assert "culture" in a.recommended_action.lower()

    def test_efaec_interpretation_mentions_enterococcus(self) -> None:
        a = _assess(efaec=True)
        assert "enterococcus" in a.interpretation.lower()

    def test_efaec_action_mentions_ampicillin(self) -> None:
        a = _assess(efaec=True)
        assert "ampicillin" in a.recommended_action.lower()

    def test_ssaph_interpretation_mentions_saprophyticus(self) -> None:
        a = _assess(ssaph=True)
        assert "saprophyticus" in a.interpretation.lower()

    def test_ecoli_kpneu_co_interpretation_mentions_esbl(self) -> None:
        a = _assess(ecoli=True, kpneu=True)
        assert "esbl" in a.interpretation.lower() or "esbl" in a.recommended_action.lower()

    def test_all_non_empty_strings(self) -> None:
        for ecoli, kpneu, efaec, ssaph in product([False, True], repeat=4):
            a = _assess(ecoli=ecoli, kpneu=kpneu, efaec=efaec, ssaph=ssaph)
            assert a.interpretation
            assert a.recommended_action


# ---------------------------------------------------------------------------
# Serialisation: to_dict
# ---------------------------------------------------------------------------


class TestToDict:
    """to_dict must include all required keys with correct types."""

    def test_required_keys_present(self) -> None:
        d = _assess(ecoli=True).to_dict()
        required = {
            "alert_level",
            "alert_score",
            "active_targets",
            "interpretation",
            "recommended_action",
            "gram_negative_detected",
            "extended_spectrum_resistance_risk",
            "antibiotic_guidance_required",
            "panel_flags",
        }
        assert required <= d.keys()

    def test_alert_level_is_string(self) -> None:
        d = _assess(ecoli=True).to_dict()
        assert isinstance(d["alert_level"], str)
        assert d["alert_level"] == "HIGH"

    def test_panel_flags_keys(self) -> None:
        d = _assess(ecoli=True, kpneu=False, efaec=True, ssaph=False).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf == {"ecoli": True, "kpneu": False, "efaec": True, "ssaph": False}

    def test_active_targets_is_list(self) -> None:
        d = _assess(ecoli=True, kpneu=True).to_dict()
        assert isinstance(d["active_targets"], list)
        assert "E.coli" in d["active_targets"]

    def test_negative_alert_level_string(self) -> None:
        d = _assess().to_dict()
        assert d["alert_level"] == "NEGATIVE"

    def test_bool_fields_are_bools(self) -> None:
        d = _assess(kpneu=True).to_dict()
        assert isinstance(d["gram_negative_detected"], bool)
        assert isinstance(d["extended_spectrum_resistance_risk"], bool)
        assert isinstance(d["antibiotic_guidance_required"], bool)


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    """flags_from_dict must parse dicts correctly and tolerate missing keys."""

    def test_all_true(self) -> None:
        f = flags_from_dict({"ecoli": True, "kpneu": True, "efaec": True, "ssaph": True})
        assert f == UtiPanelFlags(ecoli=True, kpneu=True, efaec=True, ssaph=True)

    def test_all_false_explicit(self) -> None:
        f = flags_from_dict({"ecoli": False, "kpneu": False, "efaec": False, "ssaph": False})
        assert f == UtiPanelFlags(ecoli=False, kpneu=False, efaec=False, ssaph=False)

    def test_empty_dict_defaults_to_false(self) -> None:
        f = flags_from_dict({})
        assert f == UtiPanelFlags(ecoli=False, kpneu=False, efaec=False, ssaph=False)

    def test_partial_dict(self) -> None:
        f = flags_from_dict({"ecoli": True})
        assert f == UtiPanelFlags(ecoli=True, kpneu=False, efaec=False, ssaph=False)

    def test_unknown_keys_ignored(self) -> None:
        f = flags_from_dict({"ecoli": True, "unknown_key": True})
        assert f == UtiPanelFlags(ecoli=True, kpneu=False, efaec=False, ssaph=False)

    def test_truthy_int_converted(self) -> None:
        f = flags_from_dict({"ecoli": 1, "kpneu": 0, "efaec": 1, "ssaph": 0})
        assert f.ecoli is True
        assert f.kpneu is False
        assert f.efaec is True
        assert f.ssaph is False


# ---------------------------------------------------------------------------
# File output: write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    """write_assessment_json must produce valid, parseable JSON."""

    def test_json_round_trip(self, tmp_path: Path) -> None:
        assessment = _assess(ecoli=True, efaec=True)
        out = tmp_path / "out.json"
        write_assessment_json(assessment, out)
        assert out.exists()
        data = json.loads(out.read_text(encoding="utf-8"))
        assert data["alert_level"] == "HIGH"
        assert data["panel_flags"]["ecoli"] is True
        assert data["panel_flags"]["efaec"] is True

    def test_json_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(_assess(), out)
        assert out.exists()

    def test_json_ends_with_newline(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        write_assessment_json(_assess(kpneu=True), out)
        content = out.read_bytes()
        assert content.endswith(b"\n")

    def test_json_negative_result(self, tmp_path: Path) -> None:
        out = tmp_path / "negative.json"
        write_assessment_json(_assess(), out)
        data = json.loads(out.read_text(encoding="utf-8"))
        assert data["alert_level"] == "NEGATIVE"
        assert data["alert_score"] == 0


# ---------------------------------------------------------------------------
# File output: write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    """write_assessment_csv must produce a valid two-column key-value CSV."""

    def test_csv_has_header(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(ssaph=True), out)
        with out.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert header == ["field", "value"]

    def test_csv_contains_alert_level(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(efaec=True), out)
        content = out.read_text(encoding="utf-8")
        assert "MODERATE" in content

    def test_csv_target_flags_positive_negative(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(ecoli=True, kpneu=False, efaec=False, ssaph=False), out)
        content = out.read_text(encoding="utf-8")
        assert "positive" in content
        assert "negative" in content

    def test_csv_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "deep" / "nested" / "out.csv"
        write_assessment_csv(_assess(), out)
        assert out.exists()

    def test_csv_high_alert_for_ecoli(self, tmp_path: Path) -> None:
        out = tmp_path / "ecoli.csv"
        write_assessment_csv(_assess(ecoli=True), out)
        content = out.read_text(encoding="utf-8")
        assert "HIGH" in content


# ---------------------------------------------------------------------------
# Return type and immutability
# ---------------------------------------------------------------------------


class TestReturnTypeAndImmutability:
    """assess_uti_risk must return a frozen UtiRiskAssessment."""

    def test_returns_correct_type(self) -> None:
        a = _assess(ecoli=True)
        assert isinstance(a, UtiRiskAssessment)

    def test_frozen_dataclass(self) -> None:
        a = _assess(kpneu=True)
        with pytest.raises((AttributeError, TypeError)):
            a.alert_score = 0  # type: ignore[misc]

    def test_flags_field_matches_input(self) -> None:
        flags = _flags(ecoli=False, kpneu=True, efaec=False, ssaph=False)
        a = assess_uti_risk(flags)
        assert a.flags is flags

    def test_exhaustive_combinations_no_exception(self) -> None:
        for ecoli, kpneu, efaec, ssaph in product([False, True], repeat=4):
            assess_uti_risk(_flags(ecoli=ecoli, kpneu=kpneu, efaec=efaec, ssaph=ssaph))
