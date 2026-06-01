"""Tests for the swine neonatal enteric LAMP panel risk assessment engine.

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

from lamp_forge.swine_enteric_risk import (
    SwineEntericAlertLevel,
    SwineEntericAssessment,
    SwineEntericFlags,
    assess_swine_enteric_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> SwineEntericFlags:
    defaults: dict[str, bool] = {"pedv": False, "pdcov": False, "rota_a": False}
    defaults.update(kwargs)
    return SwineEntericFlags(**defaults)


def _assess(**kwargs: bool) -> SwineEntericAssessment:
    return assess_swine_enteric_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """Each pathogen must map to the correct alert level and score."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is SwineEntericAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_pedv_alone_is_critical(self) -> None:
        a = _assess(pedv=True)
        assert a.alert_level is SwineEntericAlertLevel.CRITICAL
        assert a.alert_score == 50

    def test_pdcov_alone_is_high(self) -> None:
        a = _assess(pdcov=True)
        assert a.alert_level is SwineEntericAlertLevel.HIGH
        assert a.alert_score == 30

    def test_rota_a_alone_is_moderate(self) -> None:
        a = _assess(rota_a=True)
        assert a.alert_level is SwineEntericAlertLevel.MODERATE
        assert a.alert_score == 15


# ---------------------------------------------------------------------------
# Co-detection scoring
# ---------------------------------------------------------------------------


class TestCoDetectionScores:
    """Combined flags must sum weights correctly (capped at 100)."""

    def test_pedv_pdcov_score(self) -> None:
        a = _assess(pedv=True, pdcov=True)
        assert a.alert_score == 80
        assert a.alert_level is SwineEntericAlertLevel.CRITICAL

    def test_pedv_rota_a_score(self) -> None:
        a = _assess(pedv=True, rota_a=True)
        assert a.alert_score == 65
        assert a.alert_level is SwineEntericAlertLevel.CRITICAL

    def test_pdcov_rota_a_score(self) -> None:
        a = _assess(pdcov=True, rota_a=True)
        assert a.alert_score == 45
        assert a.alert_level is SwineEntericAlertLevel.HIGH

    def test_all_positive_score_capped(self) -> None:
        a = _assess(pedv=True, pdcov=True, rota_a=True)
        assert a.alert_score == 95
        assert a.alert_level is SwineEntericAlertLevel.CRITICAL

    def test_score_never_exceeds_100(self) -> None:
        for pedv, pdcov, rota_a in product([False, True], repeat=3):
            a = _assess(pedv=pedv, pdcov=pdcov, rota_a=rota_a)
            assert 0 <= a.alert_score <= 100


# ---------------------------------------------------------------------------
# Quarantine and biosecurity flags
# ---------------------------------------------------------------------------


class TestQuarantineFlags:
    """immediate_quarantine_required and biosecurity_lockdown semantics."""

    def test_pedv_triggers_quarantine_and_lockdown(self) -> None:
        a = _assess(pedv=True)
        assert a.immediate_quarantine_required is True
        assert a.biosecurity_lockdown is True

    def test_pdcov_no_quarantine_but_lockdown(self) -> None:
        a = _assess(pdcov=True)
        assert a.immediate_quarantine_required is False
        assert a.biosecurity_lockdown is True

    def test_rota_a_no_quarantine_no_lockdown(self) -> None:
        a = _assess(rota_a=True)
        assert a.immediate_quarantine_required is False
        assert a.biosecurity_lockdown is False

    def test_all_negative_no_flags(self) -> None:
        a = _assess()
        assert a.immediate_quarantine_required is False
        assert a.biosecurity_lockdown is False

    def test_pedv_pdcov_co_quarantine_and_lockdown(self) -> None:
        a = _assess(pedv=True, pdcov=True)
        assert a.immediate_quarantine_required is True
        assert a.biosecurity_lockdown is True

    def test_pdcov_rota_a_lockdown_no_quarantine(self) -> None:
        a = _assess(pdcov=True, rota_a=True)
        assert a.immediate_quarantine_required is False
        assert a.biosecurity_lockdown is True


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    """active_targets must list positive detections in canonical order."""

    def test_all_negative_empty(self) -> None:
        assert _assess().active_targets == ()

    def test_pedv_only(self) -> None:
        assert _assess(pedv=True).active_targets == ("PEDV",)

    def test_pdcov_only(self) -> None:
        assert _assess(pdcov=True).active_targets == ("PDCoV",)

    def test_rota_a_only(self) -> None:
        assert _assess(rota_a=True).active_targets == ("RotaA",)

    def test_all_positive_canonical_order(self) -> None:
        a = _assess(pedv=True, pdcov=True, rota_a=True)
        assert a.active_targets == ("PEDV", "PDCoV", "RotaA")

    def test_pdcov_rota_a_order(self) -> None:
        a = _assess(pdcov=True, rota_a=True)
        assert a.active_targets == ("PDCoV", "RotaA")


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestInterpretationText:
    """Key phrases must appear in interpretation/action strings."""

    def test_negative_interpretation_mentions_no_pathogens(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_pedv_interpretation_mentions_pedv(self) -> None:
        a = _assess(pedv=True)
        assert "PEDV" in a.interpretation or "epidemic diarrhea" in a.interpretation.lower()

    def test_pedv_action_mentions_quarantine(self) -> None:
        a = _assess(pedv=True)
        assert "quarantine" in a.recommended_action.lower()

    def test_pdcov_interpretation_mentions_pdcov(self) -> None:
        a = _assess(pdcov=True)
        assert "PDCoV" in a.interpretation or "deltacoronavirus" in a.interpretation.lower()

    def test_rota_a_interpretation_mentions_rotavirus(self) -> None:
        a = _assess(rota_a=True)
        assert "rotavirus" in a.interpretation.lower()

    def test_pedv_pdcov_co_interpretation_mentions_both(self) -> None:
        a = _assess(pedv=True, pdcov=True)
        assert "PDCoV" in a.interpretation

    def test_pedv_rota_a_co_interpretation_mentions_rotavirus(self) -> None:
        a = _assess(pedv=True, rota_a=True)
        assert "Rotavirus" in a.interpretation


# ---------------------------------------------------------------------------
# Serialisation: to_dict
# ---------------------------------------------------------------------------


class TestToDict:
    """to_dict must include all required keys with correct types."""

    def test_required_keys_present(self) -> None:
        d = _assess(pedv=True).to_dict()
        required = {
            "alert_level",
            "alert_score",
            "active_targets",
            "interpretation",
            "recommended_action",
            "immediate_quarantine_required",
            "biosecurity_lockdown",
            "panel_flags",
        }
        assert required <= d.keys()

    def test_alert_level_is_string(self) -> None:
        d = _assess(pedv=True).to_dict()
        assert isinstance(d["alert_level"], str)
        assert d["alert_level"] == "CRITICAL"

    def test_panel_flags_keys(self) -> None:
        d = _assess(pedv=True, pdcov=False, rota_a=True).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf == {"pedv": True, "pdcov": False, "rota_a": True}

    def test_active_targets_is_list(self) -> None:
        d = _assess(pedv=True, rota_a=True).to_dict()
        assert isinstance(d["active_targets"], list)
        assert "PEDV" in d["active_targets"]


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    """flags_from_dict must parse dicts correctly and tolerate missing keys."""

    def test_all_true(self) -> None:
        f = flags_from_dict({"pedv": True, "pdcov": True, "rota_a": True})
        assert f == SwineEntericFlags(pedv=True, pdcov=True, rota_a=True)

    def test_all_false_explicit(self) -> None:
        f = flags_from_dict({"pedv": False, "pdcov": False, "rota_a": False})
        assert f == SwineEntericFlags(pedv=False, pdcov=False, rota_a=False)

    def test_empty_dict_defaults_to_false(self) -> None:
        f = flags_from_dict({})
        assert f == SwineEntericFlags(pedv=False, pdcov=False, rota_a=False)

    def test_partial_dict(self) -> None:
        f = flags_from_dict({"pedv": True})
        assert f == SwineEntericFlags(pedv=True, pdcov=False, rota_a=False)

    def test_unknown_keys_ignored(self) -> None:
        f = flags_from_dict({"pedv": True, "unknown_key": True})
        assert f == SwineEntericFlags(pedv=True, pdcov=False, rota_a=False)

    def test_truthy_int_converted(self) -> None:
        f = flags_from_dict({"pedv": 1, "pdcov": 0, "rota_a": 1})
        assert f.pedv is True
        assert f.pdcov is False
        assert f.rota_a is True


# ---------------------------------------------------------------------------
# File output: write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    """write_assessment_json must produce valid, parseable JSON."""

    def test_json_round_trip(self, tmp_path: Path) -> None:
        assessment = _assess(pedv=True, rota_a=True)
        out = tmp_path / "out.json"
        write_assessment_json(assessment, out)
        assert out.exists()
        data = json.loads(out.read_text(encoding="utf-8"))
        assert data["alert_level"] == "CRITICAL"
        assert data["panel_flags"]["pedv"] is True
        assert data["panel_flags"]["rota_a"] is True

    def test_json_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(_assess(), out)
        assert out.exists()

    def test_json_ends_with_newline(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        write_assessment_json(_assess(pdcov=True), out)
        content = out.read_bytes()
        assert content.endswith(b"\n")


# ---------------------------------------------------------------------------
# File output: write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    """write_assessment_csv must produce a valid two-column key-value CSV."""

    def test_csv_has_header(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(rota_a=True), out)
        with out.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert header == ["field", "value"]

    def test_csv_contains_alert_level(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(rota_a=True), out)
        content = out.read_text(encoding="utf-8")
        assert "MODERATE" in content

    def test_csv_target_flags_positive_negative(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(pedv=True, pdcov=False, rota_a=False), out)
        content = out.read_text(encoding="utf-8")
        assert "positive" in content
        assert "negative" in content

    def test_csv_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "deep" / "nested" / "out.csv"
        write_assessment_csv(_assess(), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Return type and immutability
# ---------------------------------------------------------------------------


class TestReturnTypeAndImmutability:
    """assess_swine_enteric_risk must return a frozen SwineEntericAssessment."""

    def test_returns_correct_type(self) -> None:
        a = _assess(pedv=True)
        assert isinstance(a, SwineEntericAssessment)

    def test_frozen_dataclass(self) -> None:
        a = _assess(pdcov=True)
        with pytest.raises((AttributeError, TypeError)):
            a.alert_score = 0  # type: ignore[misc]

    def test_flags_field_matches_input(self) -> None:
        flags = _flags(pedv=False, pdcov=True, rota_a=False)
        a = assess_swine_enteric_risk(flags)
        assert a.flags is flags

    def test_exhaustive_combinations_no_exception(self) -> None:
        for pedv, pdcov, rota_a in product([False, True], repeat=3):
            assess_swine_enteric_risk(_flags(pedv=pedv, pdcov=pdcov, rota_a=rota_a))
