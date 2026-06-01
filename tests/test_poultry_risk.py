"""Tests for the poultry biosecurity LAMP panel risk assessment engine.

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

from lamp_forge.poultry_risk import (
    PoultryAlertLevel,
    PoultryAssessment,
    PoultryFlags,
    assess_poultry_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> PoultryFlags:
    defaults: dict[str, bool] = {"aiv": False, "ndv": False, "ibdv": False, "ibv": False}
    defaults.update(kwargs)
    return PoultryFlags(**defaults)


def _assess(**kwargs: bool) -> PoultryAssessment:
    return assess_poultry_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """Each pathogen must map to the correct alert level and score."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is PoultryAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_aiv_alone_is_critical(self) -> None:
        a = _assess(aiv=True)
        assert a.alert_level is PoultryAlertLevel.CRITICAL
        assert a.alert_score == 60

    def test_ndv_alone_is_high(self) -> None:
        a = _assess(ndv=True)
        assert a.alert_level is PoultryAlertLevel.HIGH
        assert a.alert_score == 35

    def test_ibdv_alone_is_moderate(self) -> None:
        a = _assess(ibdv=True)
        assert a.alert_level is PoultryAlertLevel.MODERATE
        assert a.alert_score == 20

    def test_ibv_alone_is_low(self) -> None:
        a = _assess(ibv=True)
        assert a.alert_level is PoultryAlertLevel.LOW
        assert a.alert_score == 10


# ---------------------------------------------------------------------------
# Co-detection scoring
# ---------------------------------------------------------------------------


class TestCoDetectionScores:
    """Combined flags must sum weights correctly (capped at 100)."""

    def test_aiv_ndv_score(self) -> None:
        a = _assess(aiv=True, ndv=True)
        assert a.alert_score == 95
        assert a.alert_level is PoultryAlertLevel.CRITICAL

    def test_aiv_ibdv_score(self) -> None:
        a = _assess(aiv=True, ibdv=True)
        assert a.alert_score == 80
        assert a.alert_level is PoultryAlertLevel.CRITICAL

    def test_aiv_ibv_score(self) -> None:
        a = _assess(aiv=True, ibv=True)
        assert a.alert_score == 70
        assert a.alert_level is PoultryAlertLevel.CRITICAL

    def test_ndv_ibdv_score(self) -> None:
        a = _assess(ndv=True, ibdv=True)
        assert a.alert_score == 55
        assert a.alert_level is PoultryAlertLevel.HIGH

    def test_ndv_ibv_score(self) -> None:
        a = _assess(ndv=True, ibv=True)
        assert a.alert_score == 45
        assert a.alert_level is PoultryAlertLevel.HIGH

    def test_ibdv_ibv_score(self) -> None:
        a = _assess(ibdv=True, ibv=True)
        assert a.alert_score == 30
        assert a.alert_level is PoultryAlertLevel.MODERATE

    def test_all_positive_score_capped_at_100(self) -> None:
        a = _assess(aiv=True, ndv=True, ibdv=True, ibv=True)
        assert a.alert_score == 100
        assert a.alert_level is PoultryAlertLevel.CRITICAL

    def test_score_never_exceeds_100(self) -> None:
        for aiv, ndv, ibdv, ibv in product([False, True], repeat=4):
            a = _assess(aiv=aiv, ndv=ndv, ibdv=ibdv, ibv=ibv)
            assert 0 <= a.alert_score <= 100

    def test_score_never_below_zero(self) -> None:
        for aiv, ndv, ibdv, ibv in product([False, True], repeat=4):
            a = _assess(aiv=aiv, ndv=ndv, ibdv=ibdv, ibv=ibv)
            assert a.alert_score >= 0


# ---------------------------------------------------------------------------
# Regulatory and risk flags
# ---------------------------------------------------------------------------


class TestRegulatoryFlags:
    """immediate_report_required, zoonotic_risk, depopulation_risk semantics."""

    def test_aiv_triggers_all_flags(self) -> None:
        a = _assess(aiv=True)
        assert a.immediate_report_required is True
        assert a.zoonotic_risk is True
        assert a.depopulation_risk is True

    def test_ndv_triggers_report_and_depopulation_not_zoonotic(self) -> None:
        a = _assess(ndv=True)
        assert a.immediate_report_required is True
        assert a.zoonotic_risk is False
        assert a.depopulation_risk is True

    def test_ibdv_no_regulatory_flags(self) -> None:
        a = _assess(ibdv=True)
        assert a.immediate_report_required is False
        assert a.zoonotic_risk is False
        assert a.depopulation_risk is False

    def test_ibv_no_regulatory_flags(self) -> None:
        a = _assess(ibv=True)
        assert a.immediate_report_required is False
        assert a.zoonotic_risk is False
        assert a.depopulation_risk is False

    def test_all_negative_no_flags(self) -> None:
        a = _assess()
        assert a.immediate_report_required is False
        assert a.zoonotic_risk is False
        assert a.depopulation_risk is False

    def test_aiv_ndv_co_all_flags(self) -> None:
        a = _assess(aiv=True, ndv=True)
        assert a.immediate_report_required is True
        assert a.zoonotic_risk is True
        assert a.depopulation_risk is True

    def test_ndv_ibdv_co_report_and_depopulation(self) -> None:
        a = _assess(ndv=True, ibdv=True)
        assert a.immediate_report_required is True
        assert a.zoonotic_risk is False
        assert a.depopulation_risk is True

    def test_ibdv_ibv_co_no_regulatory_flags(self) -> None:
        a = _assess(ibdv=True, ibv=True)
        assert a.immediate_report_required is False
        assert a.zoonotic_risk is False
        assert a.depopulation_risk is False


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    """active_targets must list positive detections in canonical order."""

    def test_all_negative_empty(self) -> None:
        assert _assess().active_targets == ()

    def test_aiv_only(self) -> None:
        assert _assess(aiv=True).active_targets == ("AIV",)

    def test_ndv_only(self) -> None:
        assert _assess(ndv=True).active_targets == ("NDV",)

    def test_ibdv_only(self) -> None:
        assert _assess(ibdv=True).active_targets == ("IBDV",)

    def test_ibv_only(self) -> None:
        assert _assess(ibv=True).active_targets == ("IBV",)

    def test_all_positive_canonical_order(self) -> None:
        a = _assess(aiv=True, ndv=True, ibdv=True, ibv=True)
        assert a.active_targets == ("AIV", "NDV", "IBDV", "IBV")

    def test_ndv_ibv_order(self) -> None:
        a = _assess(ndv=True, ibv=True)
        assert a.active_targets == ("NDV", "IBV")

    def test_ibdv_ibv_order(self) -> None:
        a = _assess(ibdv=True, ibv=True)
        assert a.active_targets == ("IBDV", "IBV")

    def test_aiv_ibv_order(self) -> None:
        a = _assess(aiv=True, ibv=True)
        assert a.active_targets == ("AIV", "IBV")


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestInterpretationText:
    """Key phrases must appear in interpretation/action strings."""

    def test_negative_interpretation_mentions_no_detection(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_aiv_interpretation_mentions_influenza(self) -> None:
        a = _assess(aiv=True)
        assert "influenza" in a.interpretation.lower() or "AIV" in a.interpretation

    def test_aiv_action_mentions_critical(self) -> None:
        a = _assess(aiv=True)
        assert "CRITICAL" in a.recommended_action or "critical" in a.recommended_action.lower()

    def test_aiv_action_mentions_national_authority(self) -> None:
        a = _assess(aiv=True)
        text = a.recommended_action.lower()
        assert "national" in text or "authority" in text

    def test_ndv_interpretation_mentions_newcastle(self) -> None:
        a = _assess(ndv=True)
        text = a.interpretation.lower()
        assert "newcastle" in text or "ndv" in text.lower()

    def test_ndv_action_mentions_quarantine(self) -> None:
        a = _assess(ndv=True)
        assert "quarantine" in a.recommended_action.lower()

    def test_ibdv_interpretation_mentions_bursal(self) -> None:
        a = _assess(ibdv=True)
        text = a.interpretation.lower()
        assert "bursal" in text or "gumboro" in text or "ibdv" in text

    def test_ibv_interpretation_mentions_bronchitis(self) -> None:
        a = _assess(ibv=True)
        text = a.interpretation.lower()
        assert "bronchitis" in text or "ibv" in text

    def test_aiv_ndv_co_mentions_both(self) -> None:
        a = _assess(aiv=True, ndv=True)
        assert "NDV" in a.interpretation

    def test_aiv_ibdv_co_mentions_ibdv(self) -> None:
        a = _assess(aiv=True, ibdv=True)
        assert "IBDV" in a.interpretation

    def test_ndv_ibdv_co_mentions_ibdv(self) -> None:
        a = _assess(ndv=True, ibdv=True)
        assert "IBDV" in a.interpretation

    def test_ibdv_ibv_co_mentions_ibv(self) -> None:
        a = _assess(ibdv=True, ibv=True)
        assert "IBV" in a.interpretation


# ---------------------------------------------------------------------------
# Serialisation: to_dict
# ---------------------------------------------------------------------------


class TestToDict:
    """to_dict must include all required keys with correct types."""

    def test_required_keys_present(self) -> None:
        d = _assess(aiv=True).to_dict()
        required = {
            "alert_level",
            "alert_score",
            "active_targets",
            "interpretation",
            "recommended_action",
            "immediate_report_required",
            "zoonotic_risk",
            "depopulation_risk",
            "panel_flags",
        }
        assert required <= d.keys()

    def test_alert_level_is_string(self) -> None:
        d = _assess(aiv=True).to_dict()
        assert isinstance(d["alert_level"], str)
        assert d["alert_level"] == "CRITICAL"

    def test_panel_flags_keys(self) -> None:
        d = _assess(aiv=True, ndv=False, ibdv=True, ibv=False).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf == {"aiv": True, "ndv": False, "ibdv": True, "ibv": False}

    def test_active_targets_is_list(self) -> None:
        d = _assess(aiv=True, ibv=True).to_dict()
        assert isinstance(d["active_targets"], list)
        assert "AIV" in d["active_targets"]

    def test_negative_panel_flags_all_false(self) -> None:
        d = _assess().to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert all(v is False for v in pf.values())


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    """flags_from_dict must parse dicts correctly and tolerate missing keys."""

    def test_all_true(self) -> None:
        f = flags_from_dict({"aiv": True, "ndv": True, "ibdv": True, "ibv": True})
        assert f == PoultryFlags(aiv=True, ndv=True, ibdv=True, ibv=True)

    def test_all_false_explicit(self) -> None:
        f = flags_from_dict({"aiv": False, "ndv": False, "ibdv": False, "ibv": False})
        assert f == PoultryFlags(aiv=False, ndv=False, ibdv=False, ibv=False)

    def test_empty_dict_defaults_to_false(self) -> None:
        f = flags_from_dict({})
        assert f == PoultryFlags(aiv=False, ndv=False, ibdv=False, ibv=False)

    def test_partial_dict(self) -> None:
        f = flags_from_dict({"aiv": True})
        assert f == PoultryFlags(aiv=True, ndv=False, ibdv=False, ibv=False)

    def test_unknown_keys_ignored(self) -> None:
        f = flags_from_dict({"aiv": True, "mystery_target": True})
        assert f == PoultryFlags(aiv=True, ndv=False, ibdv=False, ibv=False)

    def test_truthy_int_converted(self) -> None:
        f = flags_from_dict({"aiv": 1, "ndv": 0, "ibdv": 1, "ibv": 0})
        assert f.aiv is True
        assert f.ndv is False
        assert f.ibdv is True
        assert f.ibv is False


# ---------------------------------------------------------------------------
# File output: write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    """write_assessment_json must produce valid, parseable JSON."""

    def test_json_round_trip(self, tmp_path: Path) -> None:
        assessment = _assess(aiv=True, ibdv=True)
        out = tmp_path / "out.json"
        write_assessment_json(assessment, out)
        assert out.exists()
        data = json.loads(out.read_text(encoding="utf-8"))
        assert data["alert_level"] == "CRITICAL"
        assert data["panel_flags"]["aiv"] is True
        assert data["panel_flags"]["ibdv"] is True

    def test_json_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(_assess(), out)
        assert out.exists()

    def test_json_ends_with_newline(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        write_assessment_json(_assess(ndv=True), out)
        content = out.read_bytes()
        assert content.endswith(b"\n")

    def test_json_negative_assessment(self, tmp_path: Path) -> None:
        out = tmp_path / "neg.json"
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
        write_assessment_csv(_assess(ibv=True), out)
        with out.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert header == ["field", "value"]

    def test_csv_contains_alert_level(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(ibv=True), out)
        content = out.read_text(encoding="utf-8")
        assert "LOW" in content

    def test_csv_target_flags_positive_negative(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(aiv=True, ndv=False, ibdv=False, ibv=False), out)
        content = out.read_text(encoding="utf-8")
        assert "positive" in content
        assert "negative" in content

    def test_csv_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "deep" / "nested" / "out.csv"
        write_assessment_csv(_assess(), out)
        assert out.exists()

    def test_csv_contains_zoonotic_risk(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(aiv=True), out)
        content = out.read_text(encoding="utf-8")
        assert "zoonotic_risk" in content


# ---------------------------------------------------------------------------
# Return type and immutability
# ---------------------------------------------------------------------------


class TestReturnTypeAndImmutability:
    """assess_poultry_risk must return a frozen PoultryAssessment."""

    def test_returns_correct_type(self) -> None:
        a = _assess(aiv=True)
        assert isinstance(a, PoultryAssessment)

    def test_frozen_dataclass(self) -> None:
        a = _assess(ndv=True)
        with pytest.raises((AttributeError, TypeError)):
            a.alert_score = 0  # type: ignore[misc]

    def test_flags_field_matches_input(self) -> None:
        flags = _flags(aiv=False, ndv=True, ibdv=False, ibv=False)
        a = assess_poultry_risk(flags)
        assert a.flags is flags

    def test_exhaustive_combinations_no_exception(self) -> None:
        for aiv, ndv, ibdv, ibv in product([False, True], repeat=4):
            assess_poultry_risk(_flags(aiv=aiv, ndv=ndv, ibdv=ibdv, ibv=ibv))

    def test_interpretation_and_action_are_nonempty_strings(self) -> None:
        for aiv, ndv, ibdv, ibv in product([False, True], repeat=4):
            a = assess_poultry_risk(_flags(aiv=aiv, ndv=ndv, ibdv=ibdv, ibv=ibv))
            assert isinstance(a.interpretation, str) and len(a.interpretation) > 0
            assert isinstance(a.recommended_action, str) and len(a.recommended_action) > 0
