"""Tests for the swine respiratory / PRDC LAMP panel risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical alert levels.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

from lamp_forge.swine_respiratory_risk import (
    PrdcAlertLevel,
    PrdcAssessment,
    PrdcFlags,
    assess_prdc_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> PrdcFlags:
    defaults: dict[str, bool] = {"prrsv": False, "pcv2": False, "mhp": False}
    defaults.update(kwargs)
    return PrdcFlags(**defaults)


def _assess(**kwargs: bool) -> PrdcAssessment:
    return assess_prdc_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """Each pathogen must map to the correct alert level and score."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is PrdcAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_prrsv_alone_is_high(self) -> None:
        a = _assess(prrsv=True)
        assert a.alert_level is PrdcAlertLevel.HIGH
        assert a.alert_score == 40

    def test_pcv2_alone_is_moderate(self) -> None:
        a = _assess(pcv2=True)
        assert a.alert_level is PrdcAlertLevel.MODERATE
        assert a.alert_score == 30

    def test_mhp_alone_is_low(self) -> None:
        a = _assess(mhp=True)
        assert a.alert_level is PrdcAlertLevel.LOW
        assert a.alert_score == 20


# ---------------------------------------------------------------------------
# Alert level and score for dual-target patterns
# ---------------------------------------------------------------------------


class TestDualTargetAlertLevels:
    """Dual-positive patterns must map to the correct alert level and score."""

    def test_prrsv_pcv2_is_critical(self) -> None:
        a = _assess(prrsv=True, pcv2=True)
        assert a.alert_level is PrdcAlertLevel.CRITICAL
        assert a.alert_score == 70

    def test_prrsv_mhp_is_critical(self) -> None:
        a = _assess(prrsv=True, mhp=True)
        assert a.alert_level is PrdcAlertLevel.CRITICAL
        assert a.alert_score == 60

    def test_pcv2_mhp_is_high(self) -> None:
        a = _assess(pcv2=True, mhp=True)
        assert a.alert_level is PrdcAlertLevel.HIGH
        assert a.alert_score == 50


# ---------------------------------------------------------------------------
# Alert level and score for the full PRDC pattern (all three)
# ---------------------------------------------------------------------------


class TestFullPrdc:
    def test_all_positive_is_critical(self) -> None:
        a = _assess(prrsv=True, pcv2=True, mhp=True)
        assert a.alert_level is PrdcAlertLevel.CRITICAL
        assert a.alert_score == 90

    def test_all_positive_active_targets(self) -> None:
        a = _assess(prrsv=True, pcv2=True, mhp=True)
        assert set(a.active_targets) == {"PRRSV", "PCV2", "Mhp"}

    def test_all_positive_complex_prdc(self) -> None:
        a = _assess(prrsv=True, pcv2=True, mhp=True)
        assert a.complex_prdc is True

    def test_all_positive_movement_restriction(self) -> None:
        a = _assess(prrsv=True, pcv2=True, mhp=True)
        assert a.movement_restriction_required is True


# ---------------------------------------------------------------------------
# Decision flags
# ---------------------------------------------------------------------------


class TestDecisionFlags:
    """Verify movement_restriction_required and complex_prdc logic."""

    def test_movement_restriction_only_when_prrsv(self) -> None:
        for prrsv, pcv2, mhp in product([True, False], repeat=3):
            a = _assess(prrsv=prrsv, pcv2=pcv2, mhp=mhp)
            assert a.movement_restriction_required == prrsv

    def test_complex_prdc_when_two_or_more_positive(self) -> None:
        # single positive -> complex_prdc False
        for kw in [{"prrsv": True}, {"pcv2": True}, {"mhp": True}]:
            assert _assess(**kw).complex_prdc is False

        # dual positive -> complex_prdc True
        assert _assess(prrsv=True, pcv2=True).complex_prdc is True
        assert _assess(prrsv=True, mhp=True).complex_prdc is True
        assert _assess(pcv2=True, mhp=True).complex_prdc is True

    def test_no_movement_restriction_without_prrsv(self) -> None:
        for pcv2, mhp in product([True, False], repeat=2):
            a = _assess(prrsv=False, pcv2=pcv2, mhp=mhp)
            assert a.movement_restriction_required is False


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    def test_negative_has_no_active_targets(self) -> None:
        assert _assess().active_targets == ()

    def test_prrsv_only_active_target(self) -> None:
        a = _assess(prrsv=True)
        assert a.active_targets == ("PRRSV",)

    def test_pcv2_only_active_target(self) -> None:
        a = _assess(pcv2=True)
        assert a.active_targets == ("PCV2",)

    def test_mhp_only_active_target(self) -> None:
        a = _assess(mhp=True)
        assert a.active_targets == ("Mhp",)

    def test_canonical_order_all_three(self) -> None:
        a = _assess(prrsv=True, pcv2=True, mhp=True)
        assert a.active_targets == ("PRRSV", "PCV2", "Mhp")


# ---------------------------------------------------------------------------
# Interpretation and recommended_action are non-empty strings
# ---------------------------------------------------------------------------


class TestTextFields:
    def test_interpretation_is_nonempty_for_all_patterns(self) -> None:
        for prrsv, pcv2, mhp in product([True, False], repeat=3):
            a = _assess(prrsv=prrsv, pcv2=pcv2, mhp=mhp)
            assert isinstance(a.interpretation, str)
            assert len(a.interpretation) > 0

    def test_recommended_action_is_nonempty_for_all_patterns(self) -> None:
        for prrsv, pcv2, mhp in product([True, False], repeat=3):
            a = _assess(prrsv=prrsv, pcv2=pcv2, mhp=mhp)
            assert isinstance(a.recommended_action, str)
            assert len(a.recommended_action) > 0


# ---------------------------------------------------------------------------
# Score is capped at 100
# ---------------------------------------------------------------------------


class TestScoreCap:
    def test_max_score_is_90(self) -> None:
        a = _assess(prrsv=True, pcv2=True, mhp=True)
        assert a.alert_score == 90
        assert a.alert_score <= 100

    def test_min_score_is_0(self) -> None:
        a = _assess()
        assert a.alert_score == 0


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_false_by_default(self) -> None:
        f = flags_from_dict({})
        assert f == PrdcFlags(prrsv=False, pcv2=False, mhp=False)

    def test_truthy_values(self) -> None:
        f = flags_from_dict({"prrsv": True, "pcv2": False, "mhp": True})
        assert f.prrsv is True
        assert f.pcv2 is False
        assert f.mhp is True

    def test_unknown_keys_ignored(self) -> None:
        f = flags_from_dict({"prrsv": True, "unknown_key": "whatever"})
        assert f.prrsv is True
        assert f.pcv2 is False


# ---------------------------------------------------------------------------
# to_dict round-trip
# ---------------------------------------------------------------------------


class TestToDict:
    def test_roundtrip_prrsv_pcv2(self) -> None:
        a = _assess(prrsv=True, pcv2=True)
        d = a.to_dict()
        assert d["alert_level"] == "CRITICAL"
        assert d["alert_score"] == 70
        assert d["movement_restriction_required"] is True
        assert d["complex_prdc"] is True
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf["prrsv"] is True
        assert pf["pcv2"] is True
        assert pf["mhp"] is False

    def test_all_negative_dict(self) -> None:
        d = _assess().to_dict()
        assert d["alert_level"] == "NEGATIVE"
        assert d["alert_score"] == 0
        assert d["active_targets"] == []


# ---------------------------------------------------------------------------
# write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_json_file_is_valid(self, tmp_path: Path) -> None:
        a = _assess(prrsv=True, mhp=True)
        out = tmp_path / "prdc.json"
        write_assessment_json(a, out)
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert data["panel_flags"]["prrsv"] is True

    def test_parent_dirs_created(self, tmp_path: Path) -> None:
        out = tmp_path / "deep" / "nested" / "prdc.json"
        write_assessment_json(_assess(pcv2=True), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    def test_csv_contains_header_and_data(self, tmp_path: Path) -> None:
        a = _assess(prrsv=True)
        out = tmp_path / "prdc.csv"
        write_assessment_csv(a, out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["field", "value"]
        keys = {r[0] for r in rows[1:]}
        assert "alert_level" in keys
        assert "alert_score" in keys

    def test_csv_target_flags_row(self, tmp_path: Path) -> None:
        a = _assess(prrsv=True, pcv2=False, mhp=True)
        out = tmp_path / "prdc.csv"
        write_assessment_csv(a, out)
        data = {r[0]: r[1] for r in csv.reader(out.open()) if len(r) == 2}
        assert data["target_prrsv"] == "positive"
        assert data["target_pcv2"] == "negative"
        assert data["target_mhp"] == "positive"
