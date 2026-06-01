"""Tests for lamp_forge.brucella_risk -- Brucella IS711 LAMP risk assessment.

All tests are pure Python; no external binaries (MAFFT, BLAST) or network
calls are required.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from lamp_forge.brucella_risk import (
    BrucellaAlertLevel,
    BrucellaAssessment,
    BrucellaContext,
    BrucellaFlags,
    assess_brucella_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(positive: bool, context: BrucellaContext = BrucellaContext.CATTLE) -> BrucellaFlags:
    return BrucellaFlags(is711_positive=positive, context=context)


def _assess(
    positive: bool, context: BrucellaContext = BrucellaContext.CATTLE
) -> BrucellaAssessment:
    return assess_brucella_risk(_flags(positive, context))


# ---------------------------------------------------------------------------
# BrucellaFlags defaults
# ---------------------------------------------------------------------------


class TestBrucellaFlagsDefaults:
    def test_default_context_is_cattle(self) -> None:
        f = BrucellaFlags(is711_positive=False)
        assert f.context is BrucellaContext.CATTLE

    def test_frozen(self) -> None:
        f = BrucellaFlags(is711_positive=True)
        with pytest.raises(AttributeError):
            f.is711_positive = False  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Negative result -- all contexts
# ---------------------------------------------------------------------------


class TestNegativeResult:
    @pytest.mark.parametrize("ctx", list(BrucellaContext))
    def test_negative_gives_negative_level(self, ctx: BrucellaContext) -> None:
        assessment = _assess(False, ctx)
        assert assessment.alert_level is BrucellaAlertLevel.NEGATIVE

    @pytest.mark.parametrize("ctx", list(BrucellaContext))
    def test_negative_score_is_zero(self, ctx: BrucellaContext) -> None:
        assessment = _assess(False, ctx)
        assert assessment.alert_score == 0

    @pytest.mark.parametrize("ctx", list(BrucellaContext))
    def test_negative_not_notifiable(self, ctx: BrucellaContext) -> None:
        assessment = _assess(False, ctx)
        assert not assessment.notifiable

    @pytest.mark.parametrize("ctx", list(BrucellaContext))
    def test_negative_no_zoonotic_risk(self, ctx: BrucellaContext) -> None:
        assessment = _assess(False, ctx)
        assert not assessment.zoonotic_risk

    @pytest.mark.parametrize("ctx", list(BrucellaContext))
    def test_negative_no_immediate_report(self, ctx: BrucellaContext) -> None:
        assessment = _assess(False, ctx)
        assert not assessment.immediate_report_required


# ---------------------------------------------------------------------------
# Human context -- CRITICAL
# ---------------------------------------------------------------------------


class TestHumanContext:
    def test_positive_human_is_critical(self) -> None:
        assert _assess(True, BrucellaContext.HUMAN).alert_level is BrucellaAlertLevel.CRITICAL

    def test_positive_human_score_is_90(self) -> None:
        assert _assess(True, BrucellaContext.HUMAN).alert_score == 90

    def test_positive_human_is_notifiable(self) -> None:
        assert _assess(True, BrucellaContext.HUMAN).notifiable is True

    def test_positive_human_zoonotic_risk(self) -> None:
        assert _assess(True, BrucellaContext.HUMAN).zoonotic_risk is True

    def test_positive_human_immediate_report(self) -> None:
        assert _assess(True, BrucellaContext.HUMAN).immediate_report_required is True

    def test_human_interpretation_mentions_brucellosis(self) -> None:
        interp = _assess(True, BrucellaContext.HUMAN).interpretation
        assert "brucellosis" in interp.lower() or "brucella" in interp.lower()

    def test_human_action_mentions_therapy(self) -> None:
        action = _assess(True, BrucellaContext.HUMAN).recommended_action
        assert "doxycycline" in action.lower()


# ---------------------------------------------------------------------------
# Livestock contexts -- HIGH
# ---------------------------------------------------------------------------


class TestLivestockContexts:
    @pytest.mark.parametrize(
        "ctx",
        [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE],
    )
    def test_positive_livestock_is_high(self, ctx: BrucellaContext) -> None:
        assert _assess(True, ctx).alert_level is BrucellaAlertLevel.HIGH

    @pytest.mark.parametrize(
        "ctx",
        [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE],
    )
    def test_positive_livestock_score_is_60(self, ctx: BrucellaContext) -> None:
        assert _assess(True, ctx).alert_score == 60

    @pytest.mark.parametrize(
        "ctx",
        [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE],
    )
    def test_livestock_is_notifiable(self, ctx: BrucellaContext) -> None:
        assert _assess(True, ctx).notifiable is True

    @pytest.mark.parametrize(
        "ctx",
        [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE],
    )
    def test_livestock_zoonotic_risk(self, ctx: BrucellaContext) -> None:
        assert _assess(True, ctx).zoonotic_risk is True

    @pytest.mark.parametrize(
        "ctx",
        [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE],
    )
    def test_livestock_immediate_report(self, ctx: BrucellaContext) -> None:
        assert _assess(True, ctx).immediate_report_required is True

    def test_cattle_interpretation_mentions_abortus(self) -> None:
        interp = _assess(True, BrucellaContext.CATTLE).interpretation
        assert "abortus" in interp.lower()

    def test_cattle_action_mentions_quarantine(self) -> None:
        action = _assess(True, BrucellaContext.CATTLE).recommended_action
        assert "quarantine" in action.lower()

    def test_small_ruminant_interpretation_mentions_melitensis(self) -> None:
        interp = _assess(True, BrucellaContext.SMALL_RUMINANT).interpretation
        assert "melitensis" in interp.lower()

    def test_swine_interpretation_mentions_suis(self) -> None:
        interp = _assess(True, BrucellaContext.SWINE).interpretation
        assert "suis" in interp.lower()


# ---------------------------------------------------------------------------
# Environmental context -- LOW
# ---------------------------------------------------------------------------


class TestEnvironmentalContext:
    def test_positive_environmental_is_low(self) -> None:
        assert _assess(True, BrucellaContext.ENVIRONMENTAL).alert_level is BrucellaAlertLevel.LOW

    def test_positive_environmental_score_is_25(self) -> None:
        assert _assess(True, BrucellaContext.ENVIRONMENTAL).alert_score == 25

    def test_environmental_not_notifiable(self) -> None:
        assert _assess(True, BrucellaContext.ENVIRONMENTAL).notifiable is False

    def test_environmental_no_zoonotic_risk(self) -> None:
        assert _assess(True, BrucellaContext.ENVIRONMENTAL).zoonotic_risk is False

    def test_environmental_no_immediate_report(self) -> None:
        assert _assess(True, BrucellaContext.ENVIRONMENTAL).immediate_report_required is False

    def test_environmental_action_says_no_report_alone(self) -> None:
        action = _assess(True, BrucellaContext.ENVIRONMENTAL).recommended_action
        assert "do not report" in action.lower()


# ---------------------------------------------------------------------------
# to_dict / serialisation round-trip
# ---------------------------------------------------------------------------


class TestToDict:
    def test_to_dict_has_required_keys(self) -> None:
        d = _assess(True, BrucellaContext.CATTLE).to_dict()
        for key in (
            "alert_level",
            "alert_score",
            "interpretation",
            "recommended_action",
            "notifiable",
            "zoonotic_risk",
            "immediate_report_required",
            "is711_result",
            "context",
        ):
            assert key in d, f"Missing key: {key}"

    def test_to_dict_positive_is711_result(self) -> None:
        d = _assess(True, BrucellaContext.CATTLE).to_dict()
        assert d["is711_result"] == "positive"

    def test_to_dict_negative_is711_result(self) -> None:
        d = _assess(False, BrucellaContext.CATTLE).to_dict()
        assert d["is711_result"] == "negative"

    def test_to_dict_context_value(self) -> None:
        d = _assess(True, BrucellaContext.HUMAN).to_dict()
        assert d["context"] == "human"

    def test_to_dict_alert_level_is_string(self) -> None:
        d = _assess(True, BrucellaContext.CATTLE).to_dict()
        assert isinstance(d["alert_level"], str)

    def test_to_dict_json_serialisable(self) -> None:
        d = _assess(True, BrucellaContext.HUMAN).to_dict()
        json.dumps(d)  # must not raise


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_positive_is711(self) -> None:
        f = flags_from_dict({"is711_positive": True, "context": "cattle"})
        assert f.is711_positive is True
        assert f.context is BrucellaContext.CATTLE

    def test_missing_keys_use_defaults(self) -> None:
        f = flags_from_dict({})
        assert f.is711_positive is False
        assert f.context is BrucellaContext.CATTLE

    def test_all_valid_contexts(self) -> None:
        for ctx in BrucellaContext:
            f = flags_from_dict({"is711_positive": False, "context": ctx.value})
            assert f.context is ctx

    def test_invalid_context_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid context"):
            flags_from_dict({"context": "invalid_species"})

    def test_unknown_keys_are_ignored(self) -> None:
        f = flags_from_dict({"is711_positive": True, "context": "human", "extra": 99})
        assert f.is711_positive is True


# ---------------------------------------------------------------------------
# File output helpers
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_writes_valid_json(self, tmp_path: Path) -> None:
        assessment = _assess(True, BrucellaContext.CATTLE)
        out = tmp_path / "brucella.json"
        write_assessment_json(assessment, out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "HIGH"

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        assessment = _assess(False)
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(assessment, out)
        assert out.exists()

    def test_file_ends_with_newline(self, tmp_path: Path) -> None:
        assessment = _assess(True, BrucellaContext.HUMAN)
        out = tmp_path / "out.json"
        write_assessment_json(assessment, out)
        assert out.read_text().endswith("\n")


class TestWriteAssessmentCsv:
    def test_writes_csv_with_header(self, tmp_path: Path) -> None:
        assessment = _assess(True, BrucellaContext.SWINE)
        out = tmp_path / "brucella.csv"
        write_assessment_csv(assessment, out)
        lines = out.read_text().splitlines()
        assert lines[0] == "field,value"

    def test_csv_has_alert_level_row(self, tmp_path: Path) -> None:
        assessment = _assess(True, BrucellaContext.SWINE)
        out = tmp_path / "brucella.csv"
        write_assessment_csv(assessment, out)
        content = out.read_text()
        assert "alert_level" in content
        assert "HIGH" in content

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "a" / "b" / "brucella.csv"
        write_assessment_csv(_assess(False), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Ordering: CRITICAL > HIGH > LOW > NEGATIVE
# ---------------------------------------------------------------------------


class TestAlertLevelOrdering:
    def test_human_positive_highest_score(self) -> None:
        human_score = _assess(True, BrucellaContext.HUMAN).alert_score
        for ctx in [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE]:
            assert human_score > _assess(True, ctx).alert_score

    def test_livestock_above_environmental(self) -> None:
        env_score = _assess(True, BrucellaContext.ENVIRONMENTAL).alert_score
        for ctx in [BrucellaContext.CATTLE, BrucellaContext.SMALL_RUMINANT, BrucellaContext.SWINE]:
            assert _assess(True, ctx).alert_score > env_score

    def test_environmental_above_negative(self) -> None:
        assert _assess(True, BrucellaContext.ENVIRONMENTAL).alert_score > _assess(False).alert_score

    def test_negative_score_is_zero(self) -> None:
        assert _assess(False).alert_score == 0


# ---------------------------------------------------------------------------
# Flags stored on assessment
# ---------------------------------------------------------------------------


class TestAssessmentFlagsPreserved:
    def test_flags_is711_preserved(self) -> None:
        assessment = _assess(True, BrucellaContext.SMALL_RUMINANT)
        assert assessment.flags.is711_positive is True
        assert assessment.flags.context is BrucellaContext.SMALL_RUMINANT

    def test_assessment_is_frozen(self) -> None:
        assessment = _assess(True)
        with pytest.raises(AttributeError):
            assessment.alert_score = 0  # type: ignore[misc]
