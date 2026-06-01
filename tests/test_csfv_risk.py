"""Tests for lamp_forge.csfv_risk -- CSFV NS5B RT-LAMP risk assessment.

All tests are pure Python; no external binaries (MAFFT, BLAST) or network
calls are required.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from lamp_forge.csfv_risk import (
    CsfvAlertLevel,
    CsfvAssessment,
    CsfvContext,
    CsfvFlags,
    assess_csfv_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(positive: bool, context: CsfvContext = CsfvContext.ON_FARM) -> CsfvFlags:
    return CsfvFlags(ns5b_positive=positive, context=context)


def _assess(positive: bool, context: CsfvContext = CsfvContext.ON_FARM) -> CsfvAssessment:
    return assess_csfv_risk(_flags(positive, context))


# ---------------------------------------------------------------------------
# CsfvFlags defaults
# ---------------------------------------------------------------------------


class TestCsfvFlagsDefaults:
    def test_default_context_is_on_farm(self) -> None:
        f = CsfvFlags(ns5b_positive=False)
        assert f.context is CsfvContext.ON_FARM

    def test_frozen(self) -> None:
        f = CsfvFlags(ns5b_positive=True)
        with pytest.raises(AttributeError):
            f.ns5b_positive = False  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Negative result -- all contexts
# ---------------------------------------------------------------------------


class TestNegativeResult:
    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_negative_gives_negative_level(self, ctx: CsfvContext) -> None:
        assert _assess(False, ctx).alert_level is CsfvAlertLevel.NEGATIVE

    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_negative_score_is_zero(self, ctx: CsfvContext) -> None:
        assert _assess(False, ctx).alert_score == 0

    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_negative_not_notifiable(self, ctx: CsfvContext) -> None:
        assert not _assess(False, ctx).woah_notifiable

    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_negative_no_immediate_report(self, ctx: CsfvContext) -> None:
        assert not _assess(False, ctx).immediate_report_required

    def test_negative_interpretation_mentions_not_detected(self) -> None:
        interp = _assess(False).interpretation
        assert "not detected" in interp.lower() or "negative" in interp.lower()

    def test_negative_action_no_quarantine(self) -> None:
        action = _assess(False).recommended_action
        assert "quarantine" not in action.lower() or "no immediate quarantine" in action.lower()


# ---------------------------------------------------------------------------
# Positive result -- all contexts give CRITICAL
# ---------------------------------------------------------------------------


class TestPositiveResult:
    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_positive_always_critical(self, ctx: CsfvContext) -> None:
        assert _assess(True, ctx).alert_level is CsfvAlertLevel.CRITICAL

    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_positive_score_is_100(self, ctx: CsfvContext) -> None:
        assert _assess(True, ctx).alert_score == 100

    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_positive_is_woah_notifiable(self, ctx: CsfvContext) -> None:
        assert _assess(True, ctx).woah_notifiable is True

    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_positive_requires_immediate_report(self, ctx: CsfvContext) -> None:
        assert _assess(True, ctx).immediate_report_required is True


# ---------------------------------------------------------------------------
# Context-specific recommended actions
# ---------------------------------------------------------------------------


class TestContextActions:
    def test_on_farm_action_mentions_quarantine(self) -> None:
        action = _assess(True, CsfvContext.ON_FARM).recommended_action
        assert "quarantine" in action.lower()

    def test_on_farm_action_mentions_national_authority(self) -> None:
        action = _assess(True, CsfvContext.ON_FARM).recommended_action
        assert "national" in action.lower()

    def test_border_post_action_mentions_consignment(self) -> None:
        action = _assess(True, CsfvContext.BORDER_POST).recommended_action
        assert "consignment" in action.lower()

    def test_border_post_action_mentions_detain(self) -> None:
        action = _assess(True, CsfvContext.BORDER_POST).recommended_action
        assert "detain" in action.lower() or "halt" in action.lower()

    def test_market_action_mentions_halt(self) -> None:
        action = _assess(True, CsfvContext.MARKET).recommended_action
        assert "halt" in action.lower() or "stop" in action.lower()

    def test_market_interpretation_mentions_trace(self) -> None:
        interp = _assess(True, CsfvContext.MARKET).interpretation
        assert "market" in interp.lower()

    def test_on_farm_interpretation_mentions_asfv_differential(self) -> None:
        interp = _assess(True, CsfvContext.ON_FARM).interpretation
        assert "asfv" in interp.lower() or "african swine fever" in interp.lower()


# ---------------------------------------------------------------------------
# to_dict / serialisation round-trip
# ---------------------------------------------------------------------------


class TestToDict:
    def test_to_dict_has_required_keys(self) -> None:
        d = _assess(True, CsfvContext.ON_FARM).to_dict()
        for key in (
            "alert_level",
            "alert_score",
            "interpretation",
            "recommended_action",
            "immediate_report_required",
            "woah_notifiable",
            "ns5b_result",
            "context",
        ):
            assert key in d, f"Missing key: {key}"

    def test_to_dict_positive_ns5b_result(self) -> None:
        d = _assess(True, CsfvContext.ON_FARM).to_dict()
        assert d["ns5b_result"] == "positive"

    def test_to_dict_negative_ns5b_result(self) -> None:
        d = _assess(False, CsfvContext.ON_FARM).to_dict()
        assert d["ns5b_result"] == "negative"

    def test_to_dict_context_value(self) -> None:
        d = _assess(True, CsfvContext.BORDER_POST).to_dict()
        assert d["context"] == "border_post"

    def test_to_dict_alert_level_is_string(self) -> None:
        d = _assess(True, CsfvContext.ON_FARM).to_dict()
        assert isinstance(d["alert_level"], str)

    def test_to_dict_json_serialisable(self) -> None:
        d = _assess(True, CsfvContext.MARKET).to_dict()
        json.dumps(d)  # must not raise


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_positive_ns5b(self) -> None:
        f = flags_from_dict({"ns5b_positive": True, "context": "on_farm"})
        assert f.ns5b_positive is True
        assert f.context is CsfvContext.ON_FARM

    def test_missing_keys_use_defaults(self) -> None:
        f = flags_from_dict({})
        assert f.ns5b_positive is False
        assert f.context is CsfvContext.ON_FARM

    def test_all_valid_contexts(self) -> None:
        for ctx in CsfvContext:
            f = flags_from_dict({"ns5b_positive": False, "context": ctx.value})
            assert f.context is ctx

    def test_invalid_context_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid context"):
            flags_from_dict({"context": "cattle_farm"})

    def test_unknown_keys_are_ignored(self) -> None:
        f = flags_from_dict({"ns5b_positive": True, "context": "market", "extra": 99})
        assert f.ns5b_positive is True

    def test_border_post_context(self) -> None:
        f = flags_from_dict({"ns5b_positive": True, "context": "border_post"})
        assert f.context is CsfvContext.BORDER_POST


# ---------------------------------------------------------------------------
# File output helpers
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_writes_valid_json(self, tmp_path: Path) -> None:
        assessment = _assess(True, CsfvContext.ON_FARM)
        out = tmp_path / "csfv.json"
        write_assessment_json(assessment, out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"

    def test_negative_json_has_negative_level(self, tmp_path: Path) -> None:
        assessment = _assess(False)
        out = tmp_path / "csfv_neg.json"
        write_assessment_json(assessment, out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "NEGATIVE"

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        assessment = _assess(False)
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(assessment, out)
        assert out.exists()

    def test_file_ends_with_newline(self, tmp_path: Path) -> None:
        assessment = _assess(True, CsfvContext.BORDER_POST)
        out = tmp_path / "out.json"
        write_assessment_json(assessment, out)
        assert out.read_text().endswith("\n")


class TestWriteAssessmentCsv:
    def test_writes_csv_with_header(self, tmp_path: Path) -> None:
        assessment = _assess(True, CsfvContext.MARKET)
        out = tmp_path / "csfv.csv"
        write_assessment_csv(assessment, out)
        lines = out.read_text().splitlines()
        assert lines[0] == "field,value"

    def test_csv_has_alert_level_row(self, tmp_path: Path) -> None:
        assessment = _assess(True, CsfvContext.ON_FARM)
        out = tmp_path / "csfv.csv"
        write_assessment_csv(assessment, out)
        content = out.read_text()
        assert "alert_level" in content
        assert "CRITICAL" in content

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "a" / "b" / "csfv.csv"
        write_assessment_csv(_assess(False), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Score ordering: CRITICAL > NEGATIVE
# ---------------------------------------------------------------------------


class TestAlertLevelOrdering:
    @pytest.mark.parametrize("ctx", list(CsfvContext))
    def test_positive_above_negative(self, ctx: CsfvContext) -> None:
        assert _assess(True, ctx).alert_score > _assess(False, ctx).alert_score

    def test_negative_score_is_zero(self) -> None:
        assert _assess(False).alert_score == 0

    def test_positive_score_is_100(self) -> None:
        assert _assess(True).alert_score == 100


# ---------------------------------------------------------------------------
# Flags stored on assessment
# ---------------------------------------------------------------------------


class TestAssessmentFlagsPreserved:
    def test_flags_ns5b_preserved(self) -> None:
        assessment = _assess(True, CsfvContext.BORDER_POST)
        assert assessment.flags.ns5b_positive is True
        assert assessment.flags.context is CsfvContext.BORDER_POST

    def test_assessment_is_frozen(self) -> None:
        assessment = _assess(True)
        with pytest.raises(AttributeError):
            assessment.alert_score = 0  # type: ignore[misc]
