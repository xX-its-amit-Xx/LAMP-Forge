"""Tests for the farm-biosecurity LAMP panel risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical alert levels.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

from lamp_forge.farm_risk import (
    FarmAlertLevel,
    FarmPanelFlags,
    FarmRiskAssessment,
    assess_farm_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> FarmPanelFlags:
    defaults: dict[str, bool] = {
        "asfv": False,
        "fmdv": False,
        "aiv": False,
        "ndv": False,
        "prrsv": False,
    }
    defaults.update(kwargs)
    return FarmPanelFlags(**defaults)


def _assess(**kwargs: bool) -> FarmRiskAssessment:
    return assess_farm_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for key single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """WOAH Tier 1 notifiable pathogens must always yield CRITICAL."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is FarmAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_asfv_alone_is_critical(self) -> None:
        a = _assess(asfv=True)
        assert a.alert_level is FarmAlertLevel.CRITICAL
        assert a.alert_score >= 50

    def test_fmdv_alone_is_critical(self) -> None:
        a = _assess(fmdv=True)
        assert a.alert_level is FarmAlertLevel.CRITICAL
        assert a.alert_score >= 50

    def test_aiv_alone_is_critical(self) -> None:
        a = _assess(aiv=True)
        assert a.alert_level is FarmAlertLevel.CRITICAL
        assert a.alert_score >= 50

    def test_ndv_alone_is_high(self) -> None:
        a = _assess(ndv=True)
        assert a.alert_level is FarmAlertLevel.HIGH
        assert 25 <= a.alert_score < 50

    def test_prrsv_alone_is_moderate(self) -> None:
        a = _assess(prrsv=True)
        assert a.alert_level is FarmAlertLevel.MODERATE
        assert 20 <= a.alert_score < 25

    def test_score_increases_with_targets(self) -> None:
        s1 = _assess(prrsv=True).alert_score
        s2 = _assess(prrsv=True, ndv=True).alert_score
        s3 = _assess(prrsv=True, ndv=True, aiv=True).alert_score
        assert s1 < s2 < s3

    def test_score_capped_at_100(self) -> None:
        a = _assess(asfv=True, fmdv=True, aiv=True, ndv=True, prrsv=True)
        assert a.alert_score <= 100

    def test_score_non_negative(self) -> None:
        for combo in [{}, {"prrsv": True}, {"ndv": True}]:
            assert _assess(**combo).alert_score >= 0


# ---------------------------------------------------------------------------
# Multi-target scenarios
# ---------------------------------------------------------------------------


class TestMultiTargetPatterns:
    def test_ndv_prrsv_is_high(self) -> None:
        # NDV (25) + PRRSV (20) = 45 -> HIGH, not CRITICAL
        a = _assess(ndv=True, prrsv=True)
        assert a.alert_level is FarmAlertLevel.HIGH

    def test_aiv_prrsv_is_critical(self) -> None:
        a = _assess(aiv=True, prrsv=True)
        assert a.alert_level is FarmAlertLevel.CRITICAL

    def test_asfv_fmdv_is_critical(self) -> None:
        a = _assess(asfv=True, fmdv=True)
        assert a.alert_level is FarmAlertLevel.CRITICAL

    def test_all_five_is_critical(self) -> None:
        a = _assess(asfv=True, fmdv=True, aiv=True, ndv=True, prrsv=True)
        assert a.alert_level is FarmAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# Notifiable targets and immediate report flag
# ---------------------------------------------------------------------------


class TestNotifiableTargets:
    def test_asfv_requires_immediate_report(self) -> None:
        a = _assess(asfv=True)
        assert a.immediate_report_required is True
        assert "ASFV" in a.notifiable_targets

    def test_fmdv_requires_immediate_report(self) -> None:
        a = _assess(fmdv=True)
        assert a.immediate_report_required is True
        assert "FMDV" in a.notifiable_targets

    def test_aiv_requires_immediate_report(self) -> None:
        a = _assess(aiv=True)
        assert a.immediate_report_required is True
        assert "AIV" in a.notifiable_targets

    def test_ndv_requires_immediate_report(self) -> None:
        a = _assess(ndv=True)
        assert a.immediate_report_required is True
        assert "NDV" in a.notifiable_targets

    def test_prrsv_alone_does_not_require_immediate_report(self) -> None:
        a = _assess(prrsv=True)
        assert a.immediate_report_required is False
        assert a.notifiable_targets == ()

    def test_all_negative_no_immediate_report(self) -> None:
        a = _assess()
        assert a.immediate_report_required is False
        assert a.notifiable_targets == ()

    def test_notifiable_subset_of_active(self) -> None:
        a = _assess(asfv=True, prrsv=True)
        assert "ASFV" in a.notifiable_targets
        assert "PRRSV" not in a.notifiable_targets
        assert set(a.notifiable_targets).issubset(set(a.active_targets))


# ---------------------------------------------------------------------------
# Active targets ordering
# ---------------------------------------------------------------------------


class TestActiveTargets:
    def test_empty_when_all_negative(self) -> None:
        a = _assess()
        assert a.active_targets == ()

    def test_canonical_order_asfv_fmdv_aiv_ndv_prrsv(self) -> None:
        a = _assess(prrsv=True, asfv=True, ndv=True)
        assert a.active_targets == ("ASFV", "NDV", "PRRSV")

    def test_all_five_listed_when_all_positive(self) -> None:
        a = _assess(asfv=True, fmdv=True, aiv=True, ndv=True, prrsv=True)
        assert len(a.active_targets) == 5

    def test_only_positive_targets_listed(self) -> None:
        a = _assess(aiv=True, prrsv=True)
        assert "AIV" in a.active_targets
        assert "PRRSV" in a.active_targets
        assert "ASFV" not in a.active_targets
        assert "FMDV" not in a.active_targets
        assert "NDV" not in a.active_targets


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestTextContent:
    def test_all_negative_mentions_negative(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_asfv_interpretation_mentions_asfv(self) -> None:
        a = _assess(asfv=True)
        assert "ASFV" in a.interpretation or "african swine fever" in a.interpretation.lower()

    def test_asfv_action_mentions_quarantine(self) -> None:
        a = _assess(asfv=True)
        assert "quarantine" in a.recommended_action.lower()

    def test_fmdv_interpretation_mentions_fmdv(self) -> None:
        a = _assess(fmdv=True)
        assert "FMDV" in a.interpretation or "foot-and-mouth" in a.interpretation.lower()

    def test_aiv_interpretation_mentions_hpai_or_avian(self) -> None:
        a = _assess(aiv=True)
        text = a.interpretation.lower()
        assert "avian influenza" in text or "hpai" in text

    def test_ndv_interpretation_mentions_pathotype(self) -> None:
        a = _assess(ndv=True)
        assert "velogenic" in a.interpretation.lower() or "pathotype" in a.interpretation.lower()

    def test_prrsv_interpretation_mentions_prrsv(self) -> None:
        a = _assess(prrsv=True)
        assert "PRRSV" in a.interpretation or "prrsv" in a.interpretation.lower()

    def test_interpretation_non_empty_all_combos(self) -> None:
        keys = ["asfv", "fmdv", "aiv", "ndv", "prrsv"]
        for bits in product([False, True], repeat=5):
            flags = FarmPanelFlags(**dict(zip(keys, bits, strict=True)))
            a = assess_farm_risk(flags)
            assert a.interpretation.strip(), f"Empty interpretation for {flags}"
            assert a.recommended_action.strip(), f"Empty action for {flags}"

    def test_recommended_action_non_empty(self) -> None:
        for a in [_assess(), _assess(asfv=True), _assess(prrsv=True), _assess(ndv=True)]:
            assert len(a.recommended_action) > 20


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        d = _assess(asfv=True).to_dict()
        for key in (
            "alert_level",
            "alert_score",
            "active_targets",
            "notifiable_targets",
            "interpretation",
            "recommended_action",
            "immediate_report_required",
            "panel_flags",
        ):
            assert key in d, f"Missing key: {key}"

    def test_alert_level_is_string(self) -> None:
        d = _assess().to_dict()
        assert isinstance(d["alert_level"], str)

    def test_panel_flags_inlined(self) -> None:
        d = _assess(asfv=True, prrsv=True).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf["asfv"] is True
        assert pf["prrsv"] is True
        assert pf["fmdv"] is False

    def test_json_serialisable_critical(self) -> None:
        d = _assess(asfv=True, fmdv=True).to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "CRITICAL"
        assert recovered["immediate_report_required"] is True

    def test_json_serialisable_negative(self) -> None:
        d = _assess().to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "NEGATIVE"
        assert recovered["immediate_report_required"] is False


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_true(self) -> None:
        flags = flags_from_dict(
            {"asfv": True, "fmdv": True, "aiv": True, "ndv": True, "prrsv": True}
        )
        assert flags == FarmPanelFlags(asfv=True, fmdv=True, aiv=True, ndv=True, prrsv=True)

    def test_missing_keys_default_false(self) -> None:
        flags = flags_from_dict({"asfv": True})
        assert flags.asfv is True
        assert flags.fmdv is False
        assert flags.aiv is False

    def test_empty_dict_all_false(self) -> None:
        flags = flags_from_dict({})
        assert flags == FarmPanelFlags(asfv=False, fmdv=False, aiv=False, ndv=False, prrsv=False)

    def test_unknown_keys_ignored(self) -> None:
        flags = flags_from_dict({"asfv": True, "extra": "ignored"})
        assert flags.asfv is True


# ---------------------------------------------------------------------------
# write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "farm_assessment.json"
        write_assessment_json(_assess(asfv=True), out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        out = tmp_path / "farm_assessment.json"
        write_assessment_json(_assess(asfv=True), out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert "ASFV" in data["active_targets"]
        assert data["immediate_report_required"] is True

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
        out = tmp_path / "farm_assessment.csv"
        write_assessment_csv(_assess(asfv=True), out)
        assert out.exists()

    def test_has_header_row(self, tmp_path: Path) -> None:
        out = tmp_path / "farm_assessment.csv"
        write_assessment_csv(_assess(), out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["field", "value"]

    def test_alert_level_row_present(self, tmp_path: Path) -> None:
        out = tmp_path / "farm_assessment.csv"
        write_assessment_csv(_assess(asfv=True), out)
        rows = list(csv.DictReader(out.open()))
        level_rows = [r for r in rows if r["field"] == "alert_level"]
        assert len(level_rows) == 1
        assert level_rows[0]["value"] == "CRITICAL"

    def test_panel_flags_present_as_separate_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "farm_assessment.csv"
        write_assessment_csv(_assess(asfv=True, prrsv=True), out)
        rows = list(csv.DictReader(out.open()))
        fields = {r["field"] for r in rows}
        assert "target_asfv" in fields
        assert "target_prrsv" in fields
        assert "target_fmdv" in fields


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestFarmRiskCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["farm-risk", *args], obj={})

    def test_all_negative_exits_0(self) -> None:
        result = self._run([])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_target_table(self) -> None:
        result = self._run(["--asfv", "--prrsv"])
        output = result.output  # type: ignore[union-attr]
        assert "ASFV" in output
        assert "PRRSV" in output
        assert "[+]" in output

    def test_asfv_output_shows_critical(self) -> None:
        result = self._run(["--asfv"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_asfv_output_shows_immediate_report(self) -> None:
        result = self._run(["--asfv"])
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE" in output or "report" in output.lower()

    def test_prrsv_output_shows_moderate(self) -> None:
        result = self._run(["--prrsv"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]

    def test_ndv_output_shows_high(self) -> None:
        result = self._run(["--ndv"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_all_negative_shows_negative(self) -> None:
        result = self._run([])
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        result = self._run(["--asfv", "--out-json", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        result = self._run(["--ndv", "--prrsv", "--out-csv", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_flag(self, tmp_path: Path) -> None:
        flags_file = tmp_path / "flags.json"
        flags_file.write_text(
            '{"asfv": false, "fmdv": true, "aiv": false, "ndv": false, "prrsv": false}'
        )
        result = self._run(["--input-json", str(flags_file)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_fmdv_action_mentions_quarantine(self) -> None:
        result = self._run(["--fmdv"])
        assert "quarantine" in result.output.lower()  # type: ignore[union-attr]
