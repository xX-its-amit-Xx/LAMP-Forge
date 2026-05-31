"""Tests for the human POC LAMP panel risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical alert levels.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

from lamp_forge.poc_risk import (
    PocAlertLevel,
    PocPanelFlags,
    PocRiskAssessment,
    assess_poc_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> PocPanelFlags:
    defaults: dict[str, bool] = {
        "mtb": False,
        "cdiff": False,
        "sars2": False,
        "flu_a": False,
        "gas": False,
    }
    defaults.update(kwargs)
    return PocPanelFlags(**defaults)


def _assess(**kwargs: bool) -> PocRiskAssessment:
    return assess_poc_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for key single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """Single-target detections must map to expected alert levels."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is PocAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_mtb_alone_is_critical(self) -> None:
        a = _assess(mtb=True)
        assert a.alert_level is PocAlertLevel.CRITICAL
        assert a.alert_score >= 50

    def test_cdiff_alone_is_high(self) -> None:
        a = _assess(cdiff=True)
        assert a.alert_level is PocAlertLevel.HIGH
        assert 25 <= a.alert_score < 50

    def test_sars2_alone_is_moderate(self) -> None:
        a = _assess(sars2=True)
        assert a.alert_level is PocAlertLevel.MODERATE
        assert 15 <= a.alert_score < 25

    def test_flu_a_alone_is_moderate(self) -> None:
        a = _assess(flu_a=True)
        assert a.alert_level is PocAlertLevel.MODERATE
        assert 15 <= a.alert_score < 25

    def test_gas_alone_is_low(self) -> None:
        a = _assess(gas=True)
        assert a.alert_level is PocAlertLevel.LOW
        assert 5 <= a.alert_score < 15

    def test_score_non_negative(self) -> None:
        for combo in [{}, {"gas": True}, {"mtb": True, "cdiff": True}]:
            assert _assess(**combo).alert_score >= 0

    def test_score_capped_at_100(self) -> None:
        a = _assess(mtb=True, cdiff=True, sars2=True, flu_a=True, gas=True)
        assert a.alert_score <= 100

    def test_score_increases_with_target_count(self) -> None:
        s1 = _assess(sars2=True).alert_score
        s2 = _assess(sars2=True, flu_a=True).alert_score
        s3 = _assess(sars2=True, flu_a=True, gas=True).alert_score
        assert s1 < s2 < s3

    def test_sars2_flu_a_co_detection_is_high(self) -> None:
        a = _assess(sars2=True, flu_a=True)
        assert a.alert_level is PocAlertLevel.HIGH

    def test_mtb_cdiff_is_critical(self) -> None:
        a = _assess(mtb=True, cdiff=True)
        assert a.alert_level is PocAlertLevel.CRITICAL


# ---------------------------------------------------------------------------
# Clinical action flags
# ---------------------------------------------------------------------------


class TestClinicalActionFlags:
    def test_immediate_report_when_mtb(self) -> None:
        a = _assess(mtb=True)
        assert a.immediate_report_required is True

    def test_no_immediate_report_when_cdiff_only(self) -> None:
        a = _assess(cdiff=True)
        assert a.immediate_report_required is False

    def test_no_immediate_report_when_all_negative(self) -> None:
        a = _assess()
        assert a.immediate_report_required is False

    def test_contact_precautions_when_cdiff(self) -> None:
        a = _assess(cdiff=True)
        assert a.contact_precautions_required is True

    def test_contact_precautions_when_mtb(self) -> None:
        a = _assess(mtb=True)
        assert a.contact_precautions_required is True

    def test_no_contact_precautions_when_gas_only(self) -> None:
        a = _assess(gas=True)
        assert a.contact_precautions_required is False

    def test_antiviral_indicated_when_sars2(self) -> None:
        a = _assess(sars2=True)
        assert a.antiviral_indicated is True

    def test_antiviral_indicated_when_flu_a(self) -> None:
        a = _assess(flu_a=True)
        assert a.antiviral_indicated is True

    def test_no_antiviral_when_gas_only(self) -> None:
        a = _assess(gas=True)
        assert a.antiviral_indicated is False

    def test_no_antiviral_when_cdiff_only(self) -> None:
        a = _assess(cdiff=True)
        assert a.antiviral_indicated is False

    def test_antibiotic_indicated_when_gas(self) -> None:
        a = _assess(gas=True)
        assert a.antibiotic_indicated is True

    def test_no_antibiotic_when_sars2_only(self) -> None:
        a = _assess(sars2=True)
        assert a.antibiotic_indicated is False

    def test_antibiotic_indicated_gas_plus_sars2(self) -> None:
        a = _assess(gas=True, sars2=True)
        assert a.antibiotic_indicated is True
        assert a.antiviral_indicated is True


# ---------------------------------------------------------------------------
# Notifiable targets
# ---------------------------------------------------------------------------


class TestNotifiableTargets:
    def test_mtb_in_notifiable(self) -> None:
        a = _assess(mtb=True)
        assert "MTB" in a.notifiable_targets

    def test_cdiff_not_notifiable(self) -> None:
        a = _assess(cdiff=True)
        assert "CDIFF" not in a.notifiable_targets

    def test_gas_not_notifiable(self) -> None:
        a = _assess(gas=True)
        assert len(a.notifiable_targets) == 0

    def test_all_negative_no_notifiable(self) -> None:
        a = _assess()
        assert a.notifiable_targets == ()


# ---------------------------------------------------------------------------
# Detected targets
# ---------------------------------------------------------------------------


class TestDetectedTargets:
    def test_empty_when_all_negative(self) -> None:
        a = _assess()
        assert a.detected_targets == ()

    def test_only_positive_targets_listed(self) -> None:
        a = _assess(sars2=True, flu_a=True)
        assert "SARS2" in a.detected_targets
        assert "FLU_A" in a.detected_targets
        assert "MTB" not in a.detected_targets
        assert "GAS" not in a.detected_targets

    def test_canonical_order_mtb_cdiff_sars2_flu_a_gas(self) -> None:
        a = _assess(gas=True, sars2=True, mtb=True)
        assert a.detected_targets == ("MTB", "SARS2", "GAS")

    def test_all_five_detected_when_all_positive(self) -> None:
        a = _assess(mtb=True, cdiff=True, sars2=True, flu_a=True, gas=True)
        assert len(a.detected_targets) == 5


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestTextContent:
    def test_all_negative_mentions_no_pathogen(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower() or "no pathogen" in a.interpretation.lower()

    def test_mtb_interpretation_mentions_tuberculosis(self) -> None:
        a = _assess(mtb=True)
        text = a.interpretation.lower()
        assert "tuberculosis" in text or "mtb" in text or "tb" in text

    def test_mtb_action_mentions_isolation(self) -> None:
        a = _assess(mtb=True)
        assert "isolation" in a.recommended_action.lower()

    def test_cdiff_interpretation_mentions_cdiff_or_difficile(self) -> None:
        a = _assess(cdiff=True)
        text = a.interpretation.lower()
        assert "difficile" in text or "cdiff" in text or "cdi" in text

    def test_cdiff_action_mentions_contact_precautions(self) -> None:
        a = _assess(cdiff=True)
        assert "contact" in a.recommended_action.lower()

    def test_sars2_interpretation_mentions_sars_or_covid(self) -> None:
        a = _assess(sars2=True)
        text = a.interpretation.lower()
        assert "sars" in text or "covid" in text or "coronavirus" in text

    def test_flu_interpretation_mentions_influenza(self) -> None:
        a = _assess(flu_a=True)
        assert "influenza" in a.interpretation.lower()

    def test_gas_interpretation_mentions_streptococcus_or_gas(self) -> None:
        a = _assess(gas=True)
        text = a.interpretation.lower()
        assert "streptococcus" in text or "gas" in text or "strep" in text

    def test_gas_action_mentions_antibiotic(self) -> None:
        a = _assess(gas=True)
        assert (
            "antibiotic" in a.recommended_action.lower()
            or "penicillin" in a.recommended_action.lower()
        )

    def test_interpretation_non_empty_for_all_32_combos(self) -> None:
        keys = ["mtb", "cdiff", "sars2", "flu_a", "gas"]
        for bits in product([False, True], repeat=5):
            flags = PocPanelFlags(**dict(zip(keys, bits, strict=True)))
            a = assess_poc_risk(flags)
            assert a.interpretation.strip(), f"Empty interpretation for {flags}"
            assert a.recommended_action.strip(), f"Empty action for {flags}"

    def test_recommended_action_non_empty(self) -> None:
        for a in [
            _assess(),
            _assess(mtb=True),
            _assess(cdiff=True),
            _assess(sars2=True, flu_a=True, gas=True),
        ]:
            assert len(a.recommended_action) > 20


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        d = _assess(sars2=True).to_dict()
        for key in (
            "alert_level",
            "alert_score",
            "detected_targets",
            "notifiable_targets",
            "immediate_report_required",
            "contact_precautions_required",
            "antiviral_indicated",
            "antibiotic_indicated",
            "interpretation",
            "recommended_action",
            "panel_flags",
        ):
            assert key in d, f"Missing key: {key}"

    def test_alert_level_is_string(self) -> None:
        d = _assess().to_dict()
        assert isinstance(d["alert_level"], str)

    def test_panel_flags_inlined(self) -> None:
        d = _assess(sars2=True, gas=True).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf["sars2"] is True
        assert pf["gas"] is True
        assert pf["mtb"] is False

    def test_json_serialisable(self) -> None:
        d = _assess(mtb=True, cdiff=True).to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "CRITICAL"

    def test_detected_targets_is_list(self) -> None:
        d = _assess(flu_a=True).to_dict()
        assert isinstance(d["detected_targets"], list)
        assert "FLU_A" in d["detected_targets"]


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_true(self) -> None:
        flags = flags_from_dict(
            {"mtb": True, "cdiff": True, "sars2": True, "flu_a": True, "gas": True}
        )
        assert flags == PocPanelFlags(mtb=True, cdiff=True, sars2=True, flu_a=True, gas=True)

    def test_missing_keys_default_false(self) -> None:
        flags = flags_from_dict({"mtb": True})
        assert flags.mtb is True
        assert flags.cdiff is False
        assert flags.sars2 is False
        assert flags.flu_a is False
        assert flags.gas is False

    def test_empty_dict_all_false(self) -> None:
        flags = flags_from_dict({})
        assert flags == PocPanelFlags(mtb=False, cdiff=False, sars2=False, flu_a=False, gas=False)

    def test_unknown_keys_ignored(self) -> None:
        flags = flags_from_dict({"sars2": True, "extra_key": "ignored"})
        assert flags.sars2 is True


# ---------------------------------------------------------------------------
# write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        a = _assess(cdiff=True)
        out = tmp_path / "assessment.json"
        write_assessment_json(a, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        a = _assess(mtb=True)
        out = tmp_path / "assessment.json"
        write_assessment_json(a, out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert "MTB" in data["detected_targets"]

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "poc_result.json"
        write_assessment_json(_assess(), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(gas=True), out)
        assert out.exists()

    def test_has_header_row(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(), out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["field", "value"]

    def test_alert_level_row_present(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(cdiff=True), out)
        rows = list(csv.DictReader(out.open()))
        level_rows = [r for r in rows if r["field"] == "alert_level"]
        assert len(level_rows) == 1
        assert level_rows[0]["value"] == "HIGH"

    def test_panel_flags_present_as_separate_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(sars2=True), out)
        rows = list(csv.DictReader(out.open()))
        fields = {r["field"] for r in rows}
        assert "target_sars2" in fields
        assert "target_mtb" in fields
        assert "target_gas" in fields


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestPocRiskCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["poc-risk", *args], obj={})

    def test_all_negative_exits_0(self) -> None:
        result = self._run([])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_target_table(self) -> None:
        result = self._run(["--cdiff", "--gas"])
        output = result.output  # type: ignore[union-attr]
        assert "CDIFF" in output
        assert "GAS" in output
        assert "[+]" in output

    def test_mtb_shows_critical(self) -> None:
        result = self._run(["--mtb"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_cdiff_shows_high(self) -> None:
        result = self._run(["--cdiff"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_sars2_shows_moderate(self) -> None:
        result = self._run(["--sars2"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]

    def test_gas_alone_shows_low(self) -> None:
        result = self._run(["--gas"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "LOW" in result.output  # type: ignore[union-attr]

    def test_negative_shows_negative(self) -> None:
        result = self._run([])
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_mtb_shows_immediate_report_message(self) -> None:
        result = self._run(["--mtb"])
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE REPORT" in output or "immediate" in output.lower()

    def test_cdiff_shows_contact_precautions(self) -> None:
        result = self._run(["--cdiff"])
        output = result.output  # type: ignore[union-attr]
        assert "contact" in output.lower() or "precautions" in output.lower()

    def test_sars2_shows_antiviral_message(self) -> None:
        result = self._run(["--sars2"])
        output = result.output  # type: ignore[union-attr]
        assert "antiviral" in output.lower()

    def test_gas_shows_antibiotic_message(self) -> None:
        result = self._run(["--gas"])
        output = result.output  # type: ignore[union-attr]
        assert "antibiotic" in output.lower()

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        result = self._run(["--cdiff", "--out-json", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["alert_level"] == "HIGH"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        result = self._run(["--sars2", "--out-csv", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_flag(self, tmp_path: Path) -> None:
        flags_file = tmp_path / "flags.json"
        flags_file.write_text(
            '{"mtb": false, "cdiff": true, "sars2": false, "flu_a": false, "gas": false}'
        )
        result = self._run(["--input-json", str(flags_file)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_flu_a_flag(self) -> None:
        result = self._run(["--flu-a"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]
