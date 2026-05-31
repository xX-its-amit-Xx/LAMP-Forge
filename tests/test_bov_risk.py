"""Tests for the bovine respiratory LAMP panel risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical alert levels.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

from lamp_forge.bov_risk import (
    BovAlertLevel,
    BovPanelFlags,
    BovRiskAssessment,
    assess_bov_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> BovPanelFlags:
    defaults: dict[str, bool] = {
        "brsv": False,
        "bcov": False,
        "bvdv": False,
        "ibr": False,
        "mhae": False,
    }
    defaults.update(kwargs)
    return BovPanelFlags(**defaults)


def _assess(**kwargs: bool) -> BovRiskAssessment:
    return assess_bov_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is BovAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_brsv_alone_is_high(self) -> None:
        a = _assess(brsv=True)
        assert a.alert_level is BovAlertLevel.HIGH
        assert a.alert_score == 35

    def test_bcov_alone_is_moderate(self) -> None:
        a = _assess(bcov=True)
        assert a.alert_level is BovAlertLevel.MODERATE
        assert a.alert_score == 20

    def test_bvdv_alone_is_high(self) -> None:
        a = _assess(bvdv=True)
        assert a.alert_level is BovAlertLevel.HIGH
        assert a.alert_score == 30

    def test_ibr_alone_is_moderate(self) -> None:
        # IBR weight=25; boundary HIGH>=30, MODERATE>=20; so 25 -> MODERATE
        a = _assess(ibr=True)
        assert a.alert_level is BovAlertLevel.MODERATE
        assert a.alert_score == 25

    def test_mhae_alone_is_moderate(self) -> None:
        a = _assess(mhae=True)
        assert a.alert_level is BovAlertLevel.MODERATE
        assert a.alert_score == 20

    def test_score_increases_with_targets(self) -> None:
        s1 = _assess(bcov=True).alert_score
        s2 = _assess(bcov=True, ibr=True).alert_score
        s3 = _assess(bcov=True, ibr=True, brsv=True).alert_score
        assert s1 < s2 < s3

    def test_score_capped_at_100(self) -> None:
        a = _assess(brsv=True, bcov=True, bvdv=True, ibr=True, mhae=True)
        assert a.alert_score <= 100

    def test_score_non_negative(self) -> None:
        for kwargs in [{}, {"brsv": True}, {"mhae": True}]:
            assert _assess(**kwargs).alert_score >= 0


# ---------------------------------------------------------------------------
# Multi-target alert levels
# ---------------------------------------------------------------------------


class TestMultiTargetAlertLevels:
    def test_brsv_mhae_is_critical(self) -> None:
        # 35 + 20 = 55 -> CRITICAL
        a = _assess(brsv=True, mhae=True)
        assert a.alert_level is BovAlertLevel.CRITICAL

    def test_bvdv_mhae_is_critical(self) -> None:
        # 30 + 20 = 50 -> CRITICAL
        a = _assess(bvdv=True, mhae=True)
        assert a.alert_level is BovAlertLevel.CRITICAL

    def test_bcov_ibr_is_high(self) -> None:
        # 20 + 25 = 45 -> HIGH (not CRITICAL)
        a = _assess(bcov=True, ibr=True)
        assert a.alert_level is BovAlertLevel.HIGH
        assert a.alert_score == 45

    def test_all_five_is_critical(self) -> None:
        a = _assess(brsv=True, bcov=True, bvdv=True, ibr=True, mhae=True)
        assert a.alert_level is BovAlertLevel.CRITICAL

    def test_ibr_mhae_is_high(self) -> None:
        # IBR(25) + MHAE(20) = 45 -> HIGH (>=30, <50); no viral co-infection with BRSV/BCOV/BVDV
        a = _assess(ibr=True, mhae=True)
        assert a.alert_level is BovAlertLevel.HIGH
        assert a.alert_score == 45


# ---------------------------------------------------------------------------
# Bacterial co-infection flag
# ---------------------------------------------------------------------------


class TestBacterialCoinfection:
    def test_brsv_mhae_is_coinfection(self) -> None:
        a = _assess(brsv=True, mhae=True)
        assert a.bacterial_coinfection is True

    def test_bcov_mhae_is_coinfection(self) -> None:
        a = _assess(bcov=True, mhae=True)
        assert a.bacterial_coinfection is True

    def test_bvdv_mhae_is_coinfection(self) -> None:
        a = _assess(bvdv=True, mhae=True)
        assert a.bacterial_coinfection is True

    def test_ibr_mhae_is_coinfection(self) -> None:
        a = _assess(ibr=True, mhae=True)
        assert a.bacterial_coinfection is True

    def test_mhae_alone_not_coinfection(self) -> None:
        a = _assess(mhae=True)
        assert a.bacterial_coinfection is False

    def test_brsv_alone_not_coinfection(self) -> None:
        a = _assess(brsv=True)
        assert a.bacterial_coinfection is False

    def test_all_negative_not_coinfection(self) -> None:
        a = _assess()
        assert a.bacterial_coinfection is False


# ---------------------------------------------------------------------------
# IBR mandatory notification flag
# ---------------------------------------------------------------------------


class TestIbrMandatoryNotification:
    def test_ibr_sets_notification(self) -> None:
        a = _assess(ibr=True)
        assert a.ibr_mandatory_notification is True

    def test_ibr_not_set_when_absent(self) -> None:
        a = _assess(brsv=True, mhae=True)
        assert a.ibr_mandatory_notification is False

    def test_all_negative_no_notification(self) -> None:
        a = _assess()
        assert a.ibr_mandatory_notification is False

    def test_ibr_with_others_still_sets_notification(self) -> None:
        a = _assess(ibr=True, brsv=True, mhae=True)
        assert a.ibr_mandatory_notification is True


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    def test_empty_when_all_negative(self) -> None:
        a = _assess()
        assert a.active_targets == ()

    def test_canonical_order_brsv_bcov_bvdv_ibr_mhae(self) -> None:
        a = _assess(mhae=True, brsv=True, ibr=True)
        assert a.active_targets == ("BRSV", "IBR", "MHAE")

    def test_all_five_listed_when_all_positive(self) -> None:
        a = _assess(brsv=True, bcov=True, bvdv=True, ibr=True, mhae=True)
        assert len(a.active_targets) == 5

    def test_only_positive_targets_listed(self) -> None:
        a = _assess(bcov=True, bvdv=True)
        assert "BCOV" in a.active_targets
        assert "BVDV" in a.active_targets
        assert "BRSV" not in a.active_targets
        assert "IBR" not in a.active_targets
        assert "MHAE" not in a.active_targets


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestTextContent:
    def test_all_negative_mentions_negative(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_brsv_mhae_coinfection_text(self) -> None:
        a = _assess(brsv=True, mhae=True)
        text = a.interpretation.lower()
        assert "co-infection" in text or "mannheimia" in text

    def test_brsv_mhae_action_mentions_antimicrobial(self) -> None:
        a = _assess(brsv=True, mhae=True)
        assert "antimicrobial" in a.recommended_action.lower()

    def test_bvdv_mentions_pi_animals(self) -> None:
        a = _assess(bvdv=True)
        text = a.recommended_action.lower()
        assert "pi" in text or "persistently infected" in text

    def test_ibr_mentions_eu_or_regulation(self) -> None:
        a = _assess(ibr=True)
        text = (a.interpretation + a.recommended_action).lower()
        assert "eu" in text or "mandatory" in text or "eradication" in text

    def test_mhae_alone_mentions_antimicrobial(self) -> None:
        a = _assess(mhae=True)
        assert "antimicrobial" in a.recommended_action.lower()

    def test_all_combos_have_non_empty_text(self) -> None:
        keys = ["brsv", "bcov", "bvdv", "ibr", "mhae"]
        for bits in product([False, True], repeat=5):
            flags = BovPanelFlags(**dict(zip(keys, bits, strict=True)))
            a = assess_bov_risk(flags)
            assert a.interpretation.strip(), f"Empty interpretation for {flags}"
            assert a.recommended_action.strip(), f"Empty action for {flags}"

    def test_recommended_action_non_empty(self) -> None:
        for a in [_assess(), _assess(brsv=True), _assess(mhae=True), _assess(ibr=True)]:
            assert len(a.recommended_action) > 20


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        d = _assess(brsv=True).to_dict()
        for key in (
            "alert_level",
            "alert_score",
            "active_targets",
            "ibr_mandatory_notification",
            "bacterial_coinfection",
            "interpretation",
            "recommended_action",
            "panel_flags",
        ):
            assert key in d, f"Missing key: {key}"

    def test_alert_level_is_string(self) -> None:
        d = _assess().to_dict()
        assert isinstance(d["alert_level"], str)

    def test_panel_flags_inlined(self) -> None:
        d = _assess(brsv=True, mhae=True).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf["brsv"] is True
        assert pf["mhae"] is True
        assert pf["bcov"] is False

    def test_json_serialisable_critical(self) -> None:
        d = _assess(brsv=True, mhae=True).to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "CRITICAL"
        assert recovered["bacterial_coinfection"] is True

    def test_json_serialisable_negative(self) -> None:
        d = _assess().to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["alert_level"] == "NEGATIVE"
        assert recovered["bacterial_coinfection"] is False


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_true(self) -> None:
        flags = flags_from_dict(
            {"brsv": True, "bcov": True, "bvdv": True, "ibr": True, "mhae": True}
        )
        assert flags == BovPanelFlags(brsv=True, bcov=True, bvdv=True, ibr=True, mhae=True)

    def test_missing_keys_default_false(self) -> None:
        flags = flags_from_dict({"brsv": True})
        assert flags.brsv is True
        assert flags.bcov is False
        assert flags.bvdv is False

    def test_empty_dict_all_false(self) -> None:
        flags = flags_from_dict({})
        assert flags == BovPanelFlags(brsv=False, bcov=False, bvdv=False, ibr=False, mhae=False)

    def test_unknown_keys_ignored(self) -> None:
        flags = flags_from_dict({"brsv": True, "extra_key": "ignored"})
        assert flags.brsv is True


# ---------------------------------------------------------------------------
# write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "bov_assessment.json"
        write_assessment_json(_assess(brsv=True, mhae=True), out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        out = tmp_path / "bov_assessment.json"
        write_assessment_json(_assess(brsv=True, mhae=True), out)
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert "BRSV" in data["active_targets"]
        assert data["bacterial_coinfection"] is True

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
        out = tmp_path / "bov_assessment.csv"
        write_assessment_csv(_assess(brsv=True), out)
        assert out.exists()

    def test_has_header_row(self, tmp_path: Path) -> None:
        out = tmp_path / "bov_assessment.csv"
        write_assessment_csv(_assess(), out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["field", "value"]

    def test_alert_level_row_present(self, tmp_path: Path) -> None:
        out = tmp_path / "bov_assessment.csv"
        write_assessment_csv(_assess(brsv=True, mhae=True), out)
        rows = list(csv.DictReader(out.open()))
        level_rows = [r for r in rows if r["field"] == "alert_level"]
        assert len(level_rows) == 1
        assert level_rows[0]["value"] == "CRITICAL"

    def test_panel_flags_present_as_separate_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "bov_assessment.csv"
        write_assessment_csv(_assess(brsv=True, mhae=True), out)
        rows = list(csv.DictReader(out.open()))
        fields = {r["field"] for r in rows}
        assert "target_brsv" in fields
        assert "target_mhae" in fields
        assert "target_bcov" in fields


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestBovRiskCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["bov-risk", *args], obj={})

    def test_all_negative_exits_0(self) -> None:
        result = self._run([])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_target_table(self) -> None:
        result = self._run(["--brsv", "--mhae"])
        output = result.output  # type: ignore[union-attr]
        assert "BRSV" in output
        assert "MHAE" in output
        assert "[+]" in output

    def test_brsv_mhae_shows_critical(self) -> None:
        result = self._run(["--brsv", "--mhae"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_brsv_mhae_warns_coinfection(self) -> None:
        result = self._run(["--brsv", "--mhae"])
        output = result.output  # type: ignore[union-attr]
        assert "co-infection" in output.lower() or "coinfection" in output.lower()

    def test_ibr_warns_notification(self) -> None:
        result = self._run(["--ibr"])
        output = result.output  # type: ignore[union-attr]
        assert "notification" in output.lower() or "mandatory" in output.lower()

    def test_bvdv_shows_high(self) -> None:
        result = self._run(["--bvdv"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_bcov_shows_moderate(self) -> None:
        result = self._run(["--bcov"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]

    def test_all_negative_shows_negative(self) -> None:
        result = self._run([])
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        result = self._run(["--brsv", "--mhae", "--out-json", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["alert_level"] == "CRITICAL"
        assert data["bacterial_coinfection"] is True

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        result = self._run(["--bvdv", "--out-csv", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_flag(self, tmp_path: Path) -> None:
        flags_file = tmp_path / "flags.json"
        flags_file.write_text(
            '{"brsv": false, "bcov": false, "bvdv": true, "ibr": false, "mhae": false}'
        )
        result = self._run(["--input-json", str(flags_file)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_brsv_action_mentions_vaccination_or_biosecurity(self) -> None:
        result = self._run(["--brsv"])
        output = result.output.lower()  # type: ignore[union-attr]
        assert "vaccin" in output or "biosecurity" in output or "ventilation" in output
