"""Tests for lamp_forge.poc_respiratory_risk -- respiratory differential LAMP panel.

All tests are pure Python; no external binaries (MAFFT, BLAST) or network
calls are required.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from lamp_forge.poc_respiratory_risk import (
    PocRespAlertLevel,
    PocRespAssessment,
    PocRespFlags,
    assess_poc_resp_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(
    hrsv: bool = False,
    iav: bool = False,
    sars2: bool = False,
) -> PocRespFlags:
    return PocRespFlags(hrsv=hrsv, iav=iav, sars2=sars2)


def _assess(
    hrsv: bool = False,
    iav: bool = False,
    sars2: bool = False,
) -> PocRespAssessment:
    return assess_poc_resp_risk(_flags(hrsv=hrsv, iav=iav, sars2=sars2))


# ---------------------------------------------------------------------------
# Alert levels for each pattern
# ---------------------------------------------------------------------------


class TestAlertLevels:
    def test_all_negative_is_negative(self) -> None:
        assert _assess().alert_level is PocRespAlertLevel.NEGATIVE

    def test_hrsv_alone_is_low(self) -> None:
        assert _assess(hrsv=True).alert_level is PocRespAlertLevel.LOW

    def test_iav_alone_is_moderate(self) -> None:
        assert _assess(iav=True).alert_level is PocRespAlertLevel.MODERATE

    def test_sars2_alone_is_moderate(self) -> None:
        assert _assess(sars2=True).alert_level is PocRespAlertLevel.MODERATE

    def test_hrsv_iav_is_high(self) -> None:
        assert _assess(hrsv=True, iav=True).alert_level is PocRespAlertLevel.HIGH

    def test_hrsv_sars2_is_high(self) -> None:
        assert _assess(hrsv=True, sars2=True).alert_level is PocRespAlertLevel.HIGH

    def test_iav_sars2_is_high(self) -> None:
        assert _assess(iav=True, sars2=True).alert_level is PocRespAlertLevel.HIGH

    def test_triple_is_high(self) -> None:
        assert _assess(hrsv=True, iav=True, sars2=True).alert_level is PocRespAlertLevel.HIGH


# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------


class TestScores:
    def test_negative_score_is_zero(self) -> None:
        assert _assess().alert_score == 0

    def test_hrsv_score_is_10(self) -> None:
        assert _assess(hrsv=True).alert_score == 10

    def test_iav_score_is_20(self) -> None:
        assert _assess(iav=True).alert_score == 20

    def test_sars2_score_is_25(self) -> None:
        assert _assess(sars2=True).alert_score == 25

    def test_hrsv_iav_score_is_30(self) -> None:
        assert _assess(hrsv=True, iav=True).alert_score == 30

    def test_hrsv_sars2_score_is_35(self) -> None:
        assert _assess(hrsv=True, sars2=True).alert_score == 35

    def test_iav_sars2_score_is_45(self) -> None:
        assert _assess(iav=True, sars2=True).alert_score == 45

    def test_triple_score_is_55(self) -> None:
        assert _assess(hrsv=True, iav=True, sars2=True).alert_score == 55

    def test_score_capped_at_100(self) -> None:
        a = _assess(hrsv=True, iav=True, sars2=True)
        assert a.alert_score <= 100


# ---------------------------------------------------------------------------
# Clinical flags
# ---------------------------------------------------------------------------


class TestClinicalFlags:
    def test_antiviral_indicated_for_iav(self) -> None:
        assert _assess(iav=True).antiviral_indicated is True

    def test_antiviral_indicated_for_sars2(self) -> None:
        assert _assess(sars2=True).antiviral_indicated is True

    def test_antiviral_not_indicated_for_hrsv_alone(self) -> None:
        assert _assess(hrsv=True).antiviral_indicated is False

    def test_antiviral_not_indicated_for_negative(self) -> None:
        assert _assess().antiviral_indicated is False

    def test_rsv_indicated_for_hrsv(self) -> None:
        assert _assess(hrsv=True).rsv_indicated is True

    def test_rsv_not_indicated_for_iav_alone(self) -> None:
        assert _assess(iav=True).rsv_indicated is False

    def test_dual_infection_for_two_positives(self) -> None:
        assert _assess(hrsv=True, iav=True).dual_infection is True

    def test_dual_infection_for_triple(self) -> None:
        assert _assess(hrsv=True, iav=True, sars2=True).dual_infection is True

    def test_not_dual_infection_for_single(self) -> None:
        assert _assess(sars2=True).dual_infection is False

    def test_not_dual_infection_for_negative(self) -> None:
        assert _assess().dual_infection is False


# ---------------------------------------------------------------------------
# Detected targets
# ---------------------------------------------------------------------------


class TestDetectedTargets:
    def test_negative_has_empty_detected(self) -> None:
        assert _assess().detected_targets == ()

    def test_hrsv_alone_detected(self) -> None:
        assert _assess(hrsv=True).detected_targets == ("hRSV",)

    def test_iav_alone_detected(self) -> None:
        assert _assess(iav=True).detected_targets == ("IAV",)

    def test_sars2_alone_detected(self) -> None:
        assert _assess(sars2=True).detected_targets == ("SARS2",)

    def test_triple_detected_in_canonical_order(self) -> None:
        a = _assess(hrsv=True, iav=True, sars2=True)
        assert a.detected_targets == ("hRSV", "IAV", "SARS2")


# ---------------------------------------------------------------------------
# Text content checks
# ---------------------------------------------------------------------------


class TestTextContent:
    def test_negative_interpretation_mentions_negative(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_iav_alone_mentions_oseltamivir(self) -> None:
        a = _assess(iav=True)
        assert "oseltamivir" in a.recommended_action.lower()

    def test_sars2_alone_mentions_nirmatrelvir(self) -> None:
        a = _assess(sars2=True)
        assert "nirmatrelvir" in a.recommended_action.lower()

    def test_hrsv_alone_mentions_supportive(self) -> None:
        a = _assess(hrsv=True)
        assert "supportive" in a.recommended_action.lower()

    def test_hrsv_mentions_nirsevimab(self) -> None:
        a = _assess(hrsv=True)
        assert "nirsevimab" in a.recommended_action.lower()

    def test_iav_sars2_mentions_isolation(self) -> None:
        a = _assess(iav=True, sars2=True)
        assert "isolation" in a.recommended_action.lower()

    def test_triple_interpretation_not_empty(self) -> None:
        a = _assess(hrsv=True, iav=True, sars2=True)
        assert len(a.interpretation) > 20
        assert len(a.recommended_action) > 20


# ---------------------------------------------------------------------------
# Frozen dataclass
# ---------------------------------------------------------------------------


class TestFrozenDataclass:
    def test_flags_are_frozen(self) -> None:
        f = _flags(hrsv=True)
        with pytest.raises(AttributeError):
            f.hrsv = False  # type: ignore[misc]

    def test_assessment_is_frozen(self) -> None:
        a = _assess(sars2=True)
        with pytest.raises(AttributeError):
            a.alert_score = 0  # type: ignore[misc]


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_present(self) -> None:
        f = flags_from_dict({"hrsv": True, "iav": False, "sars2": True})
        assert f.hrsv is True
        assert f.iav is False
        assert f.sars2 is True

    def test_missing_keys_default_false(self) -> None:
        f = flags_from_dict({})
        assert f.hrsv is False
        assert f.iav is False
        assert f.sars2 is False

    def test_truthy_int_values(self) -> None:
        f = flags_from_dict({"hrsv": 1, "iav": 0, "sars2": 1})
        assert f.hrsv is True
        assert f.iav is False
        assert f.sars2 is True


# ---------------------------------------------------------------------------
# to_dict / round-trip
# ---------------------------------------------------------------------------


class TestToDict:
    def test_to_dict_keys(self) -> None:
        a = _assess(iav=True, sars2=True)
        d = a.to_dict()
        assert "alert_level" in d
        assert "alert_score" in d
        assert "detected_targets" in d
        assert "antiviral_indicated" in d
        assert "rsv_indicated" in d
        assert "dual_infection" in d
        assert "panel_flags" in d

    def test_to_dict_panel_flags_match(self) -> None:
        a = _assess(hrsv=True, sars2=True)
        d = a.to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf["hrsv"] is True
        assert pf["sars2"] is True
        assert pf["iav"] is False

    def test_round_trip_via_flags_from_dict(self) -> None:
        original = _assess(iav=True)
        d = original.to_dict()
        recovered_flags = flags_from_dict(d["panel_flags"])  # type: ignore[arg-type]
        recovered = assess_poc_resp_risk(recovered_flags)
        assert recovered.alert_level is original.alert_level
        assert recovered.alert_score == original.alert_score


# ---------------------------------------------------------------------------
# JSON / CSV I/O
# ---------------------------------------------------------------------------


class TestJsonOutput:
    def test_write_and_read_json(self, tmp_path: Path) -> None:
        a = _assess(hrsv=True, sars2=True)
        out = tmp_path / "resp_assessment.json"
        write_assessment_json(a, out)
        assert out.exists()
        loaded = json.loads(out.read_text())
        assert loaded["alert_level"] == PocRespAlertLevel.HIGH.value
        assert loaded["panel_flags"]["hrsv"] is True

    def test_write_json_creates_parent_dir(self, tmp_path: Path) -> None:
        a = _assess(iav=True)
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(a, out)
        assert out.exists()


class TestCsvOutput:
    def test_write_csv_has_header(self, tmp_path: Path) -> None:
        a = _assess(sars2=True)
        out = tmp_path / "resp.csv"
        write_assessment_csv(a, out)
        lines = out.read_text().splitlines()
        assert lines[0] == "key,value"

    def test_write_csv_contains_alert_level(self, tmp_path: Path) -> None:
        a = _assess(iav=True)
        out = tmp_path / "resp.csv"
        write_assessment_csv(a, out)
        content = out.read_text()
        assert "alert_level" in content
        assert "MODERATE" in content

    def test_write_csv_creates_parent_dir(self, tmp_path: Path) -> None:
        a = _assess(hrsv=True)
        out = tmp_path / "nested" / "resp.csv"
        write_assessment_csv(a, out)
        assert out.exists()
