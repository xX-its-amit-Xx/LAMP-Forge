"""Tests for the bovine neonatal calf diarrhea LAMP panel risk assessment engine.

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

from lamp_forge.bov_enteric_risk import (
    BovEntericAlertLevel,
    BovEntericAssessment,
    BovEntericFlags,
    assess_bov_enteric_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> BovEntericFlags:
    defaults: dict[str, bool] = {
        "etec_k99": False,
        "salmonella": False,
        "crypto": False,
        "bcov": False,
        "rota_a": False,
    }
    defaults.update(kwargs)
    return BovEntericFlags(**defaults)


def _assess(**kwargs: bool) -> BovEntericAssessment:
    return assess_bov_enteric_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Alert level and score for single-target patterns
# ---------------------------------------------------------------------------


class TestSingleTargetAlertLevels:
    """Each pathogen must map to the correct alert level and score."""

    def test_all_negative_is_negative(self) -> None:
        a = _assess()
        assert a.alert_level is BovEntericAlertLevel.NEGATIVE
        assert a.alert_score == 0

    def test_etec_k99_alone_is_critical(self) -> None:
        a = _assess(etec_k99=True)
        assert a.alert_level is BovEntericAlertLevel.CRITICAL
        assert a.alert_score == 50

    def test_salmonella_alone_is_high(self) -> None:
        a = _assess(salmonella=True)
        assert a.alert_level is BovEntericAlertLevel.HIGH
        assert a.alert_score == 35

    def test_crypto_alone_is_moderate(self) -> None:
        a = _assess(crypto=True)
        assert a.alert_level is BovEntericAlertLevel.MODERATE
        assert a.alert_score == 20

    def test_bcov_alone_is_moderate(self) -> None:
        a = _assess(bcov=True)
        assert a.alert_level is BovEntericAlertLevel.MODERATE
        assert a.alert_score == 15

    def test_rota_a_alone_is_low(self) -> None:
        a = _assess(rota_a=True)
        assert a.alert_level is BovEntericAlertLevel.LOW
        assert a.alert_score == 10


# ---------------------------------------------------------------------------
# Co-detection scoring
# ---------------------------------------------------------------------------


class TestCoDetectionScores:
    """Combined flags must sum weights correctly (capped at 100)."""

    def test_etec_salmonella_score(self) -> None:
        a = _assess(etec_k99=True, salmonella=True)
        assert a.alert_score == 85
        assert a.alert_level is BovEntericAlertLevel.CRITICAL

    def test_etec_rota_a_score(self) -> None:
        a = _assess(etec_k99=True, rota_a=True)
        assert a.alert_score == 60
        assert a.alert_level is BovEntericAlertLevel.CRITICAL

    def test_salmonella_crypto_score(self) -> None:
        a = _assess(salmonella=True, crypto=True)
        assert a.alert_score == 55
        assert a.alert_level is BovEntericAlertLevel.CRITICAL

    def test_bcov_rota_a_score(self) -> None:
        a = _assess(bcov=True, rota_a=True)
        assert a.alert_score == 25
        assert a.alert_level is BovEntericAlertLevel.MODERATE

    def test_crypto_bcov_score(self) -> None:
        a = _assess(crypto=True, bcov=True)
        assert a.alert_score == 35
        assert a.alert_level is BovEntericAlertLevel.HIGH

    def test_all_positive_score_capped(self) -> None:
        a = _assess(etec_k99=True, salmonella=True, crypto=True, bcov=True, rota_a=True)
        assert a.alert_score == 100
        assert a.alert_level is BovEntericAlertLevel.CRITICAL

    def test_score_never_exceeds_100(self) -> None:
        for combo in product([False, True], repeat=5):
            etec, salm, cry, bcov, rota = combo
            a = _assess(
                etec_k99=etec,
                salmonella=salm,
                crypto=cry,
                bcov=bcov,
                rota_a=rota,
            )
            assert 0 <= a.alert_score <= 100

    def test_score_never_below_zero(self) -> None:
        a = _assess()
        assert a.alert_score == 0


# ---------------------------------------------------------------------------
# Antibiotic and zoonotic flags
# ---------------------------------------------------------------------------


class TestDecisionFlags:
    """antibiotic_indicated and zoonotic_risk must reflect the detection pattern."""

    def test_etec_k99_antibiotic_and_no_zoonotic(self) -> None:
        a = _assess(etec_k99=True)
        assert a.antibiotic_indicated is True
        assert a.zoonotic_risk is False

    def test_salmonella_antibiotic_and_zoonotic(self) -> None:
        a = _assess(salmonella=True)
        assert a.antibiotic_indicated is True
        assert a.zoonotic_risk is True

    def test_crypto_no_antibiotic_but_zoonotic(self) -> None:
        a = _assess(crypto=True)
        assert a.antibiotic_indicated is False
        assert a.zoonotic_risk is True

    def test_bcov_no_antibiotic_no_zoonotic(self) -> None:
        a = _assess(bcov=True)
        assert a.antibiotic_indicated is False
        assert a.zoonotic_risk is False

    def test_rota_a_no_antibiotic_no_zoonotic(self) -> None:
        a = _assess(rota_a=True)
        assert a.antibiotic_indicated is False
        assert a.zoonotic_risk is False

    def test_all_negative_no_flags(self) -> None:
        a = _assess()
        assert a.antibiotic_indicated is False
        assert a.zoonotic_risk is False

    def test_etec_salmonella_co_antibiotic_zoonotic(self) -> None:
        a = _assess(etec_k99=True, salmonella=True)
        assert a.antibiotic_indicated is True
        assert a.zoonotic_risk is True

    def test_salmonella_crypto_co_antibiotic_zoonotic(self) -> None:
        a = _assess(salmonella=True, crypto=True)
        assert a.antibiotic_indicated is True
        assert a.zoonotic_risk is True


# ---------------------------------------------------------------------------
# Active targets
# ---------------------------------------------------------------------------


class TestActiveTargets:
    """active_targets must list positive detections in canonical order."""

    def test_all_negative_empty(self) -> None:
        assert _assess().active_targets == ()

    def test_etec_k99_only(self) -> None:
        assert _assess(etec_k99=True).active_targets == ("ETEC_K99",)

    def test_salmonella_only(self) -> None:
        assert _assess(salmonella=True).active_targets == ("Salmonella",)

    def test_crypto_only(self) -> None:
        assert _assess(crypto=True).active_targets == ("Crypto",)

    def test_bcov_only(self) -> None:
        assert _assess(bcov=True).active_targets == ("BCoV",)

    def test_rota_a_only(self) -> None:
        assert _assess(rota_a=True).active_targets == ("RotaA",)

    def test_all_positive_canonical_order(self) -> None:
        a = _assess(etec_k99=True, salmonella=True, crypto=True, bcov=True, rota_a=True)
        assert a.active_targets == ("ETEC_K99", "Salmonella", "Crypto", "BCoV", "RotaA")

    def test_bcov_rota_a_order(self) -> None:
        a = _assess(bcov=True, rota_a=True)
        assert a.active_targets == ("BCoV", "RotaA")

    def test_salmonella_rota_a_order(self) -> None:
        a = _assess(salmonella=True, rota_a=True)
        assert a.active_targets == ("Salmonella", "RotaA")


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestInterpretationText:
    """Key phrases must appear in interpretation/action strings."""

    def test_negative_interpretation_mentions_no_pathogens(self) -> None:
        a = _assess()
        assert "negative" in a.interpretation.lower()

    def test_etec_interpretation_mentions_etec(self) -> None:
        a = _assess(etec_k99=True)
        assert "ETEC" in a.interpretation or "K99" in a.interpretation

    def test_etec_action_mentions_antibiotic(self) -> None:
        a = _assess(etec_k99=True)
        assert "antibiotic" in a.recommended_action.lower()

    def test_salmonella_interpretation_mentions_salmonella(self) -> None:
        a = _assess(salmonella=True)
        assert "Salmonella" in a.interpretation or "salmonella" in a.interpretation.lower()

    def test_salmonella_action_mentions_antibiotic(self) -> None:
        a = _assess(salmonella=True)
        assert "antibiotic" in a.recommended_action.lower()

    def test_crypto_interpretation_mentions_cryptosporidium(self) -> None:
        a = _assess(crypto=True)
        assert "cryptosporidium" in a.interpretation.lower() or "crypto" in a.interpretation.lower()

    def test_crypto_action_mentions_halofuginone(self) -> None:
        a = _assess(crypto=True)
        assert "halofuginone" in a.recommended_action.lower()

    def test_bcov_interpretation_mentions_coronavirus(self) -> None:
        a = _assess(bcov=True)
        assert "coronavirus" in a.interpretation.lower()

    def test_rota_a_interpretation_mentions_rotavirus(self) -> None:
        a = _assess(rota_a=True)
        assert "rotavirus" in a.interpretation.lower()

    def test_etec_salmonella_co_interpretation_mentions_salmonella(self) -> None:
        a = _assess(etec_k99=True, salmonella=True)
        assert "Salmonella" in a.interpretation

    def test_etec_bcov_co_interpretation_mentions_bcov(self) -> None:
        a = _assess(etec_k99=True, bcov=True)
        assert "BCoV" in a.interpretation

    def test_salmonella_crypto_co_interpretation_mentions_crypto(self) -> None:
        a = _assess(salmonella=True, crypto=True)
        assert "Cryptosporidium" in a.interpretation or "Crypto" in a.interpretation


# ---------------------------------------------------------------------------
# Serialisation: to_dict
# ---------------------------------------------------------------------------


class TestToDict:
    """to_dict must include all required keys with correct types."""

    def test_required_keys_present(self) -> None:
        d = _assess(etec_k99=True).to_dict()
        required = {
            "alert_level",
            "alert_score",
            "active_targets",
            "interpretation",
            "recommended_action",
            "antibiotic_indicated",
            "zoonotic_risk",
            "panel_flags",
        }
        assert required <= d.keys()

    def test_alert_level_is_string(self) -> None:
        d = _assess(etec_k99=True).to_dict()
        assert isinstance(d["alert_level"], str)
        assert d["alert_level"] == "CRITICAL"

    def test_panel_flags_keys(self) -> None:
        d = _assess(
            etec_k99=True, salmonella=False, crypto=True, bcov=False, rota_a=False
        ).to_dict()
        pf = d["panel_flags"]
        assert isinstance(pf, dict)
        assert pf == {
            "etec_k99": True,
            "salmonella": False,
            "crypto": True,
            "bcov": False,
            "rota_a": False,
        }

    def test_active_targets_is_list(self) -> None:
        d = _assess(etec_k99=True, rota_a=True).to_dict()
        assert isinstance(d["active_targets"], list)
        assert "ETEC_K99" in d["active_targets"]

    def test_negative_alert_level_string(self) -> None:
        d = _assess().to_dict()
        assert d["alert_level"] == "NEGATIVE"
        assert d["alert_score"] == 0


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    """flags_from_dict must parse dicts correctly and tolerate missing keys."""

    def test_all_true(self) -> None:
        f = flags_from_dict(
            {"etec_k99": True, "salmonella": True, "crypto": True, "bcov": True, "rota_a": True}
        )
        assert f == BovEntericFlags(
            etec_k99=True, salmonella=True, crypto=True, bcov=True, rota_a=True
        )

    def test_all_false_explicit(self) -> None:
        f = flags_from_dict(
            {
                "etec_k99": False,
                "salmonella": False,
                "crypto": False,
                "bcov": False,
                "rota_a": False,
            }
        )
        assert f == BovEntericFlags(
            etec_k99=False, salmonella=False, crypto=False, bcov=False, rota_a=False
        )

    def test_empty_dict_defaults_to_false(self) -> None:
        f = flags_from_dict({})
        assert f == BovEntericFlags(
            etec_k99=False, salmonella=False, crypto=False, bcov=False, rota_a=False
        )

    def test_partial_dict_etec_only(self) -> None:
        f = flags_from_dict({"etec_k99": True})
        assert f.etec_k99 is True
        assert f.salmonella is False
        assert f.crypto is False
        assert f.bcov is False
        assert f.rota_a is False

    def test_unknown_keys_ignored(self) -> None:
        f = flags_from_dict({"etec_k99": True, "unknown_pathogen": True})
        assert f.etec_k99 is True
        assert f.salmonella is False

    def test_truthy_int_converted(self) -> None:
        f = flags_from_dict({"etec_k99": 1, "salmonella": 0, "crypto": 1, "bcov": 0, "rota_a": 0})
        assert f.etec_k99 is True
        assert f.salmonella is False
        assert f.crypto is True

    def test_all_keys_present_in_expected_flags(self) -> None:
        f = flags_from_dict({"salmonella": True, "bcov": True})
        assert f.salmonella is True
        assert f.bcov is True
        assert f.etec_k99 is False
        assert f.crypto is False
        assert f.rota_a is False


# ---------------------------------------------------------------------------
# File output: write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    """write_assessment_json must produce valid, parseable JSON."""

    def test_json_round_trip(self, tmp_path: Path) -> None:
        assessment = _assess(etec_k99=True, rota_a=True)
        out = tmp_path / "out.json"
        write_assessment_json(assessment, out)
        assert out.exists()
        data = json.loads(out.read_text(encoding="utf-8"))
        assert data["alert_level"] == "CRITICAL"
        assert data["panel_flags"]["etec_k99"] is True
        assert data["panel_flags"]["rota_a"] is True
        assert data["antibiotic_indicated"] is True
        assert data["zoonotic_risk"] is False

    def test_json_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "out.json"
        write_assessment_json(_assess(), out)
        assert out.exists()

    def test_json_ends_with_newline(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        write_assessment_json(_assess(salmonella=True), out)
        content = out.read_bytes()
        assert content.endswith(b"\n")

    def test_json_zoonotic_risk_salmonella(self, tmp_path: Path) -> None:
        out = tmp_path / "zoonotic.json"
        write_assessment_json(_assess(salmonella=True), out)
        data = json.loads(out.read_text(encoding="utf-8"))
        assert data["zoonotic_risk"] is True


# ---------------------------------------------------------------------------
# File output: write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    """write_assessment_csv must produce a valid two-column key-value CSV."""

    def test_csv_has_header(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(bcov=True), out)
        with out.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert header == ["field", "value"]

    def test_csv_contains_alert_level(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(bcov=True), out)
        content = out.read_text(encoding="utf-8")
        assert "MODERATE" in content

    def test_csv_target_flags_positive_negative(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        write_assessment_csv(_assess(etec_k99=True, salmonella=False), out)
        content = out.read_text(encoding="utf-8")
        assert "positive" in content
        assert "negative" in content

    def test_csv_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "deep" / "nested" / "out.csv"
        write_assessment_csv(_assess(), out)
        assert out.exists()

    def test_csv_critical_etec(self, tmp_path: Path) -> None:
        out = tmp_path / "critical.csv"
        write_assessment_csv(_assess(etec_k99=True), out)
        content = out.read_text(encoding="utf-8")
        assert "CRITICAL" in content


# ---------------------------------------------------------------------------
# Return type and immutability
# ---------------------------------------------------------------------------


class TestReturnTypeAndImmutability:
    """assess_bov_enteric_risk must return a frozen BovEntericAssessment."""

    def test_returns_correct_type(self) -> None:
        a = _assess(etec_k99=True)
        assert isinstance(a, BovEntericAssessment)

    def test_frozen_dataclass(self) -> None:
        a = _assess(salmonella=True)
        with pytest.raises((AttributeError, TypeError)):
            a.alert_score = 0  # type: ignore[misc]

    def test_flags_field_matches_input(self) -> None:
        flags = _flags(etec_k99=False, salmonella=True, crypto=False, bcov=False, rota_a=False)
        a = assess_bov_enteric_risk(flags)
        assert a.flags is flags

    def test_exhaustive_combinations_no_exception(self) -> None:
        for combo in product([False, True], repeat=5):
            etec, salm, cry, bcov, rota = combo
            assess_bov_enteric_risk(
                _flags(etec_k99=etec, salmonella=salm, crypto=cry, bcov=bcov, rota_a=rota)
            )


# ---------------------------------------------------------------------------
# Boundary checks
# ---------------------------------------------------------------------------


class TestBoundaryConditions:
    """Verify threshold boundaries and edge cases."""

    def test_salmonella_crypto_co_triggers_critical_score(self) -> None:
        a = _assess(salmonella=True, crypto=True)
        assert a.alert_score == 55
        assert a.alert_level is BovEntericAlertLevel.CRITICAL

    def test_bcov_rota_a_stays_moderate(self) -> None:
        a = _assess(bcov=True, rota_a=True)
        assert a.alert_level is BovEntericAlertLevel.MODERATE

    def test_crypto_bcov_pushes_to_high(self) -> None:
        a = _assess(crypto=True, bcov=True)
        assert a.alert_score == 35
        assert a.alert_level is BovEntericAlertLevel.HIGH

    def test_rota_a_score_is_low_not_moderate(self) -> None:
        a = _assess(rota_a=True)
        assert a.alert_level is BovEntericAlertLevel.LOW
        assert a.alert_score == 10

    def test_bcov_score_is_moderate_not_low(self) -> None:
        a = _assess(bcov=True)
        assert a.alert_level is BovEntericAlertLevel.MODERATE
        assert a.alert_score == 15
