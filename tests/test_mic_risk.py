"""Tests for the oilfield MIC risk assessment engine.

All tests are pure Python (no external binaries required).  The scoring model
is deterministic, so tests use exact integer comparisons or assertions on
categorical risk levels.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from lamp_forge.mic_risk import (
    GuildFlags,
    MICRiskAssessment,
    MICRiskLevel,
    assess_mic_risk,
    flags_from_dict,
    write_assessment_csv,
    write_assessment_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _flags(**kwargs: bool) -> GuildFlags:
    defaults: dict[str, bool] = {
        "srb": False,
        "mcr": False,
        "irb": False,
        "apb": False,
        "nrb": False,
    }
    defaults.update(kwargs)
    return GuildFlags(**defaults)


def _assess(**kwargs: bool) -> MICRiskAssessment:
    return assess_mic_risk(_flags(**kwargs))


# ---------------------------------------------------------------------------
# Score and risk level for key patterns from cookbook tables
# ---------------------------------------------------------------------------


class TestRiskScorePatterns:
    """Verify score and risk level for the patterns described in cookbook Recipes 15-17."""

    def test_all_negative_is_minimal(self) -> None:
        a = _assess()
        assert a.risk_level is MICRiskLevel.MINIMAL
        assert a.risk_score == 0

    def test_nrb_only_is_low(self) -> None:
        a = _assess(nrb=True)
        assert a.risk_level is MICRiskLevel.LOW
        assert a.risk_score >= 8

    def test_srb_only_is_moderate(self) -> None:
        a = _assess(srb=True)
        assert a.risk_level is MICRiskLevel.MODERATE
        assert 30 <= a.risk_score < 60

    def test_irb_only_is_low(self) -> None:
        a = _assess(irb=True)
        assert a.risk_level is MICRiskLevel.LOW

    def test_apb_only_is_low(self) -> None:
        a = _assess(apb=True)
        assert a.risk_level is MICRiskLevel.LOW

    def test_mcr_only_is_low(self) -> None:
        a = _assess(mcr=True)
        assert a.risk_level is MICRiskLevel.LOW

    def test_srb_mcr_is_moderate(self) -> None:
        # Recipe 17: SRB+MCR = sulfidogenesis + gas souring, standard souring protocol
        a = _assess(srb=True, mcr=True)
        assert a.risk_level is MICRiskLevel.MODERATE

    def test_srb_irb_is_high(self) -> None:
        # Recipe 15: SRB+IRB = FeS scale formation, HIGH risk
        a = _assess(srb=True, irb=True)
        assert a.risk_level is MICRiskLevel.HIGH

    def test_srb_irb_apb_is_high(self) -> None:
        # Recipe 16: SRB+IRB+APB = syntrophic loop + pitting, HIGH
        a = _assess(srb=True, irb=True, apb=True)
        assert a.risk_level is MICRiskLevel.HIGH

    def test_srb_mcr_irb_is_high(self) -> None:
        # Recipe 15: SRB+MCR+IRB = high MIC + souring
        a = _assess(srb=True, mcr=True, irb=True)
        assert a.risk_level is MICRiskLevel.HIGH

    def test_all_four_corrosion_guilds_is_critical(self) -> None:
        # Recipe 16: SRB+MCR+IRB+APB = maximum MIC risk
        a = _assess(srb=True, mcr=True, irb=True, apb=True)
        assert a.risk_level is MICRiskLevel.CRITICAL
        assert a.risk_score >= 90

    def test_all_five_guilds_is_critical(self) -> None:
        # Recipe 17: all positive (NRB not suppressing SRB at this guild level)
        a = _assess(srb=True, mcr=True, irb=True, apb=True, nrb=True)
        assert a.risk_level is MICRiskLevel.CRITICAL

    def test_srb_nrb_is_low_to_moderate(self) -> None:
        # Recipe 17: SRB+NRB = treatment under-dosed, NRB not suppressing
        a = _assess(srb=True, nrb=True)
        assert a.risk_level in (MICRiskLevel.LOW, MICRiskLevel.MODERATE)

    def test_score_increases_with_guild_count(self) -> None:
        s1 = _assess(srb=True).risk_score
        s2 = _assess(srb=True, irb=True).risk_score
        s3 = _assess(srb=True, irb=True, apb=True).risk_score
        assert s1 < s2 < s3

    def test_score_capped_at_100(self) -> None:
        a = _assess(srb=True, mcr=True, irb=True, apb=True, nrb=False)
        assert a.risk_score <= 100

    def test_score_non_negative(self) -> None:
        for combo in [
            {"nrb": True},
            {},
            {"srb": True, "nrb": True},
        ]:
            assert _assess(**combo).risk_score >= 0


# ---------------------------------------------------------------------------
# NRB suppression flags
# ---------------------------------------------------------------------------


class TestNrbFlags:
    def test_nrb_suppression_active_when_nrb_only(self) -> None:
        a = _assess(nrb=True)
        assert a.nrb_suppression_active is True
        assert a.treatment_underdosed is False

    def test_no_nrb_suppression_when_srb_also_positive(self) -> None:
        a = _assess(srb=True, nrb=True)
        assert a.nrb_suppression_active is False
        assert a.treatment_underdosed is True

    def test_no_nrb_flags_when_nrb_negative(self) -> None:
        a = _assess(srb=True)
        assert a.nrb_suppression_active is False
        assert a.treatment_underdosed is False

    def test_no_nrb_flags_when_all_negative(self) -> None:
        a = _assess()
        assert a.nrb_suppression_active is False
        assert a.treatment_underdosed is False


# ---------------------------------------------------------------------------
# Active guilds and corrosion drivers
# ---------------------------------------------------------------------------


class TestActiveGuilds:
    def test_empty_when_all_negative(self) -> None:
        a = _assess()
        assert a.active_guilds == ()
        assert a.corrosion_drivers == ()

    def test_only_positive_guilds_listed(self) -> None:
        a = _assess(srb=True, irb=True)
        assert "SRB" in a.active_guilds
        assert "IRB" in a.active_guilds
        assert "MCR" not in a.active_guilds
        assert "APB" not in a.active_guilds

    def test_canonical_order_srb_mcr_irb_apb_nrb(self) -> None:
        a = _assess(srb=True, nrb=True, irb=True)
        assert a.active_guilds == ("SRB", "IRB", "NRB")

    def test_all_five_guilds_listed_when_all_positive(self) -> None:
        a = _assess(srb=True, mcr=True, irb=True, apb=True, nrb=True)
        assert len(a.active_guilds) == 5

    def test_corrosion_drivers_count_matches_active_guilds(self) -> None:
        a = _assess(srb=True, irb=True, apb=True)
        assert len(a.corrosion_drivers) == len(a.active_guilds)

    def test_corrosion_driver_contains_gene_name(self) -> None:
        a = _assess(srb=True)
        assert any("dsrB" in d for d in a.corrosion_drivers)

    def test_nrb_corrosion_driver_mentions_suppression(self) -> None:
        a = _assess(nrb=True)
        assert any("suppression" in d.lower() for d in a.corrosion_drivers)


# ---------------------------------------------------------------------------
# Interpretation and recommended action text
# ---------------------------------------------------------------------------


class TestTextContent:
    def test_all_negative_interpretation_mentions_no_guild(self) -> None:
        a = _assess()
        assert "no" in a.interpretation.lower() or "not" in a.interpretation.lower()

    def test_critical_interpretation_mentions_maximum(self) -> None:
        a = _assess(srb=True, mcr=True, irb=True, apb=True)
        text = a.interpretation.lower()
        assert "maximum" in text or "all four" in text

    def test_srb_irb_interpretation_mentions_fes(self) -> None:
        a = _assess(srb=True, irb=True)
        assert "FeS" in a.interpretation or "iron sulfide" in a.interpretation.lower()

    def test_nrb_only_interpretation_mentions_treatment(self) -> None:
        a = _assess(nrb=True)
        text = a.interpretation.lower()
        assert "nitrate" in text or "treatment" in text

    def test_srb_apb_interpretation_mentions_syntrophic(self) -> None:
        a = _assess(srb=True, apb=True)
        assert "syntrophic" in a.interpretation.lower()

    def test_interpretation_is_non_empty_for_all_combos(self) -> None:
        from itertools import product

        keys = ["srb", "mcr", "irb", "apb", "nrb"]
        for bits in product([False, True], repeat=5):
            flags = GuildFlags(**dict(zip(keys, bits, strict=True)))
            a = assess_mic_risk(flags)
            assert a.interpretation.strip()
            assert a.recommended_action.strip()

    def test_recommended_action_is_non_empty(self) -> None:
        for a in [
            _assess(),
            _assess(srb=True),
            _assess(nrb=True),
            _assess(srb=True, irb=True, apb=True),
        ]:
            assert len(a.recommended_action) > 20


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


class TestToDict:
    def test_contains_required_keys(self) -> None:
        d = _assess(srb=True, irb=True).to_dict()
        for key in (
            "risk_level",
            "risk_score",
            "active_guilds",
            "corrosion_drivers",
            "interpretation",
            "recommended_action",
            "nrb_suppression_active",
            "treatment_underdosed",
            "guild_flags",
        ):
            assert key in d, f"Missing key: {key}"

    def test_risk_level_is_string(self) -> None:
        d = _assess().to_dict()
        assert isinstance(d["risk_level"], str)

    def test_guild_flags_inlined(self) -> None:
        d = _assess(srb=True, irb=True).to_dict()
        gf = d["guild_flags"]
        assert isinstance(gf, dict)
        assert gf["srb"] is True
        assert gf["irb"] is True
        assert gf["mcr"] is False

    def test_json_serialisable(self) -> None:
        d = _assess(srb=True, mcr=True, irb=True, apb=True, nrb=True).to_dict()
        serialised = json.dumps(d)
        recovered = json.loads(serialised)
        assert recovered["risk_level"] == "CRITICAL"


# ---------------------------------------------------------------------------
# flags_from_dict
# ---------------------------------------------------------------------------


class TestFlagsFromDict:
    def test_all_true(self) -> None:
        flags = flags_from_dict({"srb": True, "mcr": True, "irb": True, "apb": True, "nrb": True})
        assert flags == GuildFlags(srb=True, mcr=True, irb=True, apb=True, nrb=True)

    def test_missing_keys_default_false(self) -> None:
        flags = flags_from_dict({"srb": True})
        assert flags.srb is True
        assert flags.mcr is False
        assert flags.irb is False

    def test_empty_dict_all_false(self) -> None:
        flags = flags_from_dict({})
        assert flags == GuildFlags(srb=False, mcr=False, irb=False, apb=False, nrb=False)

    def test_unknown_keys_ignored(self) -> None:
        flags = flags_from_dict({"srb": True, "extra_key": "ignored"})
        assert flags.srb is True


# ---------------------------------------------------------------------------
# write_assessment_json
# ---------------------------------------------------------------------------


class TestWriteAssessmentJson:
    def test_creates_file(self, tmp_path: Path) -> None:
        a = _assess(srb=True)
        out = tmp_path / "assessment.json"
        write_assessment_json(a, out)
        assert out.exists()

    def test_valid_json_content(self, tmp_path: Path) -> None:
        a = _assess(srb=True, irb=True)
        out = tmp_path / "assessment.json"
        write_assessment_json(a, out)
        data = json.loads(out.read_text())
        assert data["risk_level"] == "HIGH"
        assert "SRB" in data["active_guilds"]

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "result.json"
        write_assessment_json(_assess(), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# write_assessment_csv
# ---------------------------------------------------------------------------


class TestWriteAssessmentCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(srb=True), out)
        assert out.exists()

    def test_has_header_row(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(), out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == ["field", "value"]

    def test_risk_level_row_present(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(srb=True, irb=True), out)
        rows = list(csv.DictReader(out.open()))
        risk_rows = [r for r in rows if r["field"] == "risk_level"]
        assert len(risk_rows) == 1
        assert risk_rows[0]["value"] == "HIGH"

    def test_guild_flags_present_as_separate_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "assessment.csv"
        write_assessment_csv(_assess(srb=True, irb=True), out)
        rows = list(csv.DictReader(out.open()))
        fields = {r["field"] for r in rows}
        assert "guild_srb" in fields
        assert "guild_irb" in fields
        assert "guild_mcr" in fields


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestMicRiskCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["mic-risk", *args], obj={})

    def test_all_negative_exits_0(self) -> None:
        result = self._run([])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_output_shows_guild_table(self) -> None:
        result = self._run(["--srb", "--irb"])
        output = result.output  # type: ignore[union-attr]
        assert "SRB" in output
        assert "IRB" in output
        assert "[+]" in output

    def test_output_shows_risk_level(self) -> None:
        result = self._run(["--srb", "--irb"])
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_critical_output_for_four_guilds(self) -> None:
        result = self._run(["--srb", "--mcr", "--irb", "--apb"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_minimal_for_no_flags(self) -> None:
        result = self._run([])
        assert "MINIMAL" in result.output  # type: ignore[union-attr]

    def test_nrb_only_shows_low(self) -> None:
        result = self._run(["--nrb"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "LOW" in result.output  # type: ignore[union-attr]

    def test_out_json_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.json"
        result = self._run(["--srb", "--out-json", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["risk_level"] == "MODERATE"

    def test_out_csv_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.csv"
        result = self._run(["--srb", "--irb", "--out-csv", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_flag(self, tmp_path: Path) -> None:
        flags_file = tmp_path / "flags.json"
        flags_file.write_text(
            '{"srb": true, "irb": true, "mcr": false, "apb": false, "nrb": false}'
        )
        result = self._run(["--input-json", str(flags_file)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_nrb_suppression_note_shown(self) -> None:
        result = self._run(["--nrb"])
        output = result.output  # type: ignore[union-attr]
        assert "suppression" in output.lower() or "effective" in output.lower()

    def test_treatment_underdosed_warning_shown(self) -> None:
        result = self._run(["--srb", "--nrb"])
        output = result.output  # type: ignore[union-attr]
        assert (
            "under-dosed" in output
            or "underdosed" in output.lower()
            or "treatment" in output.lower()
        )
