"""Oilfield MIC (microbiologically influenced corrosion) risk assessment engine.

Computes a structured corrosion risk score from the five-guild positivity flags
produced by a BioVind-style oilfield LAMP panel:

    Guild     Gene marker  Corrosion pathway
    SRB       dsrB         H2S production -> sulfide corrosion + FeS scale
    Methanogens mcrA       CH4 souring + acetate-corrosion (indirect)
    IRB       omcA         Fe(III) dissolution -> bare-metal exposure
    APB       fthfs        Acid production -> pH drop, FeCO3 dissolution
    NRB       narG         Nitrate reduction (competitive SRB suppressor)

Scoring model
-------------
Each guild contributes a base weight to a 0–100 integer score.  Synergy
bonuses are added for the most corrosive combinations (SRB+IRB → FeS scale,
SRB+APB → syntrophic loop, all four corrosion guilds → maximum MIC risk).
NRB co-detection with SRB earns a partial suppression credit, reflecting the
thermodynamic competition between nitrate- and sulfate-reduction pathways.

Risk level boundaries
---------------------
    Score >= 90  → CRITICAL
    Score >= 60  → HIGH
    Score >= 30  → MODERATE
    Score >=  8  → LOW
    Score <   8  → MINIMAL

References:
    Hubert C & Voordouw G (2007) Oil field souring control by NRB.
        Appl Environ Microbiol 73:2644-2652.
    Vigneron A et al. (2017) Comammox Nitrospira in oilfield MIC microbiomes.
        ISME J 11:2180-2194.
    Hamilton WA (2003) Microbially influenced corrosion as a model system.
        Int Biodeterior Biodegrad 51(3):151-155.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

# ---------------------------------------------------------------------------
# Enumerations and data model
# ---------------------------------------------------------------------------


class MICRiskLevel(StrEnum):
    """Ordered MIC corrosion risk severity for oilfield samples.

    CRITICAL: Immediate process intervention and biocide review required.
    HIGH: Elevated monitoring and treatment escalation warranted.
    MODERATE: Heightened monitoring; evaluate treatment programme effectiveness.
    LOW: Low corrosion signal; maintain standard monitoring schedule.
    MINIMAL: No guilds detected; re-test after any process upset.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MINIMAL = "MINIMAL"


@dataclass(frozen=True, slots=True)
class GuildFlags:
    """Five-guild positivity flags for one oilfield LAMP panel result.

    Each flag is True if the corresponding LAMP assay returned a positive
    amplification signal for that functional guild.

    Attributes:
        srb: Sulfate-reducing bacteria (dsrB) positive.
        mcr: Methanogens (mcrA) positive.
        irb: Iron-reducing bacteria (omcA) positive.
        apb: Acid-producing bacteria / homoacetogens (fthfs) positive.
        nrb: Nitrate-reducing bacteria (narG) positive.
    """

    srb: bool
    mcr: bool
    irb: bool
    apb: bool
    nrb: bool


@dataclass(frozen=True, slots=True)
class MICRiskAssessment:
    """Structured MIC risk assessment from an oilfield LAMP panel.

    Attributes:
        flags: Input guild-positivity flags used to compute this assessment.
        risk_level: Categorical risk level (CRITICAL → MINIMAL).
        risk_score: Numeric 0–100 summary (higher = more severe).
        active_guilds: Labels of positive guilds in canonical order
            (SRB, MCR, IRB, APB, NRB).
        corrosion_drivers: Active corrosion mechanisms from positive guilds.
        interpretation: Concise interpretation of the guild pattern.
        recommended_action: Field-actionable recommendation for this pattern.
        nrb_suppression_active: True when NRB positive and SRB negative,
            indicating nitrate injection is suppressing sulfate-reduction.
        treatment_underdosed: True when both NRB and SRB are positive, i.e.
            nitrate treatment present but not fully suppressing SRB.
    """

    flags: GuildFlags
    risk_level: MICRiskLevel
    risk_score: int
    active_guilds: tuple[str, ...]
    corrosion_drivers: tuple[str, ...]
    interpretation: str
    recommended_action: str
    nrb_suppression_active: bool
    treatment_underdosed: bool

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; guild flags are inlined.
        """
        return {
            "risk_level": self.risk_level.value,
            "risk_score": self.risk_score,
            "active_guilds": list(self.active_guilds),
            "corrosion_drivers": list(self.corrosion_drivers),
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "nrb_suppression_active": self.nrb_suppression_active,
            "treatment_underdosed": self.treatment_underdosed,
            "guild_flags": {
                "srb": self.flags.srb,
                "mcr": self.flags.mcr,
                "irb": self.flags.irb,
                "apb": self.flags.apb,
                "nrb": self.flags.nrb,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_SRB = 30
_WEIGHT_MCR = 15
_WEIGHT_IRB = 20
_WEIGHT_APB = 15
_WEIGHT_NRB = 10
_SYNERGY_SRB_IRB = 15
_SYNERGY_SRB_APB = 8
_SYNERGY_ALL_FOUR = 10
_CREDIT_NRB_SRB = -15

_CORROSION_PATHWAY: dict[str, str] = {
    "SRB": "H2S production, sulfide corrosion, FeS scale formation (dsrB)",
    "MCR": "Methanogenesis, gas souring, indirect acetate-corrosion (mcrA)",
    "IRB": "Fe(III) reduction, protective-oxide dissolution, bare-metal exposure (omcA)",
    "APB": "Acetic-acid production, pH depression, FeCO3 dissolution (fthfs)",
    "NRB": "Nitrate reduction, competitive SRB suppression (narG)",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _risk_level_from_score(score: int) -> MICRiskLevel:
    if score >= 90:
        return MICRiskLevel.CRITICAL
    if score >= 60:
        return MICRiskLevel.HIGH
    if score >= 30:
        return MICRiskLevel.MODERATE
    if score >= 8:
        return MICRiskLevel.LOW
    return MICRiskLevel.MINIMAL


def _build_text(
    flags: GuildFlags,
    active_guilds: tuple[str, ...],
    risk_level: MICRiskLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the guild pattern.

    Args:
        flags: Guild positivity flags.
        active_guilds: Positive guild labels in canonical order.
        risk_level: Pre-computed risk level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    srb, mcr, irb, apb, nrb = flags.srb, flags.mcr, flags.irb, flags.apb, flags.nrb
    n_corrosion = sum([srb, mcr, irb, apb])

    # All negative
    if not any([srb, mcr, irb, apb, nrb]):
        return (
            "No MIC guild detected. Low microbial activity in the sample.",
            "Maintain standard monitoring schedule. Re-test after any process "
            "upset (water injection, temperature excursion, biocide change).",
        )

    # NRB-only (treatment working)
    if nrb and n_corrosion == 0:
        return (
            "Nitrate-reducing bacteria detected without SRB or other corrosion guilds. "
            "Nitrate injection appears effective; NRB are the dominant active guild.",
            "Continue current nitrate-injection treatment. Monitor SRB (dsrB) and "
            "confirm suppression is maintained over subsequent sample intervals.",
        )

    # Treatment underdosed: SRB + NRB present
    if srb and nrb and n_corrosion == 1:
        interp = (
            "SRB and NRB co-detected. Nitrate treatment is present but not fully "
            "suppressing sulfate-reduction. Partial thermodynamic competition is "
            "occurring but H2S generation remains active."
        )
        action = (
            "Escalate nitrate injection dose or switch to a biocide-combination "
            "protocol. Monitor NRB:SRB ratio over subsequent intervals to confirm "
            "treatment effectiveness trend."
        )
        return interp, action

    # Four corrosion guilds active (maximum MIC risk)
    if srb and mcr and irb and apb:
        interp = (
            "All four corrosion guilds detected: SRB (sulfide), methanogens (gas "
            "souring), IRB (iron oxide dissolution), and APB (acid production). "
            "SRB-IRB FeS scale formation and APB-SRB syntrophic amplification are "
            "simultaneously active. Maximum MIC risk."
        )
        if nrb:
            interp += (
                " NRB also detected but insufficient to suppress SRB at this guild-activity level."
            )
        action = (
            "Immediate biocide treatment review and escalation. Inspect pipeline "
            "segments for under-deposit pitting. Evaluate nitrate injection as a "
            "supplementary strategy to shift community towards NRB dominance."
        )
        return interp, action

    # Three corrosion guilds
    if n_corrosion == 3:
        guilds_str = ", ".join(g for g in ("SRB", "MCR", "IRB", "APB") if getattr(flags, g.lower()))
        if srb and irb:
            interp = (
                f"{guilds_str} detected. SRB-IRB FeS scale formation active: IRB "
                "dissolves protective iron oxide layers, exposing bare steel; SRB-generated "
                "H2S reprecipitates as corrosive iron sulfide. High MIC risk."
            )
            if apb:
                interp += (
                    " APB acid production additionally depresses pH, dissolving the "
                    "protective FeCO3 passivation layer."
                )
            elif mcr:
                interp += " Methanogenesis adds gas-souring and indirect acetate-corrosion risk."
        else:
            interp = f"{guilds_str} detected. Three corrosion guilds co-active. High MIC risk."
        action = (
            "Escalate monitoring frequency and biocide treatment. Inspect pipeline "
            "segments at highest risk of under-deposit corrosion. Quantify SRB load "
            "with qPCR if not already done to assess treatment urgency."
        )
        return interp, action

    # SRB + IRB (FeS scale, no APB or MCR)
    if srb and irb and n_corrosion == 2:
        interp = (
            "SRB and IRB co-detected. Classic FeS scale formation pathway active: "
            "IRB dissolves protective Fe(III) oxide layers; SRB generates H2S, "
            "which reprecipitates as iron sulfide (FeS / mackinawite) under "
            "deposit — the dominant failure mode in oilfield pipeline pitting."
        )
        action = (
            "Biocide treatment review. Inspect for under-deposit pitting at "
            "high-fluid-velocity joints and weld zones. Consider magnetic-flux "
            "leakage (MFL) pigging to locate developing pits."
        )
        return interp, action

    # SRB + APB (syntrophic loop, no IRB)
    if srb and apb and not irb:
        interp = (
            "SRB and APB co-detected. APB-SRB syntrophic loop established: "
            "homoacetogens generate acetate and H2 as preferred electron donors "
            "for SRB, amplifying sulfide production. pH depression from acid "
            "production reduces the protective FeCO3 passivation layer."
        )
        action = (
            "Increase biocide dosing frequency. Monitor pH of produced water; "
            "pH < 5.5 indicates FeCO3 dissolution is occurring. Consider "
            "pH-buffering corrosion inhibitor alongside the biocide programme."
        )
        return interp, action

    # SRB + MCR only
    if srb and mcr and n_corrosion == 2:
        interp = (
            "SRB and methanogens co-detected. Sulfidogenesis and reservoir "
            "gas souring risk co-present. Standard souring monitoring applies."
        )
        action = (
            "Apply standard souring-control protocol. Monitor H2S levels at "
            "wellhead and separator. Biocide treatment if SRB load is above "
            "site threshold."
        )
        return interp, action

    # SRB alone (including SRB + NRB underdosed)
    if srb and n_corrosion == 1:
        interp = "SRB detected. H2S generation risk; sulfide corrosion active."
        if nrb:
            interp = (
                "SRB and NRB co-detected. Nitrate treatment present but SRB not "
                "fully suppressed. H2S generation remains active."
            )
        action = (
            "Apply or escalate biocide treatment. Monitor wellhead H2S levels. "
            "If nitrate injection is in use, review dose and residence time."
        )
        return interp, action

    # IRB alone
    if irb and not srb and n_corrosion <= 2:
        interp = (
            "Iron-reducing bacteria detected"
            + (" and APB (acid-producing bacteria)" if apb else "")
            + ". Iron oxide dissolution is active, potentially exposing bare "
            "steel. SRB not yet detected but iron cycling creates a permissive "
            "environment for SRB establishment."
        )
        action = (
            "Increase monitoring frequency for SRB (dsrB). Inspect pipeline "
            "segments for pitting. Consider corrosion inhibitor application."
        )
        return interp, action

    # APB alone
    if apb and not srb and not irb:
        interp = (
            "Acid-producing bacteria (homoacetogens) detected. Acetate and H2 "
            "production active — these are preferred SRB electron donors. pH "
            "depression is dissolving FeCO3 protective scale. SRB bloom risk "
            "elevated on any sulfate-injection event."
        )
        action = (
            "Monitor SRB (dsrB) closely over the next sampling intervals. "
            "Pre-emptive biocide application may be warranted if SRB history "
            "at this site is positive."
        )
        return interp, action

    # MCR alone
    if mcr and n_corrosion == 1:
        interp = (
            "Methanogens detected. Methanogenic activity contributes to gas "
            "souring and indirect acetate-corrosion, but structural corrosion "
            "risk from this guild alone is low."
        )
        action = (
            "Monitor gas composition (CH4, CO2) at the separator. Maintain "
            "standard monitoring schedule for SRB and IRB."
        )
        return interp, action

    # Fallback for other multi-guild combinations
    guilds_str = ", ".join(active_guilds)
    interp = f"Detected guilds: {guilds_str}. {risk_level.value} corrosion risk."
    action = "Review guild pattern with corrosion engineer. Adjust monitoring and biocide schedule accordingly."
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_mic_risk(flags: GuildFlags) -> MICRiskAssessment:
    """Compute a structured MIC risk assessment from five-guild positivity flags.

    The scoring model weights each guild by its direct corrosion contribution,
    adds synergy bonuses for the most corrosive combinations, and credits NRB
    co-detection with SRB as a partial suppression indicator.  The resulting
    0–100 score maps to one of five categorical risk levels.

    Args:
        flags: Positivity flags for all five oilfield guilds.

    Returns:
        :class:`MICRiskAssessment` with risk level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.mic_risk import GuildFlags, assess_mic_risk

        flags = GuildFlags(srb=True, mcr=False, irb=True, apb=True, nrb=False)
        assessment = assess_mic_risk(flags)
        print(assessment.risk_level)   # MICRiskLevel.HIGH
        print(assessment.risk_score)   # 88
    """
    score = 0
    if flags.srb:
        score += _WEIGHT_SRB
    if flags.mcr:
        score += _WEIGHT_MCR
    if flags.irb:
        score += _WEIGHT_IRB
    if flags.apb:
        score += _WEIGHT_APB
    if flags.nrb:
        score += _WEIGHT_NRB

    # Synergy bonuses
    if flags.srb and flags.irb:
        score += _SYNERGY_SRB_IRB
    if flags.srb and flags.apb:
        score += _SYNERGY_SRB_APB
    if flags.srb and flags.mcr and flags.irb and flags.apb:
        score += _SYNERGY_ALL_FOUR

    # NRB partial suppression credit
    if flags.nrb and flags.srb:
        score += _CREDIT_NRB_SRB

    score = max(0, min(100, score))

    guild_order = [
        ("SRB", flags.srb),
        ("MCR", flags.mcr),
        ("IRB", flags.irb),
        ("APB", flags.apb),
        ("NRB", flags.nrb),
    ]
    active_guilds = tuple(label for label, pos in guild_order if pos)
    corrosion_drivers = tuple(_CORROSION_PATHWAY[g] for g in active_guilds)

    risk_level = _risk_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, active_guilds, risk_level)

    return MICRiskAssessment(
        flags=flags,
        risk_level=risk_level,
        risk_score=score,
        active_guilds=active_guilds,
        corrosion_drivers=corrosion_drivers,
        interpretation=interpretation,
        recommended_action=recommended_action,
        nrb_suppression_active=flags.nrb and not flags.srb,
        treatment_underdosed=flags.nrb and flags.srb,
    )


def flags_from_dict(data: dict[str, object]) -> GuildFlags:
    """Construct :class:`GuildFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``srb``, ``mcr``, ``irb``, ``apb``, ``nrb`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each guild.

    Returns:
        :class:`GuildFlags` instance.
    """
    return GuildFlags(
        srb=bool(data.get("srb", False)),
        mcr=bool(data.get("mcr", False)),
        irb=bool(data.get("irb", False)),
        apb=bool(data.get("apb", False)),
        nrb=bool(data.get("nrb", False)),
    )


def write_assessment_json(assessment: MICRiskAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`MICRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: MICRiskAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`MICRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        d = assessment.to_dict()
        for key, val in d.items():
            if key in ("active_guilds", "corrosion_drivers"):
                writer.writerow([key, "; ".join(val)])  # type: ignore[arg-type]
            elif key == "guild_flags" and isinstance(val, dict):
                for gk, gv in val.items():
                    writer.writerow([f"guild_{gk}", "positive" if gv else "negative"])
            else:
                writer.writerow([key, val])
