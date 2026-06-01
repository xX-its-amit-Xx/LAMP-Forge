"""Swine neonatal enteric LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable swine-neonatal
enteric LAMP panel as a structured alert, prioritising pathogens by their
case-fatality rate in neonatal piglets and barn-level spread velocity.

This module complements the systemic/respiratory farm-biosecurity panel
(:mod:`lamp_forge.farm_risk`) by covering the three major enteric pathogens
that cause neonatal piglet diarrhoea and mass mortality in farrowing barns.
PEDV alone was documented as Recipe 27 in the cookbook; this module
provides the multi-pathogen interpretation layer.

Pathogens covered
-----------------

    Target     Gene     Pathogen                                CFR neonates
    PEDV       N gene   Porcine epidemic diarrhea virus         ~100%
    PDCOV      N gene   Porcine deltacoronavirus                40-80%
    ROTA_A     VP6      Porcine rotavirus A (group A)           10-40%

Scoring model
-------------
Base weights reflect neonatal case-fatality rate (CFR) and barn-spread
velocity -- the two operational decision variables for a portable assay
deployed in a farrowing barn:

    PEDV   : 50  near-100% CFR in piglets < 1 week; barn-to-barn fomite
                 spread within hours; 2013-14 US outbreak eliminated >10%
                 of the national sow herd
    PDCoV  : 30  40-80% neonatal CFR; lower than PEDV but clinically
                 indistinguishable; emerging pathogen with no approved
                 vaccine in most regions
    ROTA_A : 15  10-40% mortality in untreated neonates; manageable with
                 supportive care (oral rehydration); frequently co-infects
                 with the coronaviruses above, compounding severity

Alert level boundaries
----------------------
    Score >= 50  -> CRITICAL   (PEDV detected; barn lockdown; stop all movement)
    Score >= 30  -> HIGH       (PDCoV detected; intensive management; vet call)
    Score >= 15  -> MODERATE   (Rotavirus only; supportive care; monitor)
    Score >=  1  -> LOW        (any positive; review with veterinarian)
    Score ==  0  -> NEGATIVE   (all negative; confirm extraction control positive)

Immediate quarantine
--------------------
``immediate_quarantine_required`` is set True when PEDV is positive, reflecting
its catastrophic spread dynamics.  A PEDV-positive barn requires immediate
pen quarantine, cessation of all pig movement, and enhanced PPE protocol
before depopulation or treatment decisions are made.

``biosecurity_lockdown`` is set True when any coronavirus (PEDV or PDCoV)
is positive, because both spread via fomites and aerosolised faecal particles
at barn-to-barn distances.

References:
    Stevenson GW et al. (2013) Emergence of PEDV in the United States.
        J Vet Diagn Invest 25(5):649-654. doi:10.1177/1040638713501675
    Chen Q et al. (2015) Porcine deltacoronavirus (PDCoV) in the United States.
        Vet Microbiol 179(3-4):285-293. doi:10.1016/j.vetmic.2015.04.022
    Saif LJ et al. (2012) Rotaviruses. In: Diseases of Swine, 10th ed.
        Wiley-Blackwell, pp 621-634.
    He X et al. (2019) Development and evaluation of a RT-LAMP assay for
        rapid detection of PEDV. Transbound Emerg Dis 66(1):431-439.
        doi:10.1111/tbed.13044
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


class SwineEntericAlertLevel(StrEnum):
    """Swine neonatal enteric LAMP panel alert severity.

    CRITICAL: PEDV detected; near-100% neonatal mortality; barn lockdown
        and immediate cessation of pig movement required.
    HIGH: PDCoV detected without PEDV; 40-80% neonatal mortality; intensive
        supportive care and veterinary management required.
    MODERATE: Rotavirus A only; supportive care with oral rehydration therapy;
        monitor for secondary coronavirus co-infection.
    LOW: Any positive finding; review with site veterinarian.
    NEGATIVE: All targets negative.  Confirm extraction control (16S or
        internal RNA control) is positive before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class SwineEntericFlags:
    """Positivity flags for one swine neonatal enteric LAMP panel result.

    Each flag is True if the corresponding RT-LAMP or LAMP assay returned a
    positive amplification signal for that pathogen.

    Attributes:
        pedv: Porcine epidemic diarrhea virus (N gene) positive.
        pdcov: Porcine deltacoronavirus (N gene) positive.
        rota_a: Porcine rotavirus A (VP6 gene) positive.
    """

    pedv: bool
    pdcov: bool
    rota_a: bool


@dataclass(frozen=True, slots=True)
class SwineEntericAssessment:
    """Structured alert from a swine neonatal enteric LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation.
        immediate_quarantine_required: True when PEDV is positive.  PEDV
            spreads via fomites and aerosolised faeces within hours of the
            first piglet showing clinical signs; pen quarantine must begin
            before a full barn assessment is complete.
        biosecurity_lockdown: True when any coronavirus (PEDV or PDCoV) is
            positive.  Barn-to-barn spread requires enhanced vehicle and
            boot decontamination before any personnel leave the farrowing
            unit.
    """

    flags: SwineEntericFlags
    alert_level: SwineEntericAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    interpretation: str
    recommended_action: str
    immediate_quarantine_required: bool
    biosecurity_lockdown: bool

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; panel flags are inlined.
        """
        return {
            "alert_level": self.alert_level.value,
            "alert_score": self.alert_score,
            "active_targets": list(self.active_targets),
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "immediate_quarantine_required": self.immediate_quarantine_required,
            "biosecurity_lockdown": self.biosecurity_lockdown,
            "panel_flags": {
                "pedv": self.flags.pedv,
                "pdcov": self.flags.pdcov,
                "rota_a": self.flags.rota_a,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_PEDV = 50
_WEIGHT_PDCOV = 30
_WEIGHT_ROTA_A = 15


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> SwineEntericAlertLevel:
    if score >= 50:
        return SwineEntericAlertLevel.CRITICAL
    if score >= 30:
        return SwineEntericAlertLevel.HIGH
    if score >= 15:
        return SwineEntericAlertLevel.MODERATE
    if score >= 1:
        return SwineEntericAlertLevel.LOW
    return SwineEntericAlertLevel.NEGATIVE


def _build_text(
    flags: SwineEntericFlags,
    active_targets: tuple[str, ...],
    alert_level: SwineEntericAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    pedv, pdcov, rota_a = flags.pedv, flags.pdcov, flags.rota_a

    # All negative
    if not any([pedv, pdcov, rota_a]):
        return (
            "All swine neonatal enteric targets negative. No enteric pathogens detected.",
            "Confirm an internal extraction control (e.g. 16S rRNA) is positive "
            "before reporting a clean result. Re-test if watery diarrhoea or "
            "mortality in neonatal piglets is observed.",
        )

    # PEDV positive (dominant signal regardless of co-detection)
    if pedv:
        co = []
        if pdcov:
            co.append("PDCoV")
        if rota_a:
            co.append("Rotavirus A")
        co_str = "; co-detected with " + ", ".join(co) + "." if co else "."
        interp = (
            "Porcine epidemic diarrhea virus (PEDV) detected via the N-gene assay"
            + co_str
            + " PEDV causes near-100% mortality in neonatal piglets under 1 week old "
            "and spreads via fomites and aerosolised faeces within hours of the first "
            "clinical signs. The 2013-14 US outbreak eliminated >10% of the national "
            "sow herd within 12 months of introduction."
        )
        action = (
            "CRITICAL: Immediately quarantine all affected pens. Do NOT move pigs, "
            "personnel, or equipment out of the barn without decontamination. Notify "
            "the site veterinarian and the farm biosecurity coordinator. Implement "
            "barn lockdown (dedicated PPE, boot dips, one-way traffic flow). Collect "
            "intestinal content samples from acutely affected piglets for confirmatory "
            "RT-PCR at an accredited laboratory. Contact the farm's swine health "
            "programme and review gilt/sow vaccination history."
        )
        return interp, action

    # PDCoV positive (no PEDV)
    if pdcov:
        co_str = " Co-detected with Rotavirus A." if rota_a else ""
        interp = (
            "Porcine deltacoronavirus (PDCoV) detected via the N-gene assay."
            + co_str
            + " PDCoV causes 40-80% mortality in neonatal piglets and clinically "
            "indistinguishable watery diarrhoea from PEDV. It is an emerging pathogen "
            "first identified in the United States in 2014; no approved vaccine is "
            "available in most jurisdictions."
        )
        action = (
            "HIGH: Quarantine the affected farrowing unit. Notify the site veterinarian. "
            "Initiate enhanced biosecurity (dedicated PPE, one-way traffic flow, "
            "vehicle decontamination). Collect samples for confirmatory RT-PCR; "
            "PDCoV can co-circulate with PEDV and a negative PEDV result does not "
            "exclude mixed infection. Provide intensive supportive care: oral "
            "rehydration therapy and warmth for affected litters. Monitor adjacent "
            "pens for clinical signs."
        )
        return interp, action

    # Rotavirus A only
    if rota_a and not pedv and not pdcov:
        interp = (
            "Porcine rotavirus A detected via the VP6 assay. Rotavirus A is the most "
            "common cause of neonatal diarrhoea in swine, with 10-40% mortality in "
            "untreated neonates. It is manageable with supportive care but frequently "
            "co-circulates with PEDV and PDCoV; a negative coronavirus result does not "
            "exclude co-infection if clinical signs are severe."
        )
        action = (
            "Notify the site veterinarian. Initiate oral rehydration therapy (electrolyte "
            "solution) for affected litters. Ensure adequate creep-feeding heat and "
            "colostrum access. Monitor the pen and adjacent pens for 24-48 h for "
            "escalation to watery diarrhoea consistent with enteric coronavirus. "
            "Re-test with the full panel if mortality exceeds 20% in the affected litter."
        )
        return interp, action

    # Fallback
    targets_str = ", ".join(active_targets)
    return (
        f"Detected: {targets_str}. {alert_level.value} alert.",
        "Review with site veterinarian. Follow standard enteric disease protocol.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_swine_enteric_risk(flags: SwineEntericFlags) -> SwineEntericAssessment:
    """Compute a structured enteric-pathogen alert from three-target panel flags.

    Weights each pathogen by neonatal case-fatality rate and spread velocity,
    sums to a 0-100 score, and maps it to one of five categorical alert levels.
    Sets ``immediate_quarantine_required`` when PEDV is positive and
    ``biosecurity_lockdown`` when any coronavirus (PEDV or PDCoV) is positive.

    Args:
        flags: Positivity flags for all three swine neonatal enteric targets.

    Returns:
        :class:`SwineEntericAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.swine_enteric_risk import SwineEntericFlags, assess_swine_enteric_risk

        flags = SwineEntericFlags(pedv=True, pdcov=False, rota_a=False)
        assessment = assess_swine_enteric_risk(flags)
        print(assessment.alert_level)                  # SwineEntericAlertLevel.CRITICAL
        print(assessment.immediate_quarantine_required)  # True
    """
    score = 0
    if flags.pedv:
        score += _WEIGHT_PEDV
    if flags.pdcov:
        score += _WEIGHT_PDCOV
    if flags.rota_a:
        score += _WEIGHT_ROTA_A

    score = max(0, min(100, score))

    target_order = [
        ("PEDV", flags.pedv),
        ("PDCoV", flags.pdcov),
        ("RotaA", flags.rota_a),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, active_targets, alert_level)

    return SwineEntericAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        interpretation=interpretation,
        recommended_action=recommended_action,
        immediate_quarantine_required=flags.pedv,
        biosecurity_lockdown=flags.pedv or flags.pdcov,
    )


def flags_from_dict(data: dict[str, object]) -> SwineEntericFlags:
    """Construct :class:`SwineEntericFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``pedv``, ``pdcov``, ``rota_a`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`SwineEntericFlags` instance.
    """
    return SwineEntericFlags(
        pedv=bool(data.get("pedv", False)),
        pdcov=bool(data.get("pdcov", False)),
        rota_a=bool(data.get("rota_a", False)),
    )


def write_assessment_json(assessment: SwineEntericAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`SwineEntericAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: SwineEntericAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`SwineEntericAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        d = assessment.to_dict()
        for key, val in d.items():
            if key == "active_targets":
                writer.writerow([key, "; ".join(val)])  # type: ignore[arg-type]
            elif key == "panel_flags" and isinstance(val, dict):
                for pk, pv in val.items():
                    writer.writerow([f"target_{pk}", "positive" if pv else "negative"])
            else:
                writer.writerow([key, val])
