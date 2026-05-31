"""Farm-biosecurity LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable farm-biosecurity
LAMP panel as a structured alert, prioritising pathogens by their regulatory
status, case-fatality rate, and field-containment urgency.

Pathogens covered
-----------------

    Target    Gene     Pathogen                           Regulatory status
    ASFV      B646L    African swine fever virus          WOAH-listed Tier 1
    FMDV      3Dpol    Foot-and-mouth disease virus       WOAH-listed Tier 1
    AIV       M gene   Avian influenza A (all subtypes)   WOAH-listed (HPAI)
    NDV       M gene   Newcastle disease virus            WOAH-listed
    PRRSV     ORF7     Porcine reproductive & resp. syndrome virus  Not listed

The assay design (Recipes 7, 12, 14, 18, 11 in the cookbook) targets the most
conserved gene per pathogen to maximise broad-strain coverage.  The AIV M-gene
assay is pan-influenza-A; a positive should be treated as a presumptive HPAI
detection pending subtype confirmation by a reference laboratory.  The NDV
M-gene assay is pan-NDV; a positive requires pathotyping of the F-gene
cleavage site to distinguish velogenic from lentogenic strains.

Scoring model
-------------
Base weights reflect case-fatality rate (CFR), economic impact, and
contagiousness in their primary host species:

    ASFV  : 60  near-100% CFR in pigs; no widely deployed vaccine; WOAH Tier 1
    FMDV  : 55  WOAH Tier 1; 2001 UK outbreak ~GBP 8 billion; 7 antigenically
                distinct serotypes; affects cattle, pigs, sheep, goats
    AIV   : 50  HPAI H5N1 clade 2.3.4.4b catastrophic poultry losses; zoonotic
                spillover risk; WOAH-listed
    NDV   : 25  WOAH-listed; velogenic strains approach 100% flock mortality;
                lentogenic strains subclinical -- see score note below
    PRRSV : 20  not WOAH-listed; USD 664 M/yr US losses; endemic in most regions

Alert level boundaries
----------------------
    Score >= 50  -> CRITICAL   (immediate notifiable-disease response)
    Score >= 25  -> HIGH       (quarantine, veterinary attention, reporting)
    Score >= 20  -> MODERATE   (veterinary investigation, movement restriction)
    Score >=  1  -> LOW        (positive detection in a managed or non-reportable target)
    Score ==  0  -> NEGATIVE   (all negative; sample adequacy check recommended)

Immediate reporting
-------------------
``immediate_report_required`` is set True whenever any WOAH-listed pathogen
(ASFV, FMDV, AIV, NDV) is positive.  In the field, a positive result on any
of these four targets should trigger notification to the national veterinary
authority and confirmatory testing at a WOAH Reference Laboratory before
depopulation or trade restrictions are enforced.  PRRSV alone does not require
immediate national reporting in most jurisdictions.

References:
    WOAH (2024) List of OIE-notifiable diseases, infections and infestations.
        https://www.woah.org/en/what-we-do/animal-health-and-welfare/
        disease-data-collection/list-of-diseases/
    Arias M et al. (2021) African swine fever. Transbound Emerg Dis 68:1157.
    Grubman MJ & Baxt B (2004) Foot-and-mouth disease. Clin Microbiol Rev
        17(2):465-493.
    Swayne DE (2012) Impact of vaccines and vaccination on global control of
        avian influenza. Avian Dis 56(4 suppl):818-828.
    Holtkamp DJ et al. (2013) Assessment of PRRSV economic impact. J Swine
        Health Prod 21:72-84.
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


class FarmAlertLevel(StrEnum):
    """Farm-biosecurity LAMP panel alert severity.

    CRITICAL: WOAH-listed notifiable pathogen detected; immediate response
        and confirmatory testing required.  Do not move animals.
    HIGH: WOAH-listed pathogen detected or multi-target outbreak scenario;
        veterinary investigation and quarantine within 24 h.
    MODERATE: Actionable finding in a managed or non-notifiable pathogen;
        veterinary attention and movement restriction.
    LOW: Positive in a non-notifiable, usually endemic target; review
        with site veterinarian.
    NEGATIVE: All targets negative.  If a 16S sample-adequacy control is
        running, confirm it is positive before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class FarmPanelFlags:
    """Positivity flags for one farm-biosecurity LAMP panel result.

    Each flag is True if the corresponding LAMP assay returned a positive
    amplification signal for that pathogen.

    Attributes:
        asfv: African swine fever virus (B646L / p72 capsid) positive.
        fmdv: Foot-and-mouth disease virus (3Dpol RNA polymerase) positive.
        aiv: Avian influenza A virus (M gene, pan-IAV) positive.
        ndv: Newcastle disease virus (M gene, pan-NDV) positive.
        prrsv: Porcine reproductive and respiratory syndrome virus (ORF7) positive.
    """

    asfv: bool
    fmdv: bool
    aiv: bool
    ndv: bool
    prrsv: bool


@dataclass(frozen=True, slots=True)
class FarmRiskAssessment:
    """Structured alert from a farm-biosecurity LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        notifiable_targets: Subset of active_targets that are WOAH-listed.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation.
        immediate_report_required: True when any WOAH-listed pathogen is
            positive (ASFV, FMDV, AIV, NDV).
    """

    flags: FarmPanelFlags
    alert_level: FarmAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    notifiable_targets: tuple[str, ...]
    interpretation: str
    recommended_action: str
    immediate_report_required: bool

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; panel flags are inlined.
        """
        return {
            "alert_level": self.alert_level.value,
            "alert_score": self.alert_score,
            "active_targets": list(self.active_targets),
            "notifiable_targets": list(self.notifiable_targets),
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "immediate_report_required": self.immediate_report_required,
            "panel_flags": {
                "asfv": self.flags.asfv,
                "fmdv": self.flags.fmdv,
                "aiv": self.flags.aiv,
                "ndv": self.flags.ndv,
                "prrsv": self.flags.prrsv,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_ASFV = 60
_WEIGHT_FMDV = 55
_WEIGHT_AIV = 50
_WEIGHT_NDV = 25
_WEIGHT_PRRSV = 20

# WOAH-listed pathogens in this panel (excluding PRRSV)
_NOTIFIABLE = frozenset({"ASFV", "FMDV", "AIV", "NDV"})

_TARGET_DESCRIPTION: dict[str, str] = {
    "ASFV": "African swine fever virus (B646L / p72 capsid) -- WOAH Tier 1 notifiable",
    "FMDV": "Foot-and-mouth disease virus (3Dpol RNA polymerase) -- WOAH Tier 1 notifiable",
    "AIV": "Avian influenza A virus (M gene, pan-IAV) -- presumptive HPAI; confirm subtype",
    "NDV": "Newcastle disease virus (M gene, pan-NDV) -- WOAH notifiable; pathotype TBC",
    "PRRSV": "Porcine reproductive and respiratory syndrome virus (ORF7)",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> FarmAlertLevel:
    if score >= 50:
        return FarmAlertLevel.CRITICAL
    if score >= 25:
        return FarmAlertLevel.HIGH
    if score >= 20:
        return FarmAlertLevel.MODERATE
    if score >= 1:
        return FarmAlertLevel.LOW
    return FarmAlertLevel.NEGATIVE


def _build_text(
    flags: FarmPanelFlags,
    active_targets: tuple[str, ...],
    notifiable_targets: tuple[str, ...],
    alert_level: FarmAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        notifiable_targets: WOAH-listed subset of active_targets.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    asfv, fmdv, aiv, ndv, prrsv = (
        flags.asfv,
        flags.fmdv,
        flags.aiv,
        flags.ndv,
        flags.prrsv,
    )

    # All negative
    if not any([asfv, fmdv, aiv, ndv, prrsv]):
        return (
            "All farm-biosecurity targets negative. No notifiable pathogens detected.",
            "Confirm a 16S sample-adequacy control is positive before reporting a "
            "clean result. Re-test if clinical signs are observed.",
        )

    # ASFV positive
    if asfv:
        co = [t for t in ("FMDV", "AIV", "NDV", "PRRSV") if getattr(flags, t.lower())]
        co_str = ", ".join(co) if co else ""
        interp = (
            "ASFV (African swine fever virus) detected via the B646L / p72 assay. "
            "ASFV is a WOAH Tier 1 notifiable disease with near-100% case fatality "
            "in virulent genotype II strains. No widely deployed vaccine is available."
        )
        if co_str:
            interp += f" Co-detections: {co_str}."
        action = (
            "IMMEDIATE: Quarantine all pigs on the premises. Do not move any animals. "
            "Notify the national veterinary authority and request confirmatory testing "
            "at a WOAH Reference Laboratory. Contact the site veterinarian and follow "
            "national ASFV response protocol (stamping-out, biosecurity zone establishment)."
        )
        return interp, action

    # FMDV positive
    if fmdv:
        co = [t for t in ("AIV", "NDV", "PRRSV") if getattr(flags, t.lower())]
        co_str = ", ".join(co) if co else ""
        interp = (
            "FMDV (foot-and-mouth disease virus) detected via the 3Dpol assay "
            "(pan-serotype: O, A, C, SAT1, SAT2, SAT3, Asia1). FMDV is a WOAH Tier 1 "
            "notifiable disease affecting cloven-hoofed livestock. The 2001 UK outbreak "
            "cost an estimated GBP 8 billion."
        )
        if co_str:
            interp += f" Co-detections: {co_str}."
        action = (
            "IMMEDIATE: Quarantine all susceptible animals (cattle, pigs, sheep, goats). "
            "Do not move animals, vehicles, or personnel off the premises. "
            "Notify the national veterinary authority. Request serotype confirmation "
            "at a WOAH Reference Laboratory (e.g. Pirbright Institute, FAO World Reference "
            "Laboratory for FMD). This is a triage result -- do not cull without "
            "confirmatory testing."
        )
        return interp, action

    # AIV positive (no ASFV or FMDV)
    if aiv:
        co = [t for t in ("NDV", "PRRSV") if getattr(flags, t.lower())]
        co_str = ", ".join(co) if co else ""
        interp = (
            "Avian influenza A (pan-IAV M-gene) detected. This assay does not "
            "discriminate HPAI from LPAI -- a positive result is a presumptive HPAI "
            "screening result. HPAI H5N1 clade 2.3.4.4b has caused catastrophic "
            "poultry losses globally since 2021 and carries zoonotic spillover risk."
        )
        if co_str:
            interp += f" Co-detections: {co_str}."
        action = (
            "IMMEDIATE: Quarantine the affected poultry house. Do not move birds, "
            "eggs, or fomites off the premises. Notify the national veterinary "
            "authority and request subtype confirmation at a national reference "
            "laboratory. Biosafety: treat positive samples as potentially infectious "
            "for humans (PPE level appropriate for HPAI response) until LPAI is "
            "confirmed. Do not cull without regulatory guidance."
        )
        return interp, action

    # NDV positive (no ASFV, FMDV, or AIV)
    if ndv:
        interp = (
            "Newcastle disease virus (pan-NDV M-gene) detected. This assay does not "
            "discriminate velogenic (near-100% flock mortality) from lentogenic "
            "(subclinical) strains. Velogenic NDV (vNDV) is WOAH-listed and a "
            "USDA-regulated select agent in the United States."
        )
        if prrsv:
            interp += (
                " PRRSV also detected: mixed avian-porcine outbreak scenario or dual-species site."
            )
        action = (
            "Quarantine the affected flock. Notify the national veterinary authority; "
            "vNDV is a reportable disease in most jurisdictions. Request F-gene "
            "pathotyping to distinguish velogenic from lentogenic strains before "
            "depopulation decisions. Wet-lab confirmation at a WOAH Reference "
            "Laboratory (e.g. USDA APHIS NVSL) is required before regulatory action."
        )
        return interp, action

    # PRRSV alone (only non-notifiable target)
    if prrsv and not any([asfv, fmdv, aiv, ndv]):
        interp = (
            "PRRSV (porcine reproductive and respiratory syndrome virus) detected via "
            "the ORF7 assay. PRRSV is endemic in most swine-producing regions but "
            "causes significant reproductive failure and respiratory disease. It is "
            "not WOAH-listed and does not require immediate national reporting."
        )
        action = (
            "Notify the site veterinarian. Implement standard PRRSV biosecurity "
            "(restrict pig movement, consider PRRSV-status-based segregation). "
            "Review vaccination programme if applicable. Monitor for co-circulation "
            "of notifiable pathogens (ASFV, FMDV) in the herd."
        )
        return interp, action

    # Fallback
    targets_str = ", ".join(active_targets)
    interp = f"Detected targets: {targets_str}. {alert_level.value} alert level."
    action = (
        "Review detection pattern with site veterinarian. Follow national "
        "protocols for each detected pathogen."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_farm_risk(flags: FarmPanelFlags) -> FarmRiskAssessment:
    """Compute a structured farm-biosecurity alert from five-pathogen panel flags.

    Weights each pathogen by case-fatality rate, economic impact, and regulatory
    status, sums to a 0-100 score, and maps it to one of five categorical alert
    levels.  Sets ``immediate_report_required`` when any WOAH-listed pathogen
    is positive.

    Args:
        flags: Positivity flags for all five farm-biosecurity panel targets.

    Returns:
        :class:`FarmRiskAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.farm_risk import FarmPanelFlags, assess_farm_risk

        flags = FarmPanelFlags(asfv=True, fmdv=False, aiv=False, ndv=False, prrsv=False)
        assessment = assess_farm_risk(flags)
        print(assessment.alert_level)   # FarmAlertLevel.CRITICAL
        print(assessment.immediate_report_required)   # True
    """
    score = 0
    if flags.asfv:
        score += _WEIGHT_ASFV
    if flags.fmdv:
        score += _WEIGHT_FMDV
    if flags.aiv:
        score += _WEIGHT_AIV
    if flags.ndv:
        score += _WEIGHT_NDV
    if flags.prrsv:
        score += _WEIGHT_PRRSV

    score = max(0, min(100, score))

    target_order = [
        ("ASFV", flags.asfv),
        ("FMDV", flags.fmdv),
        ("AIV", flags.aiv),
        ("NDV", flags.ndv),
        ("PRRSV", flags.prrsv),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)
    notifiable_targets = tuple(t for t in active_targets if t in _NOTIFIABLE)

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(
        flags, active_targets, notifiable_targets, alert_level
    )

    return FarmRiskAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        notifiable_targets=notifiable_targets,
        interpretation=interpretation,
        recommended_action=recommended_action,
        immediate_report_required=bool(notifiable_targets),
    )


def flags_from_dict(data: dict[str, object]) -> FarmPanelFlags:
    """Construct :class:`FarmPanelFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``asfv``, ``fmdv``, ``aiv``, ``ndv``, ``prrsv`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`FarmPanelFlags` instance.
    """
    return FarmPanelFlags(
        asfv=bool(data.get("asfv", False)),
        fmdv=bool(data.get("fmdv", False)),
        aiv=bool(data.get("aiv", False)),
        ndv=bool(data.get("ndv", False)),
        prrsv=bool(data.get("prrsv", False)),
    )


def write_assessment_json(assessment: FarmRiskAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`FarmRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: FarmRiskAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`FarmRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        d = assessment.to_dict()
        for key, val in d.items():
            if key in ("active_targets", "notifiable_targets"):
                writer.writerow([key, "; ".join(val)])  # type: ignore[arg-type]
            elif key == "panel_flags" and isinstance(val, dict):
                for pk, pv in val.items():
                    writer.writerow([f"target_{pk}", "positive" if pv else "negative"])
            else:
                writer.writerow([key, val])
