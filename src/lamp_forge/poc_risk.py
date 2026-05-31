"""Human point-of-care (POC) LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable human POC LAMP
panel as a structured clinical alert, prioritising pathogens by their public
health urgency, treatment implications, and isolation requirements.

Pathogens covered
-----------------

    Target    Gene   Pathogen                               Clinical urgency
    MTB       rpoB   Mycobacterium tuberculosis             WHO/national notifiable
    CDIFF     tcdB   Clostridioides difficile (toxigenic)   HAI, contact precautions
    SARS2     N      SARS-CoV-2                             Isolation, antiviral
    FLU_A     M      Influenza A (all subtypes, pan-IAV)    Antiviral, isolation
    GAS       speB   Group A Streptococcus (S. pyogenes)    Antibiotic, stewardship

The assay designs are described in Recipes 1 (MTB/rpoB), 21 (GAS/speB),
22 (C.diff/tcdB), 2 (SARS-CoV-2/N), and 12 (Influenza A/M) of the cookbook.

Scoring model
-------------
Each pathogen is assigned a base weight reflecting public-health severity,
treatment urgency, and isolation requirements:

    MTB    : 50  WHO/national notifiable; 6-month RIPE treatment; highest
                 public-health consequence of a missed diagnosis
    CDIFF  : 25  HAI, CDI; contact precautions and specific antibiotic therapy
                 (vancomycin/fidaxomicin); ~5-10% 30-day mortality
    SARS2  : 20  Respiratory isolation required; antivirals (nirmatrelvir,
                 remdesivir) with narrow treatment window
    FLU_A  : 15  Antiviral therapy (oseltamivir) most effective within 48 h;
                 respiratory isolation and contact precautions
    GAS    :  10  Antibiotic stewardship impact; penicillin/amoxicillin; low
                  mortality but prevents rheumatic fever sequelae

Alert level boundaries
----------------------
    Score >= 50  -> CRITICAL   (MTB or very high combined burden)
    Score >= 25  -> HIGH       (CDiff, or multi-pathogen respiratory burden)
    Score >= 15  -> MODERATE   (single respiratory viral pathogen)
    Score >=  5  -> LOW        (GAS pharyngitis or low-severity single target)
    Score ==  0  -> NEGATIVE   (all targets negative)

Special clinical flags
----------------------
``immediate_report_required`` is True whenever MTB is positive.  In most
jurisdictions tuberculosis is a mandatorily notifiable disease; the national
public health authority must be notified before the patient leaves the clinic.

``contact_precautions_required`` is True when CDiff or MTB is positive:
CDiff requires contact precautions; MTB requires airborne/respiratory
isolation with contact precautions.

``antiviral_indicated`` is True when SARS-CoV-2 or Influenza A is positive,
signalling that treatment-window antivirals should be considered.

``antibiotic_indicated`` is True when GAS is positive, supporting the
antibiotic-prescribing decision for streptococcal pharyngitis.

References:
    WHO (2024) Global tuberculosis report. WHO/HTM/TB/2024.
        https://www.who.int/teams/global-tuberculosis-programme/tb-reports
    McDonald LC et al. (2018) Clinical practice guidelines for CDI, 2017
        update. Clin Infect Dis 66(7):e1-e48. doi:10.1093/cid/cix1085
    Shulman ST et al. (2012) GAS pharyngitis guidelines, IDSA 2012.
        Clin Infect Dis 55(10):e86-e102. doi:10.1093/cid/cis629
    WHO (2022) WHO guidelines for the pharmacological treatment of pandemic
        influenza A (H1N1) 2009 and other influenza viruses.
    Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
        Nucleic Acids Res 28(12):e63. doi:10.1093/nar/28.12.e63
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


class PocAlertLevel(StrEnum):
    """Human POC LAMP panel clinical alert severity.

    CRITICAL: MTB detected; immediate public-health notification and
        confirmatory testing required before patient is released.
    HIGH: Serious pathogen requiring same-day clinical action, isolation,
        and/or specific antibiotic/antiviral therapy.
    MODERATE: Actionable detection; antiviral or isolation guidance applies.
    LOW: Positive detection in a manageable target; clinical review and
        standard antibiotic or supportive care.
    NEGATIVE: All targets negative.  Confirm a sample-adequacy control
        (e.g. 16S LAMP) is positive before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class PocPanelFlags:
    """Five-target positivity flags for one human POC LAMP panel result.

    Each flag is True if the corresponding LAMP assay returned a positive
    amplification signal for that pathogen.

    Attributes:
        mtb: Mycobacterium tuberculosis (rpoB) positive.
        cdiff: Clostridioides difficile toxin B (tcdB) positive.
        sars2: SARS-CoV-2 (N gene) positive.
        flu_a: Influenza A virus (M gene, pan-IAV) positive.
        gas: Group A Streptococcus / S. pyogenes (speB) positive.
    """

    mtb: bool
    cdiff: bool
    sars2: bool
    flu_a: bool
    gas: bool


@dataclass(frozen=True, slots=True)
class PocRiskAssessment:
    """Structured clinical alert from a human POC LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        detected_targets: Labels of positive targets in canonical order.
        notifiable_targets: Subset of detected_targets requiring public-
            health notification (MTB in all jurisdictions; SARS-CoV-2
            in many).
        immediate_report_required: True when MTB is positive; tuberculosis
            is mandatorily notifiable in virtually all jurisdictions.
        contact_precautions_required: True when CDiff or MTB is positive.
        antiviral_indicated: True when SARS-CoV-2 or Influenza A positive.
        antibiotic_indicated: True when GAS is positive.
        interpretation: Concise clinical interpretation of the panel result.
        recommended_action: Actionable clinical recommendation for this
            detection pattern.
    """

    flags: PocPanelFlags
    alert_level: PocAlertLevel
    alert_score: int
    detected_targets: tuple[str, ...]
    notifiable_targets: tuple[str, ...]
    immediate_report_required: bool
    contact_precautions_required: bool
    antiviral_indicated: bool
    antibiotic_indicated: bool
    interpretation: str
    recommended_action: str

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; panel flags are inlined.
        """
        return {
            "alert_level": self.alert_level.value,
            "alert_score": self.alert_score,
            "detected_targets": list(self.detected_targets),
            "notifiable_targets": list(self.notifiable_targets),
            "immediate_report_required": self.immediate_report_required,
            "contact_precautions_required": self.contact_precautions_required,
            "antiviral_indicated": self.antiviral_indicated,
            "antibiotic_indicated": self.antibiotic_indicated,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "panel_flags": {
                "mtb": self.flags.mtb,
                "cdiff": self.flags.cdiff,
                "sars2": self.flags.sars2,
                "flu_a": self.flags.flu_a,
                "gas": self.flags.gas,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_MTB = 50
_WEIGHT_CDIFF = 25
_WEIGHT_SARS2 = 20
_WEIGHT_FLU_A = 15
_WEIGHT_GAS = 10

# Mandatorily notifiable in virtually all jurisdictions
_NOTIFIABLE = frozenset({"MTB"})

_TARGET_DESCRIPTION: dict[str, str] = {
    "MTB": "Mycobacterium tuberculosis (rpoB) -- WHO/national notifiable disease",
    "CDIFF": "Clostridioides difficile toxin B (tcdB) -- HAI, contact precautions required",
    "SARS2": "SARS-CoV-2 (N gene) -- respiratory isolation and antiviral treatment",
    "FLU_A": "Influenza A virus (M gene, pan-IAV) -- antiviral within 48 h, isolation",
    "GAS": "Group A Streptococcus (speB) -- antibiotic therapy, prevents sequelae",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> PocAlertLevel:
    if score >= 50:
        return PocAlertLevel.CRITICAL
    if score >= 25:
        return PocAlertLevel.HIGH
    if score >= 15:
        return PocAlertLevel.MODERATE
    if score >= 5:
        return PocAlertLevel.LOW
    return PocAlertLevel.NEGATIVE


def _build_text(
    flags: PocPanelFlags,
    detected_targets: tuple[str, ...],
    alert_level: PocAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        detected_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    mtb, cdiff, sars2, flu_a, gas = (
        flags.mtb,
        flags.cdiff,
        flags.sars2,
        flags.flu_a,
        flags.gas,
    )

    # All negative
    if not any([mtb, cdiff, sars2, flu_a, gas]):
        return (
            "All POC panel targets negative. No pathogen detected.",
            "Confirm a sample-adequacy control (e.g. 16S LAMP) is positive before "
            "reporting a clean result. If clinical suspicion remains, consider "
            "repeat testing or alternative diagnostics.",
        )

    # MTB positive -- always CRITICAL regardless of co-detections
    if mtb:
        co = [t for t in ("CDIFF", "SARS2", "FLU_A", "GAS") if getattr(flags, t.lower())]
        co_str = ", ".join(co) if co else ""
        interp = (
            "Mycobacterium tuberculosis (rpoB LAMP) detected. TB is a mandatorily "
            "notifiable disease in virtually all jurisdictions. Confirmatory testing "
            "(sputum smear, culture, Xpert MTB/RIF) is required before initiating "
            "RIPE treatment. The rpoB target also co-localises with rifampicin-"
            "resistance mutations; an rpoB LAMP positive should be followed by "
            "drug-susceptibility testing."
        )
        if co_str:
            interp += f" Panel co-detections: {co_str}."
        action = (
            "IMMEDIATE: Place patient in airborne isolation (negative-pressure room "
            "or equivalent; N95/FFP2 respirators for HCW). Notify the national "
            "public health authority per local TB notification protocol. Collect "
            "three sputum specimens for smear microscopy and culture. Do not start "
            "RIPE therapy until confirmatory diagnostics are in progress."
        )
        return interp, action

    # CDiff positive (no MTB)
    if cdiff:
        co_resp = [t for t in ("SARS2", "FLU_A") if getattr(flags, t.lower())]
        co_str = ", ".join(co_resp) if co_resp else ""
        interp = (
            "Toxigenic Clostridioides difficile (tcdB LAMP) detected. CDI is the "
            "leading cause of hospital-acquired infectious diarrhoea. The tcdB "
            "target detects all toxigenic strains including the A-B+ phenotype "
            "(ribotypes RT017, RT033) that a tcdA-only assay would miss."
        )
        if gas:
            interp += (
                " GAS also detected: dual infection is unusual; review clinical "
                "presentation to confirm both results."
            )
        if co_str:
            interp += f" Respiratory co-detections: {co_str}."
        action = (
            "Implement contact precautions immediately (gloves, gown, single "
            "room or cohort isolation; soap-and-water hand hygiene -- alcohol "
            "gel does not eliminate C. diff spores). Initiate therapy per local "
            "CDI guideline (vancomycin PO or fidaxomicin for non-severe CDI; "
            "vancomycin PO for severe; consider bezlotoxumab for recurrence "
            "prevention). Avoid proton-pump inhibitors and unnecessary antibiotics. "
            "Do not use LAMP for test-of-cure: tcdB can remain positive for "
            "weeks after clinical resolution."
        )
        return interp, action

    # Respiratory co-detection: SARS-CoV-2 + Flu A (no MTB or CDiff)
    if sars2 and flu_a:
        interp = (
            "SARS-CoV-2 (N gene) and Influenza A (M gene) co-detected. Dual "
            "respiratory virus infection (flurona) is clinically significant: "
            "both antivirals may be warranted and respiratory isolation applies "
            "for both pathogens simultaneously."
        )
        if gas:
            interp += (
                " GAS also detected: bacterial pharyngitis co-infection; "
                "antibiotic therapy for GAS is additionally indicated."
            )
        action = (
            "Respiratory isolation (droplet + contact precautions). Consider "
            "antivirals for both SARS-CoV-2 (nirmatrelvir/ritonavir or remdesivir "
            "within 5 days of symptom onset for high-risk patients) and Influenza A "
            "(oseltamivir within 48 h of symptom onset). Report to local public "
            "health surveillance per applicable protocol."
        )
        return interp, action

    # SARS-CoV-2 alone (no MTB, CDiff, or Flu A)
    if sars2:
        interp = (
            "SARS-CoV-2 (N gene LAMP) detected. Respiratory isolation is required. "
            "Antiviral treatment (nirmatrelvir/ritonavir, remdesivir, or molnupiravir) "
            "should be considered for patients at high risk of progression to severe "
            "disease; the treatment window is within 5 days of symptom onset."
        )
        if gas:
            interp += (
                " GAS also detected: bacterial pharyngitis co-infection; "
                "antibiotic therapy for GAS is additionally indicated."
            )
        action = (
            "Respiratory and contact isolation. Assess eligibility for outpatient "
            "antivirals per local COVID-19 treatment protocol. Notify local "
            "surveillance if variant-of-concern guidance requires reporting. "
            "Advise patient on isolation duration per current public-health guidance."
        )
        return interp, action

    # Flu A alone (no MTB, CDiff, or SARS-CoV-2)
    if flu_a:
        interp = (
            "Influenza A (M-gene LAMP, pan-IAV) detected. This assay detects all "
            "influenza A subtypes (H1N1, H3N2, H5N1, etc.) without subtype "
            "discrimination. In a patient with poultry or swine exposure, "
            "zoonotic influenza (HPAI H5N1, H1N2v) should be considered and "
            "reported to public health authorities."
        )
        if gas:
            interp += (
                " GAS co-detected: viral-bacterial co-infection; antibiotic "
                "therapy for GAS is additionally indicated."
            )
        action = (
            "Respiratory and contact precautions. Initiate oseltamivir within "
            "48 h of symptom onset (75 mg twice daily for 5 days) for high-risk "
            "patients or those with severe disease. If avian/swine exposure is "
            "reported, notify public health immediately for zoonotic influenza "
            "investigation."
        )
        return interp, action

    # GAS alone (only non-respiratory, non-notifiable target)
    if gas and not any([mtb, cdiff, sars2, flu_a]):
        interp = (
            "Group A Streptococcus (speB LAMP) detected. GAS pharyngitis is the "
            "most common bacterial cause of sore throat and the primary indication "
            "for antibiotic prescribing in primary care. Treatment prevents "
            "rheumatic fever, rheumatic heart disease, and post-streptococcal "
            "glomerulonephritis sequelae."
        )
        action = (
            "Prescribe antibiotics per local GAS pharyngitis guideline "
            "(penicillin V or amoxicillin 10-day course; erythromycin or "
            "cefalexin for penicillin allergy). Advise patient on symptom "
            "duration, return precautions, and antibiotic adherence. No "
            "routine test-of-cure is required for uncomplicated GAS pharyngitis."
        )
        return interp, action

    # Fallback for any other combination
    targets_str = ", ".join(detected_targets)
    interp = f"Detected targets: {targets_str}. {alert_level.value} alert level."
    action = (
        "Review detection pattern with the treating clinician. Follow applicable "
        "clinical guidelines for each detected pathogen."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_poc_risk(flags: PocPanelFlags) -> PocRiskAssessment:
    """Compute a structured clinical alert from five-target human POC panel flags.

    Weights each pathogen by public-health severity, treatment urgency, and
    isolation requirements, sums to a 0-100 score, and maps to one of five
    categorical alert levels.  Sets clinical action flags for isolation,
    antiviral, and antibiotic guidance.

    Args:
        flags: Positivity flags for all five human POC panel targets.

    Returns:
        :class:`PocRiskAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.poc_risk import PocPanelFlags, assess_poc_risk

        flags = PocPanelFlags(
            mtb=False, cdiff=True, sars2=False, flu_a=False, gas=False
        )
        assessment = assess_poc_risk(flags)
        print(assessment.alert_level)              # PocAlertLevel.HIGH
        print(assessment.contact_precautions_required)  # True
    """
    score = 0
    if flags.mtb:
        score += _WEIGHT_MTB
    if flags.cdiff:
        score += _WEIGHT_CDIFF
    if flags.sars2:
        score += _WEIGHT_SARS2
    if flags.flu_a:
        score += _WEIGHT_FLU_A
    if flags.gas:
        score += _WEIGHT_GAS

    score = max(0, min(100, score))

    target_order = [
        ("MTB", flags.mtb),
        ("CDIFF", flags.cdiff),
        ("SARS2", flags.sars2),
        ("FLU_A", flags.flu_a),
        ("GAS", flags.gas),
    ]
    detected_targets = tuple(label for label, pos in target_order if pos)
    notifiable_targets = tuple(t for t in detected_targets if t in _NOTIFIABLE)

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, detected_targets, alert_level)

    return PocRiskAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        detected_targets=detected_targets,
        notifiable_targets=notifiable_targets,
        immediate_report_required=flags.mtb,
        contact_precautions_required=flags.cdiff or flags.mtb,
        antiviral_indicated=flags.sars2 or flags.flu_a,
        antibiotic_indicated=flags.gas,
        interpretation=interpretation,
        recommended_action=recommended_action,
    )


def flags_from_dict(data: dict[str, object]) -> PocPanelFlags:
    """Construct :class:`PocPanelFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``mtb``, ``cdiff``, ``sars2``, ``flu_a``, ``gas`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`PocPanelFlags` instance.
    """
    return PocPanelFlags(
        mtb=bool(data.get("mtb", False)),
        cdiff=bool(data.get("cdiff", False)),
        sars2=bool(data.get("sars2", False)),
        flu_a=bool(data.get("flu_a", False)),
        gas=bool(data.get("gas", False)),
    )


def write_assessment_json(assessment: PocRiskAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`PocRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: PocRiskAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`PocRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        d = assessment.to_dict()
        for key, val in d.items():
            if key in ("detected_targets", "notifiable_targets"):
                writer.writerow([key, "; ".join(val)])  # type: ignore[arg-type]
            elif key == "panel_flags" and isinstance(val, dict):
                for pk, pv in val.items():
                    writer.writerow([f"target_{pk}", "positive" if pv else "negative"])
            else:
                writer.writerow([key, val])
