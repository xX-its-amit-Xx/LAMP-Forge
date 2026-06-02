"""Swine Respiratory Disease Complex (PRDC) LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable swine respiratory
LAMP panel as a structured alert, prioritising pathogens by their individual
severity and their synergistic effect in the Porcine Respiratory Disease
Complex (PRDC).

PRDC is not a single-pathogen disease -- it is a multi-pathogen syndrome in
which two or three co-infecting agents potentiate each other.  PRRSV
immunosuppresses the alveolar macrophages that normally clear secondary
pathogens; PCV2 further dysregulates innate immunity; Mhp destroys the
mucociliary escalator and provides attachment sites for Pasteurella multocida
and Actinobacillus pleuropneumoniae.  The combination of all three pathogens
in a growing pig produces mortality and culling rates far exceeding what any
single pathogen causes alone.

Pathogens covered
-----------------

    Target  Gene   Pathogen                              NA type  Economic impact
    PRRSV   ORF7   Porcine reprod. & resp. syndrome      RNA      USD 664 M/yr (US)
    PCV2    cap    Porcine circovirus 2 (PCVAD)          DNA      EUR 385 M/yr (EU)
    Mhp     p97    Mycoplasma hyopneumoniae (SEP)        DNA      USD >1 B/yr (global)

Scoring model
-------------
Base weights reflect individual economic impact and the severity of the
immunosuppressive contribution each pathogen makes to the PRDC complex:

    PRRSV : 40  The most costly endemic swine pathogen in North America and
                Europe; USD 664 M/yr US losses (Holtkamp et al. 2013).  Causes
                reproductive failure in sows (abortion storms, "blue ear")
                and severe respiratory disease in nursery and grow-finish pigs.
                Destroys alveolar macrophages, enabling all secondary pathogens.
                PRRSV alone justifies immediate veterinary assessment.
    PCV2  : 30  Ubiquitous ssDNA virus; the primary cause of PCVAD (PMWS,
                PDNS).  EUR 385 M/yr EU losses (Segalés et al. 2005).  Acts
                synergistically with PRRSV; co-infection multiplies lung lesion
                scores and mortality.  The "immunosuppressive core" of complex
                PRDC.
    Mhp   : 20  Aetiological agent of Swine Enzootic Pneumonia (SEP); nearly
                universal in commercial herds.  Direct losses via ADG reduction
                (~15%) and FCR deterioration (~10%); indirect losses via
                secondary bacterial co-infections.  The "gateway" pathogen:
                destroys mucociliary escalator and provides Pasteurella /
                Actinobacillus attachment sites.

Alert level boundaries
----------------------
    Score >= 60  -> CRITICAL   (PRRSV + at least one co-pathogen: full PRDC)
    Score >= 40  -> HIGH       (PRRSV alone; or PCV2 + Mhp dual infection)
    Score >= 30  -> MODERATE   (PCV2 alone; requires veterinary assessment)
    Score >=  1  -> LOW        (Mhp alone; endemic management protocol)
    Score ==  0  -> NEGATIVE   (all negative; confirm extraction control)

Special decision flags
----------------------
``movement_restriction_required`` is set True when PRRSV is positive.  PRRSV
is a highly contagious RNA virus spread via aerosol, direct pig contact, and
fomites.  Pigs in affected units must not be moved to negative premises until
the herd has completed a monitored PRRSV-elimination or PRRSV-stabilisation
protocol.

``complex_prdc`` is set True when two or more targets are positive
simultaneously.  The multi-pathogen synergy in complex PRDC substantially
increases mortality risk and antibiotic treatment requirements beyond what
individual-pathogen models predict.  A combined veterinary and production
team review is warranted before making any treatment or culling decisions.

References:
    Holtkamp DJ et al. (2013) Assessment of the economic impact of PRRSV
        on United States pork producers. J Swine Health Prod 21(2):72-84.
    Segalés J, Allan GM, Domingo M (2005) Porcine circovirus diseases.
        Anim Health Res Rev 6(2):119-142. doi:10.1079/AHR2005106
    Maes D et al. (2018) Swine enzootic pneumonia: a review. Vet Rec Open
        5(1):e000228. doi:10.1136/vetreco-2017-000228
    Opriessnig T et al. (2011) Porcine circovirus type 2 and coinfecting
        pathogens in field cases: a retrospective study. J Vet Diagn Invest
        23(3):449-457. doi:10.1177/1040638711407356
    Harms PA et al. (2002) Experimental reproduction of PRDC with PRRSV,
        Mhp, and PCV2. J Swine Health Prod 10(3):105-111.
    Zhang X et al. (2019) Loop-mediated isothermal amplification for
        Mycoplasma hyopneumoniae detection. J Vet Diagn Invest
        31(4):571-576. doi:10.1177/1040638719843372
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


class PrdcAlertLevel(StrEnum):
    """Swine respiratory / PRDC LAMP panel alert severity.

    Attributes:
        CRITICAL: PRRSV detected alongside PCV2 and/or Mhp -- full complex
            PRDC presentation.  Mortality and production losses are
            synergistically elevated; immediate veterinary management
            and movement restriction required.
        HIGH: PRRSV detected alone; or PCV2 and Mhp detected together
            without PRRSV.  Either scenario warrants same-day veterinary
            contact and movement restriction evaluation.
        MODERATE: PCV2 detected alone; PCVAD risk requires veterinary
            assessment.  Monitor for emerging PRRSV or Mhp co-infection.
        LOW: Mhp detected alone; endemic enzootic pneumonia management
            protocol.  Review metaphylaxis and vaccination schedule.
        NEGATIVE: All three targets negative.  Confirm extraction control
            (e.g. 16S rRNA) is positive before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class PrdcFlags:
    """Positivity flags for one swine respiratory / PRDC LAMP panel result.

    Each flag is True if the corresponding LAMP or RT-LAMP assay returned a
    positive amplification signal for that pathogen.

    Attributes:
        prrsv: PRRSV ORF7 (nucleocapsid N; RT-LAMP, RNA target) positive.
        pcv2: PCV2 cap gene (ORF2; standard LAMP, DNA target) positive.
        mhp: Mycoplasma hyopneumoniae p97 adhesin (standard LAMP, DNA
            target) positive.
    """

    prrsv: bool
    pcv2: bool
    mhp: bool


@dataclass(frozen=True, slots=True)
class PrdcAssessment:
    """Structured alert from a swine respiratory / PRDC LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation.
        movement_restriction_required: True when PRRSV is positive.  PRRSV
            spreads by aerosol and direct contact; pigs must not leave the
            affected unit to negative premises until elimination or
            stabilisation criteria are met.
        complex_prdc: True when two or more targets are simultaneously
            positive.  Multi-pathogen synergy substantially increases
            mortality and production losses beyond single-pathogen models.
    """

    flags: PrdcFlags
    alert_level: PrdcAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    interpretation: str
    recommended_action: str
    movement_restriction_required: bool
    complex_prdc: bool

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; panel flags are inlined under
            the ``panel_flags`` key.
        """
        return {
            "alert_level": self.alert_level.value,
            "alert_score": self.alert_score,
            "active_targets": list(self.active_targets),
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "movement_restriction_required": self.movement_restriction_required,
            "complex_prdc": self.complex_prdc,
            "panel_flags": {
                "prrsv": self.flags.prrsv,
                "pcv2": self.flags.pcv2,
                "mhp": self.flags.mhp,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_PRRSV = 40
_WEIGHT_PCV2 = 30
_WEIGHT_MHP = 20


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> PrdcAlertLevel:
    if score >= 60:
        return PrdcAlertLevel.CRITICAL
    if score >= 40:
        return PrdcAlertLevel.HIGH
    if score >= 30:
        return PrdcAlertLevel.MODERATE
    if score >= 1:
        return PrdcAlertLevel.LOW
    return PrdcAlertLevel.NEGATIVE


def _build_text(
    flags: PrdcFlags,
    active_targets: tuple[str, ...],
    alert_level: PrdcAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    prrsv, pcv2, mhp = flags.prrsv, flags.pcv2, flags.mhp

    if not any([prrsv, pcv2, mhp]):
        return (
            "All swine respiratory / PRDC targets negative. No PRRSV, PCV2, "
            "or Mycoplasma hyopneumoniae detected.",
            "Confirm an internal extraction control (e.g. 16S rRNA) is positive "
            "before reporting a clean result. Continue routine surveillance. Re-test "
            "if coughing, dyspnoea, or mortality above baseline is observed in "
            "nursery or grow-finish pigs.",
        )

    # All three -- full complex PRDC
    if prrsv and pcv2 and mhp:
        return (
            "Full Porcine Respiratory Disease Complex (PRDC) detected: PRRSV ORF7, "
            "PCV2 cap gene, and Mycoplasma hyopneumoniae p97 all positive. This "
            "triple-positive pattern is the classic full PRDC presentation. PRRSV "
            "destroys alveolar macrophages; PCV2 further dysregulates innate immunity; "
            "Mhp eliminates the mucociliary escalator. Combined mortality and culling "
            "rates in triple-positive grow-finish pigs substantially exceed the sum of "
            "individual pathogen effects. Rapid deterioration in the coming 5-14 days "
            "is expected without intervention.",
            "CRITICAL: Immediately notify the site veterinarian and production manager. "
            "Restrict all pig movement from the affected unit. Implement full PRDC "
            "management protocol: evaluate enrofloxacin or doxycycline metaphylaxis for "
            "Mhp, assess PRRSV vaccination status and consider modified-live vaccine "
            "exposure, and review PCV2 vaccination compliance. Increase pen-by-pen "
            "monitoring to detect early mortality. Contact the herd health veterinarian "
            "for a whole-herd respiratory health assessment. Submit additional samples "
            "to an accredited laboratory for secondary pathogen culture (Pasteurella "
            "multocida, Actinobacillus pleuropneumoniae, Streptococcus suis) to guide "
            "antibiotic selection.",
        )

    # PRRSV + PCV2 (no Mhp)
    if prrsv and pcv2 and not mhp:
        return (
            "PRRSV and PCV2 co-detected (Mhp negative). PRRSV-PCV2 co-infection is "
            "the most common two-pathogen PRDC presentation and carries a substantially "
            "elevated risk of PCVAD (wasting, lymphadenopathy, interstitial pneumonia). "
            "The combination dysregulates both innate and adaptive immunity, enabling "
            "secondary bacterial pneumonia even without Mhp.",
            "HIGH: Notify the site veterinarian. Restrict pig movement from the affected "
            "unit until PRRSV status is clarified. Assess PCV2 vaccination compliance in "
            "the affected cohort; emergency booster vaccination may be indicated for "
            "unvaccinated pigs. Evaluate PRRSV herd-level status and consider whole-herd "
            "exposure or MLV vaccination protocol under veterinary supervision. Monitor "
            "for Mhp emergence and secondary bacterial infection over the next 2-4 weeks.",
        )

    # PRRSV + Mhp (no PCV2)
    if prrsv and mhp and not pcv2:
        return (
            "PRRSV and Mycoplasma hyopneumoniae co-detected (PCV2 negative). PRRSV-Mhp "
            "co-infection potentiates respiratory disease: PRRSV impairs alveolar "
            "macrophage function while Mhp destroys the mucociliary escalator, together "
            "creating conditions for severe secondary bacterial pneumonia. Lung lesion "
            "scores in co-infected pigs significantly exceed either pathogen alone.",
            "HIGH: Notify the site veterinarian. Restrict movement from the affected "
            "unit. Evaluate antibiotic metaphylaxis targeting Mhp (enrofloxacin, "
            "doxycycline, or tiamulin; note Mhp is intrinsically resistant to "
            "beta-lactams). Assess PRRSV vaccination status and implement a PRRSV "
            "management protocol. Monitor daily for worsening dyspnoea, increased "
            "mortality, or clinical signs consistent with Actinobacillus "
            "pleuropneumoniae co-infection. Re-test in 2 weeks to assess whether "
            "PCV2 has emerged.",
        )

    # PCV2 + Mhp (no PRRSV)
    if pcv2 and mhp and not prrsv:
        return (
            "PCV2 and Mycoplasma hyopneumoniae co-detected (PRRSV negative). PCV2-Mhp "
            "co-infection without PRRSV is a common endemic-herd scenario. Both pathogens "
            "contribute to respiratory disease; PCV2 may cause PCVAD (wasting, "
            "lymphadenopathy) in unvaccinated cohorts, and Mhp degrades mucociliary "
            "clearance. The absence of PRRSV is favourable, but the dual detection "
            "warrants veterinary attention and vaccination review.",
            "HIGH: Notify the site veterinarian. Review PCV2 vaccination compliance -- "
            "emergency booster vaccination may be indicated for unvaccinated pigs. "
            "Implement antibiotic metaphylaxis targeting Mhp under veterinary "
            "prescription. Monitor the herd for PRRSV emergence; PRRSV in this context "
            "would escalate to CRITICAL. Continue routine PRRSV surveillance.",
        )

    # PRRSV alone
    if prrsv and not pcv2 and not mhp:
        return (
            "PRRSV detected via the ORF7 assay (PCV2 and Mhp negative). PRRSV alone "
            "causes significant reproductive failure in sows and severe respiratory "
            "disease in nursery and grow-finish pigs (USD 664 M/yr US economic impact). "
            "A PRRSV-positive grow-finish unit without PCV2 or Mhp co-infection is a "
            "more manageable situation, but PRRSV immunosuppression makes the unit highly "
            "vulnerable to PCV2 and secondary bacterial respiratory pathogens. The panel "
            "should be repeated within 2-4 weeks to detect emerging co-infections.",
            "HIGH: Notify the site veterinarian immediately. Restrict pig movement from "
            "the affected unit to PRRSV-negative premises. Assess PRRSV vaccination "
            "status of the affected cohort and the breeding herd. Consider whole-herd "
            "modified-live vaccine exposure (MLV) under veterinary supervision. Implement "
            "enhanced biosecurity (nose-to-nose contact prevention, air filtration if "
            "available). Re-test with the full PRDC panel in 2-4 weeks.",
        )

    # PCV2 alone
    if pcv2 and not prrsv and not mhp:
        return (
            "PCV2 detected via the cap-gene assay (PRRSV and Mhp negative). PCV2 alone "
            "causes PCVAD -- including PMWS (wasting, lymphadenopathy, jaundice) and "
            "interstitial pneumonia -- particularly in unvaccinated or incompletely "
            "vaccinated pig cohorts. The absence of PRRSV reduces the risk of acute "
            "complex PRDC, but PCV2 alone can cause significant mortality (10-15% in "
            "unvaccinated PMWS outbreaks) and represents the immunosuppressive substrate "
            "for PRRSV co-infection if PRRSV is introduced.",
            "MODERATE: Notify the site veterinarian. Review PCV2 vaccination compliance "
            "in all age groups and consider emergency booster vaccination for unvaccinated "
            "or under-vaccinated cohorts. Monitor for PMWS clinical signs (wasting, "
            "jaundice, enlarged lymph nodes). Verify that PRRSV status remains negative "
            "by repeating the full PRDC panel within 4 weeks. Maintain biosecurity to "
            "prevent PRRSV introduction.",
        )

    # Mhp alone
    if mhp and not prrsv and not pcv2:
        return (
            "Mycoplasma hyopneumoniae (Mhp) detected via the p97 assay (PRRSV and PCV2 "
            "negative). Mhp alone causes Swine Enzootic Pneumonia (SEP), characterised "
            "by chronic dry coughing, reduced ADG (~15%), and poor FCR. Nearly all "
            "commercial herds are endemically infected; an Mhp-positive result without "
            "PRRSV or PCV2 indicates a manageable endemic situation rather than an "
            "acute disease crisis. However, Mhp destroys the mucociliary escalator and "
            "predisposes pigs to secondary bacterial pneumonia.",
            "LOW: Review the metaphylaxis and vaccination programme for Mhp with the "
            "site veterinarian. Antibiotic metaphylaxis (enrofloxacin, doxycycline, or "
            "tiamulin; NOT beta-lactams -- Mhp is intrinsically resistant) at the "
            "expected colonisation window (typically 4-6 weeks post-weaning) reduces "
            "lung lesion scores. Confirm PCV2 and PRRSV remain negative; emergence of "
            "either would escalate the clinical situation significantly. Consider the "
            "all-in/all-out management to reduce cumulative Mhp challenge pressure.",
        )

    # Fallback
    targets_str = ", ".join(active_targets)
    return (
        f"Detected: {targets_str}. {alert_level.value} alert.",
        "Review with site veterinarian. Follow standard respiratory disease protocol.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_prdc_risk(flags: PrdcFlags) -> PrdcAssessment:
    """Compute a structured PRDC alert from three-target panel flags.

    Weights each pathogen by its individual economic impact and its
    contribution to multi-pathogen synergy in the Porcine Respiratory Disease
    Complex, sums to a 0-100 score, and maps it to one of five categorical
    alert levels.

    Sets ``movement_restriction_required`` when PRRSV is positive and
    ``complex_prdc`` when two or more targets are simultaneously positive.

    Args:
        flags: Positivity flags for all three swine respiratory targets.

    Returns:
        :class:`PrdcAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.swine_respiratory_risk import PrdcFlags, assess_prdc_risk

        flags = PrdcFlags(prrsv=True, pcv2=True, mhp=False)
        assessment = assess_prdc_risk(flags)
        print(assessment.alert_level)     # PrdcAlertLevel.CRITICAL
        print(assessment.complex_prdc)    # True
    """
    score = 0
    if flags.prrsv:
        score += _WEIGHT_PRRSV
    if flags.pcv2:
        score += _WEIGHT_PCV2
    if flags.mhp:
        score += _WEIGHT_MHP

    score = max(0, min(100, score))

    target_order = [
        ("PRRSV", flags.prrsv),
        ("PCV2", flags.pcv2),
        ("Mhp", flags.mhp),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    n_positive = sum(1 for _, pos in target_order if pos)
    complex_prdc = n_positive >= 2

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, active_targets, alert_level)

    return PrdcAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        interpretation=interpretation,
        recommended_action=recommended_action,
        movement_restriction_required=flags.prrsv,
        complex_prdc=complex_prdc,
    )


def flags_from_dict(data: dict[str, object]) -> PrdcFlags:
    """Construct :class:`PrdcFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``prrsv``, ``pcv2``, ``mhp`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`PrdcFlags` instance.
    """
    return PrdcFlags(
        prrsv=bool(data.get("prrsv", False)),
        pcv2=bool(data.get("pcv2", False)),
        mhp=bool(data.get("mhp", False)),
    )


def write_assessment_json(assessment: PrdcAssessment, path: Path) -> None:
    """Write the PRDC assessment to a JSON file.

    Args:
        assessment: :class:`PrdcAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: PrdcAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the PRDC assessment.

    Args:
        assessment: :class:`PrdcAssessment` to serialise.
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
                writer.writerow([key, str(val)])
