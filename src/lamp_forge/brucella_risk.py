"""Brucella IS711 LAMP risk assessment engine.

Interprets a pan-Brucella IS711 LAMP result as a structured alert,
contextualised by the sample type (livestock species or human).

IS711 (also annotated IS6501 / ISBm2) is a Brucella-specific insertion
sequence present in 5-40 copies per genome across all nine Brucella species.
Its multi-copy nature gives exceptional LAMP sensitivity -- 10-100x lower
analytical LOD than single-copy targets -- making it the gold-standard
molecular target for brucellosis diagnostics in field-deployable settings
(Abd El Wahed et al. 2015 PLoS NTD; Al-Dahouk et al. 2010 Ann Clin Microbiol).

A single IS711 LAMP assay detects any Brucella species (B. abortus, B.
melitensis, B. suis, B. ovis, B. canis, B. neotomae, and the newer marine
species) but cannot discriminate between species.  The alert level and
recommended action therefore depend critically on the *context* -- which
animal species (or human) the sample came from -- because species identity
determines:

  - Probable Brucella species (B. abortus in cattle, B. melitensis in small
    ruminants, B. suis in swine, and all three in human brucellosis).
  - Regulatory reporting pathway (WOAH-listed for all livestock species; public
    health notification for human cases).
  - Zoonotic risk tier (B. melitensis poses the highest human risk globally).
  - Recommended clinical or veterinary action.

For environmental samples (water, soil, fomite wipes), a positive IS711
result signals recent Brucella contamination but does not confirm active
infection; follow-up animal testing is required before regulatory reporting.

BioVind vertical applicability
-------------------------------
Farm (animal biosecurity):
    Bovine brucellosis (B. abortus) is WOAH-listed and covered by mandatory
    national eradication programmes in the US, EU, and many other
    jurisdictions.  B. melitensis in small ruminants is the most significant
    cause of human brucellosis globally.  B. suis biovars 1/3 are zoonotic
    and economically significant in swine herds.  A positive IS711 result
    from a livestock sample triggers immediate herd quarantine and reporting
    to the national veterinary authority.

Human POC:
    Human brucellosis (undulant fever) is a febrile illness requiring
    doxycycline + rifampin combination therapy.  LAMP on blood or serum
    is a fast, field-appropriate screening tool ahead of confirmatory
    serology or culture.

References:
    Abd El Wahed A et al. (2015) Reverse transcription recombinase polymerase
        amplification assay for the detection of Middle East respiratory
        syndrome coronavirus. PLoS Curr 7. doi:10.1371/currents.outbreaks
    Abd El Wahed A et al. (2015) Loop-mediated isothermal amplification on a
        chip for detection of Brucella spp. in milk samples.
        PLoS Negl Trop Dis 9(10):e0004124. doi:10.1371/journal.pntd.0004124
    Al-Dahouk S et al. (2010) Evaluation of Brucella IS711-based loop-mediated
        isothermal amplification for human brucellosis.
        Ann Clin Microbiol Antimicrob 9:15. doi:10.1186/1476-0711-9-15
    OIE/WOAH (2021). Bovine brucellosis. Terrestrial Manual, Chapter 3.1.2.
    Godfroid J et al. (2011). Brucellosis at the animal/ecosystem/human
        interface at the beginning of the 21st century.
        Prev Vet Med 102(2):118-131. doi:10.1016/j.prevetmed.2011.04.007
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


class BrucellaContext(StrEnum):
    """Sample context for IS711 LAMP result interpretation.

    Attributes:
        CATTLE: Blood, milk, lymph node, or tissue from cattle or bison.
            Probable species: Brucella abortus.
        SMALL_RUMINANT: Sample from sheep or goats.
            Probable species: Brucella melitensis (highest human-risk species).
        SWINE: Sample from domestic or feral pigs.
            Probable species: Brucella suis (biovars 1/3 are zoonotic).
        HUMAN: Blood, serum, or tissue from a human patient.
            Undulant fever presentation; B. melitensis and B. abortus most
            common causes of human brucellosis globally.
        ENVIRONMENTAL: Water, soil, fomite wipe, or aborted fetal material.
            Positive IS711 indicates Brucella DNA contamination; active
            infection in the source animal requires follow-up testing.
    """

    CATTLE = "cattle"
    SMALL_RUMINANT = "small_ruminant"
    SWINE = "swine"
    HUMAN = "human"
    ENVIRONMENTAL = "environmental"


class BrucellaAlertLevel(StrEnum):
    """Categorical alert level for an IS711 LAMP result.

    Attributes:
        CRITICAL: Positive result in a human sample.  Brucellosis is a
            notifiable zoonosis requiring immediate antibiotic therapy and
            public health reporting.
        HIGH: Positive result in a livestock (cattle, small ruminant, or swine)
            sample.  WOAH-listed disease; mandatory herd quarantine and
            veterinary authority reporting required.
        LOW: Positive result in an environmental sample.  Indicates recent
            Brucella DNA contamination; follow-up animal testing required
            before regulatory reporting.
        NEGATIVE: IS711 LAMP negative.  Brucella not detected.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


_ALERT_SCORE: dict[BrucellaAlertLevel, int] = {
    BrucellaAlertLevel.CRITICAL: 90,
    BrucellaAlertLevel.HIGH: 60,
    BrucellaAlertLevel.LOW: 25,
    BrucellaAlertLevel.NEGATIVE: 0,
}


@dataclass(frozen=True, slots=True)
class BrucellaFlags:
    """IS711 LAMP result and sample context for one Brucella assessment.

    Attributes:
        is711_positive: True if the IS711 LAMP assay returned a positive
            amplification signal.
        context: Sample type and likely Brucella species context.  Defaults
            to CATTLE if not specified.
    """

    is711_positive: bool
    context: BrucellaContext = BrucellaContext.CATTLE


@dataclass(frozen=True, slots=True)
class BrucellaAssessment:
    """Structured alert from a Brucella IS711 LAMP result.

    Attributes:
        flags: IS711 result and context flags used to generate this assessment.
        alert_level: Categorical alert level (CRITICAL / HIGH / LOW / NEGATIVE).
        alert_score: Numeric severity, 0-100 (higher is more urgent).
        interpretation: Explanation of the detection result and its significance.
        recommended_action: Field-actionable recommendation.
        notifiable: True when a positive result triggers a mandatory regulatory
            report (livestock contexts; WOAH-listed bovine, caprine/ovine, and
            swine brucellosis; public health for human cases).
        zoonotic_risk: True when the context implies direct human exposure risk
            (livestock handling, veterinary procedures, or human POC).
        immediate_report_required: True when positive result requires same-day
            notification of the national veterinary or public health authority.
    """

    flags: BrucellaFlags
    alert_level: BrucellaAlertLevel
    alert_score: int
    interpretation: str
    recommended_action: str
    notifiable: bool
    zoonotic_risk: bool
    immediate_report_required: bool

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; flags are inlined under
            ``is711_result`` and ``context``.
        """
        return {
            "alert_level": self.alert_level.value,
            "alert_score": self.alert_score,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "notifiable": self.notifiable,
            "zoonotic_risk": self.zoonotic_risk,
            "immediate_report_required": self.immediate_report_required,
            "is711_result": "positive" if self.flags.is711_positive else "negative",
            "context": self.flags.context.value,
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _build_text(flags: BrucellaFlags) -> tuple[str, str]:
    """Return (interpretation, recommended_action) for the IS711 result.

    Args:
        flags: IS711 result and context.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    if not flags.is711_positive:
        return (
            "IS711 LAMP negative: Brucella genus not detected in this sample. "
            "IS711 is a multi-copy insertion sequence with exceptional analytical "
            "sensitivity; a negative result is strong evidence against active Brucella "
            "infection in the tested material.",
            "No immediate action required. Confirm an extraction control "
            "(e.g. 16S rRNA internal control) is positive before reporting a clean "
            "result. Re-test if clinical signs consistent with brucellosis are observed.",
        )

    ctx = flags.context

    if ctx is BrucellaContext.HUMAN:
        interp = (
            "IS711 LAMP positive: Brucella genus detected in a human sample. "
            "Human brucellosis (undulant fever) is a notifiable zoonotic disease. "
            "IS711 is present in all nine Brucella species; probable species depends "
            "on occupational and animal contact history "
            "(B. melitensis from small-ruminant contact is the most common cause "
            "of human brucellosis globally). "
            "Confirmatory serology (Rose Bengal plate test + standard tube agglutination "
            "or Wright titre) and/or blood culture at a BSL-2/3 reference laboratory "
            "is required before initiating therapy."
        )
        action = (
            "Initiate doxycycline (100 mg twice daily) + rifampin (600-900 mg once daily) "
            "combination therapy for 6 weeks after collection of confirmatory samples. "
            "Report to the national public health authority per local notifiable-disease "
            "regulations. "
            "Interview patient for animal contact history (livestock handling, "
            "consumption of unpasteurised dairy products, veterinary/abattoir work) "
            "to identify the exposure source and guide public health investigation. "
            "Do not delay therapy pending culture results; culture sensitivity for "
            "Brucella from blood is only 50-80%."
        )
        return interp, action

    if ctx is BrucellaContext.CATTLE:
        interp = (
            "IS711 LAMP positive: Brucella genus detected in a cattle or bison sample. "
            "B. abortus (bovine brucellosis) is the most likely species. "
            "Bovine brucellosis is a WOAH-listed disease subject to mandatory national "
            "eradication programmes in the US (USDA APHIS), EU, and most other "
            "jurisdictions. "
            "All farm workers and veterinarians with recent animal contact are at "
            "zoonotic risk."
        )
        action = (
            "Immediately quarantine the herd: do not move animals off the premises. "
            "Notify the national veterinary authority (USDA APHIS Veterinary Services "
            "or national equivalent) and the state/provincial veterinarian today. "
            "Submit blood samples from the presumptive positive animal(s) for "
            "confirmatory serology (RBPT + CFT or i-ELISA) and culture at an "
            "accredited state or national reference laboratory. "
            "Advise all farm workers and veterinarians of potential zoonotic exposure; "
            "recommend medical evaluation for anyone with mucous-membrane or wound "
            "contact with birth materials, aborted fetuses, or vaginal discharge. "
            "Do not resample or breed the herd until confirmed negative."
        )
        return interp, action

    if ctx is BrucellaContext.SMALL_RUMINANT:
        interp = (
            "IS711 LAMP positive: Brucella genus detected in a small-ruminant "
            "(sheep or goat) sample. "
            "B. melitensis is the most likely species and is the most important "
            "cause of human brucellosis globally -- more infectious and virulent "
            "in humans than B. abortus. "
            "Ovine/caprine brucellosis (B. melitensis) is WOAH-listed and "
            "notifiable under all major national animal-disease programmes."
        )
        action = (
            "Immediately quarantine the flock: cease all animal movement. "
            "Notify the national veterinary authority today. "
            "B. melitensis carries a higher zoonotic risk to farm workers, "
            "shepherds, and veterinarians than B. abortus; all personnel with "
            "recent contact with aborted material, placenta, or birth fluids must "
            "be referred immediately for occupational health assessment and "
            "preventive treatment as indicated. "
            "Confirmatory serology (RBPT + CFT or i-ELISA) and culture at an "
            "accredited reference laboratory. "
            "Do not distribute milk or dairy products from the flock until cleared."
        )
        return interp, action

    if ctx is BrucellaContext.SWINE:
        interp = (
            "IS711 LAMP positive: Brucella genus detected in a porcine sample. "
            "B. suis is the most likely species; biovars 1 and 3 are zoonotic "
            "and are notifiable in most jurisdictions. "
            "Porcine brucellosis is WOAH-listed and clinically presents as "
            "reproductive failure, orchitis, and lameness in swine. "
            "Feral-pig contact is a known route of introduction to domestic herds."
        )
        action = (
            "Immediately quarantine the herd. "
            "Notify the national veterinary authority today (USDA APHIS VS or equivalent). "
            "Depopulation and disinfection of the premises are required in most "
            "jurisdictions for confirmed B. suis farms. "
            "Submit blood samples for confirmatory serology (RBPT + CFT) and "
            "culture at an accredited reference laboratory. "
            "Advise farm workers of zoonotic exposure risk; "
            "recommend medical evaluation for personnel with direct contact with "
            "aborted material or reproductive tissues. "
            "Investigate feral-pig contact as a potential source if the herd has "
            "outdoor access."
        )
        return interp, action

    # ENVIRONMENTAL
    interp = (
        "IS711 LAMP positive: Brucella genus DNA detected in an environmental sample. "
        "IS711 is a multi-copy insertion sequence with very high analytical sensitivity; "
        "environmental DNA from Brucella can persist in soil, water, or fomites after "
        "animal movement, parturition events, or aborted fetuses. "
        "A positive environmental result does not confirm active infection in the "
        "current animal population on the premises."
    )
    action = (
        "Review the animal history on the premises: recent livestock movements, "
        "parturitions, abortions, or feral-animal contact. "
        "Test associated animals (blood serology, milk ring test for dairy cattle, "
        "lymph node aspirate) using the IS711 LAMP assay or confirmatory serology "
        "to determine whether active infection is present. "
        "Do not report to the national veterinary authority on the basis of an "
        "environmental result alone -- regulatory reporting requires a positive "
        "result in an animal or human sample, confirmed by a reference laboratory. "
        "Repeat environmental sampling after animal testing to monitor contamination."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_brucella_risk(flags: BrucellaFlags) -> BrucellaAssessment:
    """Return a structured risk assessment for a Brucella IS711 LAMP result.

    Maps the IS711 positivity and sample context to an alert level, numeric
    score, contextualised interpretation, and field-actionable recommendation.

    Alert level logic:

    - Negative result: NEGATIVE (score 0) regardless of context.
    - Positive + HUMAN context: CRITICAL (score 90); zoonotic notifiable disease.
    - Positive + livestock context (CATTLE / SMALL_RUMINANT / SWINE):
      HIGH (score 60); WOAH-listed notifiable animal disease.
    - Positive + ENVIRONMENTAL context: LOW (score 25); follow-up animal
      testing required before regulatory reporting.

    Args:
        flags: IS711 result and sample context.

    Returns:
        :class:`BrucellaAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.brucella_risk import BrucellaContext, BrucellaFlags
        from lamp_forge.brucella_risk import assess_brucella_risk

        flags = BrucellaFlags(is711_positive=True, context=BrucellaContext.CATTLE)
        assessment = assess_brucella_risk(flags)
        print(assessment.alert_level)            # BrucellaAlertLevel.HIGH
        print(assessment.immediate_report_required)  # True
    """
    if not flags.is711_positive:
        level = BrucellaAlertLevel.NEGATIVE
        notifiable = False
        zoonotic_risk = False
        immediate_report = False
    elif flags.context is BrucellaContext.HUMAN:
        level = BrucellaAlertLevel.CRITICAL
        notifiable = True
        zoonotic_risk = True
        immediate_report = True
    elif flags.context in (
        BrucellaContext.CATTLE,
        BrucellaContext.SMALL_RUMINANT,
        BrucellaContext.SWINE,
    ):
        level = BrucellaAlertLevel.HIGH
        notifiable = True
        zoonotic_risk = True
        immediate_report = True
    else:  # ENVIRONMENTAL
        level = BrucellaAlertLevel.LOW
        notifiable = False
        zoonotic_risk = False
        immediate_report = False

    interpretation, recommended_action = _build_text(flags)

    return BrucellaAssessment(
        flags=flags,
        alert_level=level,
        alert_score=_ALERT_SCORE[level],
        interpretation=interpretation,
        recommended_action=recommended_action,
        notifiable=notifiable,
        zoonotic_risk=zoonotic_risk,
        immediate_report_required=immediate_report,
    )


def flags_from_dict(data: dict[str, object]) -> BrucellaFlags:
    """Construct :class:`BrucellaFlags` from a plain dict (e.g. parsed JSON).

    Expected keys:

    - ``is711_positive`` (bool): IS711 LAMP result.
    - ``context`` (str): One of ``cattle``, ``small_ruminant``, ``swine``,
      ``human``, ``environmental``.  Defaults to ``cattle`` if absent or invalid.

    Unknown keys are ignored.

    Args:
        data: Dict with optional field values.

    Returns:
        :class:`BrucellaFlags` instance.

    Raises:
        ValueError: If ``context`` is present but not a valid :class:`BrucellaContext`.
    """
    is711 = bool(data.get("is711_positive", False))
    raw_ctx = data.get("context", BrucellaContext.CATTLE.value)
    try:
        ctx = BrucellaContext(str(raw_ctx))
    except ValueError as exc:
        valid = ", ".join(c.value for c in BrucellaContext)
        raise ValueError(f"Invalid context {raw_ctx!r}. Valid values: {valid}") from exc
    return BrucellaFlags(is711_positive=is711, context=ctx)


def write_assessment_json(assessment: BrucellaAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`BrucellaAssessment` to serialise.
        path: Output file path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: BrucellaAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`BrucellaAssessment` to serialise.
        path: Output file path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        for key, val in assessment.to_dict().items():
            writer.writerow([key, val])
