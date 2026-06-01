"""CSFV NS5B RT-LAMP risk assessment engine.

Interprets a Classical Swine Fever Virus (CSFV) NS5B RT-LAMP result as a
structured alert, contextualised by the deployment setting.

CSFV (hog cholera; European swine fever) is a WOAH-listed Tier 1 notifiable
disease caused by Pestivirus A (*Flaviviridae*).  It produces near-100% case
fatality in naive pig herds and is clinically indistinguishable from African
Swine Fever Virus (ASFV) on presentation alone -- both cause acute hemorrhagic
fever, high mortality, and rapid spread.  A BioVind portable platform running
CSFV (NS5B) and ASFV (B646L/p72) assays in the same cartridge resolves this
differential at the farm gate or border checkpoint in under 60 minutes.

NS5B (RNA-directed RNA polymerase, ~1.8 kb) is the most conserved ORF across
all three CSFV genotype groups (1/2/3) and is the basis of published CSFV
RT-PCR and RT-LAMP assays (Wang et al. 2018; Rodriguez-Prieto et al. 2017).
CSFV is a positive-sense ssRNA virus; the NS5B assay requires one-step RT-LAMP
at 63-65 degC with NEB RTx or equivalent reverse transcriptase.

Unlike Brucella, CSFV is strictly a pig pathogen (Sus scrofa and wild boar
*Sus scrofa ferus*) with no established human infection under natural conditions.
The alert level is CRITICAL for any positive NS5B result regardless of context;
the deployment context (on-farm, border post, or livestock market) governs the
recommended containment action.

SPECIFICITY NOTE: BVDV (Bovine Viral Diarrhea Virus, Pestivirus B) is the
closest relative in the off-target panel, sharing ~50-60% NS5B nucleotide
identity with CSFV.  NS5B assays with `min_identity_threshold >= 0.85` and
BVDV in the off-target panel reliably distinguish the two pestiviruses.

BioVind vertical applicability
-------------------------------
Farm (animal biosecurity):
    CSFV is the primary differential diagnosis for ASFV.  Both cause
    hemorrhagic fever in swine with overlapping clinical signs; molecular
    discrimination is required before any depopulation or trade decision.
    CSFV is endemic in parts of Asia, Central America, and sub-Saharan
    Africa; eradicated from North America, the EU, and Oceania.

Field deployment context:
    On-farm: standard herd surveillance, post-morbidity investigation.
    Border post: live-pig or pork-product consignment inspection.
    Livestock market: aggregation-point screening before sale or movement.

References:
    Wang A, Cai X, Chen H et al. (2018) Rapid and sensitive detection of
        classical swine fever virus by reverse transcription loop-mediated
        isothermal amplification.
        Vet Microbiol 216:66-71. doi:10.1016/j.vetmic.2018.02.003
    Rodriguez-Prieto V, Vicente-Rubiano M, Sanchez-Matamoros A et al. (2017)
        Real-time and conventional RT-LAMP for detection of classical swine
        fever virus. J Virol Methods 250:7-12.
        doi:10.1016/j.jviromet.2017.09.014
    WOAH (2022). Classical swine fever. Terrestrial Manual, Chapter 3.8.3.
    Moennig V et al. (2003) Clinical signs and epidemiology of classical
        swine fever: a review of new knowledge. Vet J 165(1):11-20.
        doi:10.1016/S1090-0233(02)00112-0
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


class CsfvContext(StrEnum):
    """Deployment context for CSFV NS5B RT-LAMP result interpretation.

    Attributes:
        ON_FARM: Herd surveillance or post-morbidity investigation on a pig
            farm.  The affected herd must be quarantined pending confirmation.
        BORDER_POST: Consignment screening at a border inspection post or
            port of entry.  A positive blocks the consignment.
        MARKET: Screening at a livestock market or auction site before
            sale or inter-farm movement.  A positive halts all transactions
            and triggers enhanced contact tracing for recent buyers.
    """

    ON_FARM = "on_farm"
    BORDER_POST = "border_post"
    MARKET = "market"


class CsfvAlertLevel(StrEnum):
    """Categorical alert level for a CSFV NS5B RT-LAMP result.

    Attributes:
        CRITICAL: NS5B positive in any context.  CSFV is a WOAH-listed
            Tier 1 notifiable disease with near-100% case fatality in
            naive herds.  Immediate notification of the national veterinary
            authority is mandatory before any movement or trade decision.
        NEGATIVE: NS5B negative.  CSFV not detected in this sample.
    """

    CRITICAL = "CRITICAL"
    NEGATIVE = "NEGATIVE"


_ALERT_SCORE: dict[CsfvAlertLevel, int] = {
    CsfvAlertLevel.CRITICAL: 100,
    CsfvAlertLevel.NEGATIVE: 0,
}


@dataclass(frozen=True, slots=True)
class CsfvFlags:
    """NS5B RT-LAMP result and deployment context for one CSFV assessment.

    Attributes:
        ns5b_positive: True if the NS5B RT-LAMP assay returned a positive
            amplification signal.
        context: Deployment context governing the recommended action.
            Defaults to ON_FARM if not specified.
    """

    ns5b_positive: bool
    context: CsfvContext = CsfvContext.ON_FARM


@dataclass(frozen=True, slots=True)
class CsfvAssessment:
    """Structured alert from a CSFV NS5B RT-LAMP result.

    Attributes:
        flags: NS5B result and deployment context used to generate this
            assessment.
        alert_level: Categorical alert level (CRITICAL / NEGATIVE).
        alert_score: Numeric severity, 0-100 (higher is more urgent).
            100 for any positive result; 0 for negative.
        interpretation: Explanation of the NS5B result and its significance
            in the deployment context.
        recommended_action: Field-actionable containment and reporting
            recommendation tailored to the deployment context.
        immediate_report_required: True when a positive result requires
            same-day notification of the national veterinary authority.
            Always True for CRITICAL results; False for NEGATIVE.
        woah_notifiable: True when the result triggers a WOAH notification
            obligation.  CSFV is a WOAH-listed disease in all jurisdictions.
    """

    flags: CsfvFlags
    alert_level: CsfvAlertLevel
    alert_score: int
    interpretation: str
    recommended_action: str
    immediate_report_required: bool
    woah_notifiable: bool

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all assessment fields; flags are inlined under
            ``ns5b_result`` and ``context``.
        """
        return {
            "alert_level": self.alert_level.value,
            "alert_score": self.alert_score,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "immediate_report_required": self.immediate_report_required,
            "woah_notifiable": self.woah_notifiable,
            "ns5b_result": "positive" if self.flags.ns5b_positive else "negative",
            "context": self.flags.context.value,
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _build_text(flags: CsfvFlags) -> tuple[str, str]:
    """Return (interpretation, recommended_action) for the NS5B result.

    Args:
        flags: NS5B result and deployment context.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    if not flags.ns5b_positive:
        return (
            "NS5B RT-LAMP negative: Classical Swine Fever Virus (CSFV) not detected "
            "in this sample. A negative result in an appropriate sample type (blood, "
            "serum, tonsil swab) taken from acutely febrile animals is strong evidence "
            "against CSFV as the cause of the clinical signs. "
            "If ASFV (African Swine Fever Virus) is also negative on a co-run assay, "
            "consider other causes of hemorrhagic fever in swine (PCV2, Erysipelas, "
            "PRRSV-associated vasculitis).",
            "No immediate quarantine action based on CSFV LAMP alone. "
            "Confirm that an internal extraction control (e.g. swine actin or universal "
            "16S) is positive before reporting a clean result. "
            "If clinical signs persist or worsen, collect fresh samples from additional "
            "animals and repeat testing. Consult the national veterinary authority if "
            "the hemorrhagic disease cluster is not explained by the LAMP panel results.",
        )

    ctx = flags.context

    if ctx is CsfvContext.ON_FARM:
        interp = (
            "NS5B RT-LAMP positive: Classical Swine Fever Virus (CSFV) detected in "
            "a pig-farm sample. CSFV (hog cholera) is a WOAH-listed Tier 1 notifiable "
            "disease causing near-100% case fatality in naive herds within days of "
            "onset. All affected and in-contact pigs are presumptively exposed. "
            "This result is the critical differential for African Swine Fever Virus "
            "(ASFV): if an ASFV assay was run simultaneously and is negative, CSFV is "
            "the presumptive diagnosis; if ASFV is also positive, a mixed infection or "
            "assay error must be investigated by the reference laboratory."
        )
        action = (
            "IMMEDIATELY quarantine the farm: prohibit all live-pig, pork-product, "
            "and vehicle movement on and off the premises. "
            "Notify the national veterinary authority (e.g. USDA APHIS VS, national "
            "competent authority) and the state or provincial veterinarian TODAY -- "
            "CSFV is a WOAH-listed disease and same-day reporting is mandatory. "
            "Do not euthanase or depopulate animals before receiving instructions from "
            "the national authority; confirmatory testing at an approved reference "
            "laboratory (OIE Reference Laboratory for CSF) is required to trigger the "
            "official response. "
            "Collect fresh EDTA blood and tonsil swabs from all sick and in-contact "
            "animals and submit for RT-PCR, virus isolation, and serology. "
            "Implement enhanced biosecurity: restrict personnel, disinfect all entry "
            "points, and document all recent pig movements and contacts."
        )
        return interp, action

    if ctx is CsfvContext.BORDER_POST:
        interp = (
            "NS5B RT-LAMP positive: Classical Swine Fever Virus (CSFV) detected in "
            "a consignment sample at a border inspection post. "
            "CSFV is a WOAH-listed Tier 1 notifiable disease; introduction into a "
            "CSFV-free country through an infected consignment would constitute an "
            "acute disease event requiring immediate national response. "
            "The consignment must be regarded as infected pending confirmation."
        )
        action = (
            "IMMEDIATELY detain the entire consignment: halt all animal or "
            "pork-product movement pending official instructions. "
            "Notify the national veterinary authority and border inspection authority "
            "TODAY; CSFV introduction at a border point is a notifiable emergency "
            "under WOAH guidelines. "
            "Isolate the consignment in a secure holding area; do not allow contact "
            "with domestic pigs or other susceptible species. "
            "Submit confirmatory samples (blood, swab, tissue) to the national or "
            "OIE Reference Laboratory for CSF; official disposal or repatriation of "
            "the consignment must follow national competent-authority instructions. "
            "Document all other consignments that shared transport or handling "
            "facilities with this shipment for trace-forward investigation."
        )
        return interp, action

    # MARKET context
    interp = (
        "NS5B RT-LAMP positive: Classical Swine Fever Virus (CSFV) detected in a "
        "sample from a livestock market or auction. "
        "CSFV spreads rapidly via nose-to-nose contact, contaminated vehicles, and "
        "shared equipment; a market detection represents a high-risk amplification "
        "event -- multiple farms may already be exposed through animals bought and "
        "sold at the same venue. "
        "CSFV is a WOAH-listed Tier 1 notifiable disease requiring immediate action."
    )
    action = (
        "IMMEDIATELY halt all pig transactions and movements at the market site. "
        "Quarantine all pigs currently on the premises -- no animals may leave. "
        "Notify the national veterinary authority TODAY; a market detection triggers "
        "a structured trace-back/trace-forward investigation under national emergency "
        "response protocols. "
        "Obtain a complete record of all pigs, buyers, sellers, and vehicles "
        "involved in the past 21 days (the CSFV incubation period); trace-forward "
        "contact of all recent buyers is essential to contain spread. "
        "Submit confirmatory samples from all suspect animals to the national or "
        "OIE Reference Laboratory for CSF before releasing any animal from the site. "
        "Disinfect all facilities, vehicles, and equipment before market activity resumes."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_csfv_risk(flags: CsfvFlags) -> CsfvAssessment:
    """Return a structured risk assessment for a CSFV NS5B RT-LAMP result.

    Maps the NS5B positivity and deployment context to an alert level, numeric
    score, contextualised interpretation, and field-actionable recommendation.

    Alert level logic:

    - Negative NS5B result: NEGATIVE (score 0) regardless of context.
    - Positive NS5B result: CRITICAL (score 100) in all contexts.
      The recommended action varies by deployment context (on-farm,
      border post, or livestock market).

    Args:
        flags: NS5B result and deployment context.

    Returns:
        :class:`CsfvAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.csfv_risk import CsfvContext, CsfvFlags
        from lamp_forge.csfv_risk import assess_csfv_risk

        flags = CsfvFlags(ns5b_positive=True, context=CsfvContext.ON_FARM)
        assessment = assess_csfv_risk(flags)
        print(assessment.alert_level)              # CsfvAlertLevel.CRITICAL
        print(assessment.immediate_report_required)  # True
    """
    if flags.ns5b_positive:
        level = CsfvAlertLevel.CRITICAL
        immediate_report = True
        woah_notifiable = True
    else:
        level = CsfvAlertLevel.NEGATIVE
        immediate_report = False
        woah_notifiable = False

    interpretation, recommended_action = _build_text(flags)

    return CsfvAssessment(
        flags=flags,
        alert_level=level,
        alert_score=_ALERT_SCORE[level],
        interpretation=interpretation,
        recommended_action=recommended_action,
        immediate_report_required=immediate_report,
        woah_notifiable=woah_notifiable,
    )


def flags_from_dict(data: dict[str, object]) -> CsfvFlags:
    """Construct :class:`CsfvFlags` from a plain dict (e.g. parsed JSON).

    Expected keys:

    - ``ns5b_positive`` (bool): NS5B RT-LAMP result.
    - ``context`` (str): One of ``on_farm``, ``border_post``, ``market``.
      Defaults to ``on_farm`` if absent.

    Unknown keys are ignored.

    Args:
        data: Dict with optional field values.

    Returns:
        :class:`CsfvFlags` instance.

    Raises:
        ValueError: If ``context`` is present but not a valid
            :class:`CsfvContext` value.
    """
    ns5b = bool(data.get("ns5b_positive", False))
    raw_ctx = data.get("context", CsfvContext.ON_FARM.value)
    try:
        ctx = CsfvContext(str(raw_ctx))
    except ValueError as exc:
        valid = ", ".join(c.value for c in CsfvContext)
        raise ValueError(f"Invalid context {raw_ctx!r}. Valid values: {valid}") from exc
    return CsfvFlags(ns5b_positive=ns5b, context=ctx)


def write_assessment_json(assessment: CsfvAssessment, path: Path) -> None:
    """Write the CSFV assessment to a JSON file.

    Args:
        assessment: :class:`CsfvAssessment` to serialise.
        path: Output file path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: CsfvAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the CSFV assessment.

    Args:
        assessment: :class:`CsfvAssessment` to serialise.
        path: Output file path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        for key, val in assessment.to_dict().items():
            writer.writerow([key, val])
