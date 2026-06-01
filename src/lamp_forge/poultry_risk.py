"""Poultry-biosecurity LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable poultry-
biosecurity LAMP panel as a structured alert, prioritising pathogens by their
regulatory status, flock case-fatality rate, and zoonotic spillover risk.

This module provides a dedicated four-target interpretation layer for
commercial broiler, layer, and breeder-flock surveillance.  It complements the
mixed-livestock farm-biosecurity panel (:mod:`lamp_forge.farm_risk`), which
covers systemic swine and cattle pathogens alongside AIV and NDV.

Pathogens covered
-----------------

    Target  Gene   Pathogen                                  Regulatory status
    AIV     M      Avian influenza A (pan-IAV)               WOAH-listed (HPAI Tier 1)
    NDV     M      Newcastle disease virus (pan-NDV)         WOAH-listed
    IBDV    VP2    Infectious bursal disease virus (Gumboro) WOAH-listed (some regions)
    IBV     N      Infectious bronchitis virus (pan-IBV)     Not globally listed

Assay notes
-----------
* The AIV M-gene assay is pan-influenza-A; it detects all subtypes including
  HPAI H5N1 clade 2.3.4.4b and LPAI strains equally.  A positive result must
  be treated as a presumptive HPAI event pending HA/NA subtype confirmation by
  a WOAH Reference Laboratory.  Do not delay quarantine while awaiting subtype.
* The NDV M-gene assay is pan-NDV.  Velogenic strains cause near-100% flock
  mortality; lentogenic strains are subclinical and frequently vaccine-origin.
  Pathotyping of the F-gene cleavage site is required at a reference
  laboratory to distinguish velogenic from lentogenic strains.
* IBDV VP2 detects both classic and very virulent (vvIBDV) strains.  Very
  virulent IBDV causes 60-70% mortality; classic strains cause immunosuppression
  that predisposes flocks to secondary pathogens (E. coli, Marek's, IBV).
* IBV N-gene assay is pan-IBV and does not resolve serotypes.  Serotype-specific
  discrimination (Massachusetts, 793B, QX, etc.) requires S1 gene sequencing.

Scoring model
-------------
Base weights reflect flock case-fatality rate (CFR), regulatory status, and
zoonotic risk:

    AIV  : 60  WOAH Tier 1 notifiable; HPAI flock CFR ~100%; zoonotic
               spillover risk (H5N1, H7N9); 2022-24 HPAI wave eliminated
               hundreds of millions of commercial birds globally
    NDV  : 35  WOAH-listed; velogenic strains ~100% flock mortality; Tier 1
               global trade restriction trigger; zoonotic conjunctivitis risk
               but no sustained human transmission
    IBDV : 20  major immunosuppressive disease; vvIBDV up to 60-70%
               direct flock mortality; classic strains low direct CFR but
               massive secondary-infection and vaccine-interference losses;
               estimated USD 2 billion/yr global economic impact
    IBV  : 10  most prevalent avian respiratory pathogen; respiratory signs
               and egg production losses; usually 0-10% direct mortality;
               manageable with ring vaccination

Alert level boundaries
----------------------
    Score >= 60  -> CRITICAL   (AIV detected; presumptive HPAI; immediate reporting)
    Score >= 35  -> HIGH       (NDV alone or multi-pathogen co-detection)
    Score >= 20  -> MODERATE   (IBDV alone or IBDV + IBV; vet call)
    Score >= 10  -> LOW        (IBV alone; respiratory management)
    Score ==  0  -> NEGATIVE   (all negative; confirm extraction control positive)

Immediate reporting
-------------------
``immediate_report_required`` is set True whenever AIV or NDV is positive.
Both are WOAH-listed pathogens; notification to the national veterinary
authority is required in virtually all jurisdictions before any depopulation
or flock movement.

``zoonotic_risk`` is set True when AIV is positive.  Highly pathogenic avian
influenza A viruses can infect humans in close contact with infected birds;
all personnel entering an affected house must wear N95 respirators, eye
protection, and gloves before the subtype is confirmed.

``depopulation_risk`` is set True when AIV or NDV is positive.  These
WOAH-listed diseases commonly require emergency depopulation under national
response plans; this flag signals that the farm owner and national authority
should be consulted before any birds are moved.

References:
    Swayne DE (2012) Impact of vaccines and vaccination on global control of
        avian influenza. Avian Dis 56(4 suppl):818-828.
        doi:10.1637/10138-091511-Review.1
    Alexander DJ (2000) Newcastle disease and other avian paramyxoviruses.
        Rev Sci Tech OIE 19(2):443-462. doi:10.20506/rst.19.2.1231
    van den Berg TP (2000) Acute infectious bursal disease in poultry:
        a review. Avian Pathol 29(3):175-194.
        doi:10.1080/030794500412385
    Cavanagh D (2007) Coronavirus avian infectious bronchitis virus.
        Vet Res 38(2):281-297. doi:10.1051/vetres:2006055
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


class PoultryAlertLevel(StrEnum):
    """Poultry biosecurity LAMP panel alert severity.

    CRITICAL: AIV detected; presumptive HPAI event; all personnel must use
        enhanced PPE; immediate notification of the national veterinary
        authority required before any birds are moved or culled.
    HIGH: NDV detected, or NDV co-detected with IBDV or IBV; WOAH-listed
        pathogen; flock quarantine and confirmatory testing required.
    MODERATE: IBDV detected without AIV or NDV; immunosuppression risk;
        veterinary assessment and vaccination programme review required.
    LOW: IBV alone; respiratory management, biosecurity audit, and
        vaccination review with flock health veterinarian.
    NEGATIVE: All targets negative.  Confirm an internal extraction control
        (e.g. 16S LAMP or Avian Beta-Actin control) is positive before
        reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class PoultryFlags:
    """Positivity flags for one poultry biosecurity LAMP panel result.

    Each flag is True if the corresponding LAMP or RT-LAMP assay returned a
    positive amplification signal for that pathogen.

    Attributes:
        aiv: Avian influenza A (M gene, pan-IAV) positive.
        ndv: Newcastle disease virus (M gene, pan-NDV) positive.
        ibdv: Infectious bursal disease virus / Gumboro (VP2 gene) positive.
        ibv: Infectious bronchitis virus (N gene, pan-IBV) positive.
    """

    aiv: bool
    ndv: bool
    ibdv: bool
    ibv: bool


@dataclass(frozen=True, slots=True)
class PoultryAssessment:
    """Structured alert from a poultry biosecurity LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation.
        immediate_report_required: True when AIV or NDV is positive.  Both
            are WOAH-listed; national veterinary authority notification is
            legally required before depopulation or movement in most
            jurisdictions.
        zoonotic_risk: True when AIV is positive.  All personnel entering
            the affected house must wear N95 respirators, eye protection,
            and gloves irrespective of pathotype confirmation.
        depopulation_risk: True when AIV or NDV is positive.  WOAH-listed
            diseases commonly trigger statutory emergency depopulation;
            consult the national authority before acting.
    """

    flags: PoultryFlags
    alert_level: PoultryAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    interpretation: str
    recommended_action: str
    immediate_report_required: bool
    zoonotic_risk: bool
    depopulation_risk: bool

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
            "immediate_report_required": self.immediate_report_required,
            "zoonotic_risk": self.zoonotic_risk,
            "depopulation_risk": self.depopulation_risk,
            "panel_flags": {
                "aiv": self.flags.aiv,
                "ndv": self.flags.ndv,
                "ibdv": self.flags.ibdv,
                "ibv": self.flags.ibv,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_AIV = 60
_WEIGHT_NDV = 35
_WEIGHT_IBDV = 20
_WEIGHT_IBV = 10


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> PoultryAlertLevel:
    if score >= 60:
        return PoultryAlertLevel.CRITICAL
    if score >= 35:
        return PoultryAlertLevel.HIGH
    if score >= 20:
        return PoultryAlertLevel.MODERATE
    if score >= 10:
        return PoultryAlertLevel.LOW
    return PoultryAlertLevel.NEGATIVE


def _build_text(
    flags: PoultryFlags,
    active_targets: tuple[str, ...],
    alert_level: PoultryAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    aiv, ndv, ibdv, ibv = flags.aiv, flags.ndv, flags.ibdv, flags.ibv

    if not any([aiv, ndv, ibdv, ibv]):
        return (
            "All poultry biosecurity targets negative. "
            "No influenza A, Newcastle disease, IBDV, or IBV detected.",
            "Confirm an internal extraction control (e.g. Avian Beta-Actin LAMP "
            "or 16S channel) is positive before reporting a clean result. "
            "Re-test if clinical respiratory signs, increased mortality, or "
            "drop in egg production is observed in the flock.",
        )

    # AIV positive -- dominant signal regardless of co-detection
    if aiv:
        co: list[str] = []
        if ndv:
            co.append("NDV")
        if ibdv:
            co.append("IBDV")
        if ibv:
            co.append("IBV")
        co_str = "; also detected: " + ", ".join(co) + "." if co else "."
        interp = (
            "Avian influenza A virus (AIV) detected via the M-gene pan-IAV assay"
            + co_str
            + " This result must be treated as a presumptive highly pathogenic avian "
            "influenza (HPAI) event until the HA and NA subtypes are confirmed by a "
            "WOAH Reference Laboratory. HPAI H5N1 clade 2.3.4.4b causes near-100% "
            "flock mortality within 24-48 h and carries a zoonotic spillover risk for "
            "personnel in direct contact with infected birds."
        )
        action = (
            "CRITICAL: Immediately restrict all movement of birds, personnel, "
            "equipment, and vehicles on the premises. Notify the national veterinary "
            "authority (notifiable disease obligation) and the farm's accredited "
            "veterinarian. All personnel entering the affected house must wear N95 "
            "respirators, safety goggles, gloves, and disposable coveralls before the "
            "subtype is confirmed. Collect cloacal and oropharyngeal swabs from five "
            "to ten acutely affected birds for confirmatory RT-PCR and virus isolation "
            "at an approved WOAH Reference Laboratory. Do NOT depopulate before "
            "receiving official veterinary authority instruction, as premature "
            "depopulation may compromise confirmatory testing and indemnity claims."
        )
        return interp, action

    # NDV positive (no AIV)
    if ndv:
        co2: list[str] = []
        if ibdv:
            co2.append("IBDV")
        if ibv:
            co2.append("IBV")
        co_str2 = " Co-detected with " + ", ".join(co2) + "." if co2 else ""
        interp = (
            "Newcastle disease virus (NDV) detected via the M-gene pan-NDV assay."
            + co_str2
            + " Velogenic NDV strains cause near-100% flock mortality with acute "
            "neurological signs and haemorrhagic enteritis; lentogenic strains are "
            "subclinical and may be vaccine-origin. Pathotyping of the F-gene "
            "cleavage site at a reference laboratory is required to distinguish "
            "velogenic (IVPI > 1.2) from mesogenic and lentogenic strains. "
            "NDV is a WOAH-listed notifiable disease."
        )
        action = (
            "HIGH: Immediately quarantine the affected house and restrict all bird "
            "and personnel movement. Notify the national veterinary authority (NDV "
            "is notifiable in most jurisdictions) and the accredited veterinarian. "
            "Collect oropharyngeal and cloacal swabs from five to ten acutely "
            "affected birds for confirmatory RT-PCR and F-gene cleavage-site "
            "sequencing. Review the flock vaccination history (La Sota, Clone 30, "
            "or killed-oil-emulsion vaccines). If IBDV co-detected, immunosuppression "
            "may have compromised vaccine-induced NDV immunity -- reassess the entire "
            "flock vaccination programme with the flock veterinarian."
        )
        return interp, action

    # IBDV positive (no AIV, no NDV)
    if ibdv:
        co_str3 = " Co-detected with IBV." if ibv else ""
        interp = (
            "Infectious bursal disease virus (IBDV / Gumboro) detected via the "
            "VP2 assay." + co_str3 + " IBDV causes severe immunosuppression of the "
            "bursa of Fabricius in 3-6 week old chickens. Very virulent IBDV (vvIBDV) "
            "strains cause direct mortality of 60-70%; classic strains cause subclinical "
            "immunosuppression that amplifies losses from secondary E. coli infections, "
            "Marek's disease, and respiratory viruses including IBV."
        )
        action = (
            "Notify the flock health veterinarian. Assess birds aged 2-8 weeks for "
            "bursal swelling, vent pecking, ruffled feathers, and watery faeces. "
            "Collect bursae from acutely affected birds for histopathology and virus "
            "isolation to distinguish classic from vvIBDV strains; vvIBDV requires "
            "notification to the national authority in some jurisdictions. Review "
            "maternal antibody levels and the flock vaccination programme (Lukert, "
            "intermediate, intermediate-plus, or vector vaccines). "
            "If IBV co-detected, secondary respiratory disease is likely; increase "
            "monitoring frequency and review ventilation management."
        )
        return interp, action

    # IBV positive only
    if ibv and not aiv and not ndv and not ibdv:
        interp = (
            "Infectious bronchitis virus (IBV) detected via the N-gene pan-IBV assay. "
            "IBV is the most prevalent avian respiratory coronavirus; it causes tracheal "
            "rales, nasal discharge, and reduced egg production (quality and quantity). "
            "The N-gene assay does not resolve serotype; serotype-specific "
            "discrimination (Massachusetts, 793B, QX, IS/1494/06, etc.) requires "
            "S1 gene sequencing for targeted vaccination selection."
        )
        action = (
            "Notify the flock health veterinarian. Assess the flock for respiratory "
            "signs, mortality rate, egg quality (watery albumen, misshapen eggs), "
            "and feed and water consumption. Collect tracheal swabs from five to ten "
            "acutely affected birds for IBV RT-PCR and S1 gene sequencing to "
            "identify the circulating serotype and match vaccine strains. Review the "
            "current vaccination programme (live attenuated and/or killed-oil-emulsion "
            "vaccines); serotype mismatch is the most common cause of IBV vaccine "
            "failure. Improve ventilation and reduce stocking density if pen humidity "
            "and ammonia levels are high."
        )
        return interp, action

    # Fallback (should not be reached with correct flag logic)
    targets_str = ", ".join(active_targets)
    return (
        f"Detected: {targets_str}. {alert_level.value} alert.",
        "Review with flock health veterinarian. Follow standard biosecurity protocol.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_poultry_risk(flags: PoultryFlags) -> PoultryAssessment:
    """Compute a structured biosecurity alert from a four-target poultry panel.

    Weights each pathogen by flock case-fatality rate, regulatory status, and
    zoonotic risk; sums to a 0-100 score; and maps it to one of five categorical
    alert levels.  Sets ``immediate_report_required`` when AIV or NDV is
    positive, ``zoonotic_risk`` when AIV is positive, and ``depopulation_risk``
    when AIV or NDV is positive.

    Args:
        flags: Positivity flags for all four poultry biosecurity targets.

    Returns:
        :class:`PoultryAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.poultry_risk import PoultryFlags, assess_poultry_risk

        flags = PoultryFlags(aiv=True, ndv=False, ibdv=False, ibv=False)
        assessment = assess_poultry_risk(flags)
        print(assessment.alert_level)              # PoultryAlertLevel.CRITICAL
        print(assessment.immediate_report_required)  # True
        print(assessment.zoonotic_risk)            # True
    """
    score = 0
    if flags.aiv:
        score += _WEIGHT_AIV
    if flags.ndv:
        score += _WEIGHT_NDV
    if flags.ibdv:
        score += _WEIGHT_IBDV
    if flags.ibv:
        score += _WEIGHT_IBV

    score = max(0, min(100, score))

    target_order = [
        ("AIV", flags.aiv),
        ("NDV", flags.ndv),
        ("IBDV", flags.ibdv),
        ("IBV", flags.ibv),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, active_targets, alert_level)

    return PoultryAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        interpretation=interpretation,
        recommended_action=recommended_action,
        immediate_report_required=flags.aiv or flags.ndv,
        zoonotic_risk=flags.aiv,
        depopulation_risk=flags.aiv or flags.ndv,
    )


def flags_from_dict(data: dict[str, object]) -> PoultryFlags:
    """Construct :class:`PoultryFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``aiv``, ``ndv``, ``ibdv``, ``ibv`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`PoultryFlags` instance.
    """
    return PoultryFlags(
        aiv=bool(data.get("aiv", False)),
        ndv=bool(data.get("ndv", False)),
        ibdv=bool(data.get("ibdv", False)),
        ibv=bool(data.get("ibv", False)),
    )


def write_assessment_json(assessment: PoultryAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`PoultryAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: PoultryAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`PoultryAssessment` to serialise.
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
