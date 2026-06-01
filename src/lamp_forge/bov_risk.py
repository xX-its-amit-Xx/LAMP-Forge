"""Bovine respiratory LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable bovine respiratory
LAMP panel as a structured alert, prioritising pathogens by their regulatory
status, case-fatality rate, and field-containment urgency.

Bovine respiratory disease (BRD) -- colloquially "shipping fever" -- is the
leading cause of morbidity and mortality in beef feedlots and a significant
cause of dairy production loss.  A single portable multiplex LAMP panel that
simultaneously interrogates the five major BRD agents lets veterinarians make
treatment decisions on arrival, before culture results are available.

Pathogens covered
-----------------

    Target   Gene       Pathogen                              Regulatory status
    BRSV     N gene     Bovine respiratory syncytial virus    Not WOAH-listed
    BCOV     N gene     Bovine coronavirus                    Not WOAH-listed
    BVDV     5'-UTR     Bovine viral diarrhea virus (1+2)     Not listed; BVDV2 CAN
    IBR      gB/UL27    Inf. bovine rhinotracheitis (BoHV-1)  Not WOAH-listed; some
                        mandatory control programmes
    MHAE     lktA       Mannheimia haemolytica (leukotoxin A) Not listed; notifiable
                        in some jurisdictions as secondary BRD agent

Notes:
-----
* BRSV, BCoV, and BVDV are **RNA viruses** -- detection requires one-step
  RT-LAMP (Bst 2.0 WarmStart + NEB RTx at 63-65 degC).
* IBR (BoHV-1) is a **DNA herpesvirus** -- standard DNA-LAMP (no RT step).
* *Mannheimia haemolytica* is the primary **bacterial** secondary pathogen in
  BRD.  Its leukotoxin (*lktA*) gene encodes the RTX toxin responsible for
  the alveolar exudate and the characteristic lesion.  This is a chromosomal
  DNA target (standard LAMP, no RT).

Scoring model
-------------
Base weights reflect case-fatality rate, production impact, treatment
complexity, and the probability of co-infection synergy in a BRD outbreak:

    BRSV  : 35  Most common primary trigger of BRD outbreaks; immunosuppresses
                the respiratory tract, predisposing to bacterial secondary
                infection; no approved vaccine in all regions
    BCOV   : 20  Causes calf scours and winter dysentery in addition to
                respiratory disease; ubiquitous; combined with BRSV = high risk
    BVDV   : 30  Marked immunosuppression (persistently infected animals shed
                continuously); BVDV2 causes type-II BVD with severe
                haemorrhagic syndrome; mandatory control programmes in EU/UK
    IBR    : 25  BoHV-1 latency allows reservoir persistence; many countries
                have mandatory notifiable control programmes (EU, Scandinavia,
                UK)
    MHAE   : 20  Primary cause of BRD mortality when bacterial co-infection
                occurs; lktA-positive strains are uniformly virulent

Alert level boundaries
----------------------
    Score >= 50  -> CRITICAL   (immediate veterinary response; quarantine)
    Score >= 30  -> HIGH       (prompt veterinary attention; movement restriction)
    Score >= 20  -> MODERATE   (veterinary review; enhanced monitoring)
    Score >=  1  -> LOW        (single low-weight detection; monitor)
    Score ==  0  -> NEGATIVE   (all negative; confirm sample adequacy)

IBR mandatory reporting
-----------------------
``ibr_mandatory_notification`` is True when IBR (BoHV-1) is positive.  Many
jurisdictions (EU Council Directive 64/432/EEC, UK, Scandinavia, Switzerland)
run mandatory BoHV-1 eradication or control programmes that require notifying
the national competent authority within 24-48 h.  Check local regulations
before reporting.

References:
----------
Snowder GD et al. (2006) Bovine respiratory disease in feedlot cattle: economic
    impact and preventive measures. J Anim Sci 84(4):860-870.
    doi:10.2527/2006.844860x
Grissett GP et al. (2015) Structured literature review of responses of cattle
    to viral and bacterial pathogens causing bovine respiratory disease complex.
    J Vet Intern Med 29(3):771-794. doi:10.1111/jvim.12597
Baker JC (1995) The clinical manifestations of bovine viral diarrhea infection.
    Vet Clin North Am Food Anim Pract 11(3):425-445.
    doi:10.1016/S0749-0720(15)30460-6
Ackermann M & Engels M (2006) Pro and contra IBR eradication. Vet Microbiol
    113(3-4):293-302. doi:10.1016/j.vetmic.2005.11.043
Dabo SM et al. (2008) Pasteurella multocida and bovine respiratory disease.
    Anim Health Res Rev 9(2):245-259. doi:10.1017/S1466252308001557
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


class BovAlertLevel(StrEnum):
    """Bovine respiratory LAMP panel alert severity.

    CRITICAL: Multiple or highly virulent pathogens detected; immediate
        veterinary response required and quarantine of affected animals.
    HIGH: Single high-impact pathogen or moderate co-infection; prompt
        veterinary attention and movement restriction.
    MODERATE: Single lower-weight detection; veterinary review.
    LOW: Single very-low-weight detection; continue monitoring.
    NEGATIVE: All targets negative.  Confirm a 16S sample-adequacy control
        is positive before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class BovPanelFlags:
    """Positivity flags for one bovine respiratory LAMP panel result.

    Each flag is True if the corresponding LAMP assay returned a positive
    amplification signal for that pathogen.

    Attributes:
        brsv: Bovine respiratory syncytial virus (N gene) positive.
        bcov: Bovine coronavirus (N gene) positive.
        bvdv: Bovine viral diarrhea virus (5-UTR, types 1 and 2) positive.
        ibr: Infectious bovine rhinotracheitis / BoHV-1 (gB / UL27) positive.
        mhae: Mannheimia haemolytica (lktA leukotoxin A) positive.
        mbovis: Mycoplasma bovis (uvrC DNA repair gene) positive.
            6th BRD channel; pan-resistant strains require susceptibility
            testing before antimicrobial therapy.  Default False for
            backward compatibility with 5-channel panel outputs.
    """

    brsv: bool
    bcov: bool
    bvdv: bool
    ibr: bool
    mhae: bool
    mbovis: bool = False


@dataclass(frozen=True, slots=True)
class BovRiskAssessment:
    """Structured alert from a bovine respiratory LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        ibr_mandatory_notification: True when IBR (BoHV-1) is positive;
            many EU/UK/Scandinavian jurisdictions mandate notification.
        bacterial_coinfection: True when MHAE is positive and at least one
            viral target is also positive -- the highest-mortality BRD pattern.
        mbovis_pan_resistant: True when M. bovis is positive.  Field strains
            are intrinsically resistant to beta-lactams, aminoglycosides, and
            rifampin; susceptibility testing is required before therapy.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation.
    """

    flags: BovPanelFlags
    alert_level: BovAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    ibr_mandatory_notification: bool
    bacterial_coinfection: bool
    mbovis_pan_resistant: bool
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
            "active_targets": list(self.active_targets),
            "ibr_mandatory_notification": self.ibr_mandatory_notification,
            "bacterial_coinfection": self.bacterial_coinfection,
            "mbovis_pan_resistant": self.mbovis_pan_resistant,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "panel_flags": {
                "brsv": self.flags.brsv,
                "bcov": self.flags.bcov,
                "bvdv": self.flags.bvdv,
                "ibr": self.flags.ibr,
                "mhae": self.flags.mhae,
                "mbovis": self.flags.mbovis,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_BRSV = 35
_WEIGHT_BCOV = 20
_WEIGHT_BVDV = 30
_WEIGHT_IBR = 25
_WEIGHT_MHAE = 20
_WEIGHT_MBOVIS = 20

# Targets where a positive triggers regulatory notification obligations
# (jurisdiction-dependent -- flag and advise user to check local rules)
_IBR_NOTIFICATION_REQUIRED = True  # IBR/BoHV-1 mandatory programmes in EU/UK

_TARGET_DESCRIPTION: dict[str, str] = {
    "BRSV": "Bovine respiratory syncytial virus (N gene) -- primary viral BRD trigger",
    "BCOV": "Bovine coronavirus (N gene) -- respiratory + enteric; ubiquitous",
    "BVDV": "Bovine viral diarrhea virus (5-UTR) -- immunosuppressive; BVDV2 haemorrhagic",
    "IBR": "Infectious bovine rhinotracheitis / BoHV-1 (gB) -- latent; mandatory control in EU/UK",
    "MHAE": "Mannheimia haemolytica (lktA) -- primary bacterial BRD pathogen; high mortality",
    "MBOVIS": (
        "Mycoplasma bovis (uvrC) -- chronic pneumonia, mastitis, arthritis; "
        "pan-resistant to beta-lactams; susceptibility testing required"
    ),
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> BovAlertLevel:
    if score >= 50:
        return BovAlertLevel.CRITICAL
    if score >= 30:
        return BovAlertLevel.HIGH
    if score >= 20:
        return BovAlertLevel.MODERATE
    if score >= 1:
        return BovAlertLevel.LOW
    return BovAlertLevel.NEGATIVE


def _build_text(
    flags: BovPanelFlags,
    active_targets: tuple[str, ...],
    alert_level: BovAlertLevel,
    bacterial_coinfection: bool,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.
        bacterial_coinfection: True when MHAE co-detected with a viral target.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    any_positive = any([flags.brsv, flags.bcov, flags.bvdv, flags.ibr, flags.mhae, flags.mbovis])

    if not any_positive:
        return (
            "All bovine respiratory targets negative. No BRD pathogens detected.",
            "Confirm a 16S sample-adequacy control is positive before reporting a "
            "clean result. Re-test if clinical signs (nasal discharge, cough, "
            "pyrexia, increased respiratory rate) are observed.",
        )

    target_str = ", ".join(active_targets)
    n_viral = sum([flags.brsv, flags.bcov, flags.bvdv, flags.ibr])
    n_total = len(active_targets)

    # Highest-severity pattern: bacterial co-infection with viral trigger(s)
    if bacterial_coinfection:
        viral_detected = [t for t in ("BRSV", "BCOV", "BVDV", "IBR") if t in active_targets]
        viral_str = ", ".join(viral_detected)
        mbovis_note = (
            " MBOVIS also detected: M. bovis is intrinsically resistant to beta-lactams "
            "and aminoglycosides -- do NOT use penicillin/amoxicillin-based protocols; "
            "request susceptibility testing before therapy."
            if flags.mbovis
            else ""
        )
        interp = (
            f"Viral-bacterial BRD co-infection detected: {viral_str} + Mannheimia haemolytica "
            f"(lktA). This is the highest-mortality pattern in bovine respiratory disease -- "
            f"M. haemolytica exploits viral immunosuppression to colonise the lower respiratory "
            f"tract, producing leukotoxin (LktA) that destroys alveolar macrophages and "
            f"neutrophils. Combined viral + M. haemolytica BRD carries case-fatality rates of "
            f"5-20% in untreated feedlot cattle.{mbovis_note}"
        )
        action = (
            "IMMEDIATE: Initiate antimicrobial therapy (florfenicol, tulathromycin, or enrofloxacin "
            "per local formulary and susceptibility data). Implement metaphylaxis for in-contact "
            "animals showing early signs. Isolate affected animals. Notify site veterinarian. "
            "If BVDV is positive, screen for persistently infected (PI) animals -- a PI animal "
            "shedding virus is the index-case driver and must be removed from the herd."
        )
        return interp, action

    # M. bovis alone or with MHAE (no viral co-detection)
    if flags.mbovis and n_viral == 0:
        mhae_note = (
            " Concurrent M. haemolytica (lktA) detection raises the probability of "
            "polymicrobial bacterial pneumonia; initiate broad-spectrum antimicrobials "
            "that cover MHAE (florfenicol, tulathromycin) while awaiting M. bovis "
            "susceptibility data."
            if flags.mhae
            else ""
        )
        interp = (
            "Mycoplasma bovis (uvrC) detected. M. bovis is the leading mycoplasmal cause "
            "of bovine respiratory disease, mastitis, and polyarthritis. It is intrinsically "
            "resistant to all beta-lactam antibiotics (penicillins, cephalosporins), "
            "aminoglycosides, and rifampin -- standard BRD antimicrobial protocols are "
            "INEFFECTIVE. Susceptibility testing (MIC panel: florfenicol, enrofloxacin, "
            "tetracyclines, macrolides) is required before selecting therapy. Chronic "
            "pneumonia and arthritis sequelae carry high production loss even in survivors."
            f"{mhae_note}"
        )
        action = (
            "DO NOT use beta-lactam or aminoglycoside antibiotics for M. bovis. "
            "Submit samples for mycoplasma culture and antimicrobial susceptibility testing. "
            "While awaiting results, florfenicol or an approved tetracycline may be used "
            "empirically per local veterinary guidance. Segregate affected animals to "
            "prevent horizontal transmission (M. bovis spreads via nasal secretions and "
            "milk). Review colostrum sourcing -- infected dams are a key reservoir. "
            "Notify site veterinarian."
        )
        return interp, action

    # BVDV alone or dominant
    if flags.bvdv and not flags.mhae:
        co = [t for t in ("BRSV", "BCOV", "IBR") if t in active_targets]
        co_str = f" Co-detections: {', '.join(co)}." if co else ""
        interp = (
            "BVDV (bovine viral diarrhea virus) detected via the 5-UTR assay "
            "(pan-genotypic: BVDV-1 and BVDV-2). BVDV causes profound T-lymphocyte "
            "immunosuppression, predisposing to BRD and secondary bacterial pneumonia. "
            "BVDV-2 can cause haemorrhagic syndrome with high mortality in naive herds."
            f"{co_str}"
        )
        action = (
            "Screen for persistently infected (PI) animals using individual antigen ELISA or "
            "ear-notch RT-PCR -- a single PI animal continuously shedding BVDV at high titre "
            "is the most likely source. Quarantine and remove PI animals. Notify veterinarian. "
            "Review vaccination programme. Monitor for secondary bacterial BRD in the following "
            "2-3 weeks."
        )
        return interp, action

    # IBR alone or dominant
    if flags.ibr and not flags.mhae:
        co = [t for t in ("BRSV", "BCOV", "BVDV") if t in active_targets]
        co_str = f" Co-detections: {', '.join(co)}." if co else ""
        interp = (
            "IBR (infectious bovine rhinotracheitis / BoHV-1) detected via the gB assay. "
            "BoHV-1 establishes latency in trigeminal ganglia; seropositive animals may "
            "reactivate and shed under transport or weaning stress without clinical signs. "
            "IBR is subject to mandatory eradication or control programmes in many EU "
            "member states, the United Kingdom, and Scandinavia."
            f"{co_str}"
        )
        action = (
            "Notify the site veterinarian. Check local regulatory requirements for IBR "
            "notification (mandatory in many EU/UK jurisdictions). Restrict animal movement "
            "and check vaccination status of the affected group. If animals are unvaccinated, "
            "consider emergency ring vaccination with marker vaccine (DIVA strategy). "
            "Report to the national competent authority if mandatory control programme applies."
        )
        return interp, action

    # BRSV alone or with BCoV (no MHAE, no BVDV, no IBR)
    if flags.brsv and not any([flags.bvdv, flags.ibr, flags.mhae]):
        bcov_note = " Co-detection of BCoV." if flags.bcov else ""
        interp = (
            "BRSV (bovine respiratory syncytial virus) detected via the N-gene assay. "
            "BRSV is the most common primary viral trigger of BRD outbreaks in beef feedlots "
            "and dairy calves.  It causes airway epithelial damage and suppresses mucociliary "
            "clearance, creating the window for secondary Mannheimia haemolytica or Pasteurella "
            "multocida pneumonia.  No confirmed M. haemolytica co-detection at this time."
            f"{bcov_note}"
        )
        action = (
            "Implement BRSV biosecurity: reduce animal density, improve ventilation, minimise "
            "stress-inducing procedures during the outbreak. Initiate supportive care "
            "(NSAIDs for pyrexia). Monitor for secondary bacterial pneumonia and re-test "
            "within 3-5 days if respiratory signs progress.  Review BRSV vaccination schedule "
            "(modified-live vaccine series for naives)."
        )
        return interp, action

    # BCoV alone
    if flags.bcov and not any([flags.brsv, flags.bvdv, flags.ibr, flags.mhae]):
        interp = (
            "BCoV (bovine coronavirus) detected via the N-gene assay. BCoV causes neonatal "
            "calf diarrhoea ('white scours'), winter dysentery in adults, and respiratory "
            "disease in feedlot calves.  Respiratory BCoV predisposes to secondary bacterial "
            "BRD in a similar manner to BRSV."
        )
        action = (
            "Isolate affected animals. Implement supportive care (oral electrolytes for "
            "diarrhoeic calves; NSAIDs for pyrexia). Improve colostrum management for "
            "neonates -- passive transfer failure is the primary risk factor. Monitor for "
            "bacterial secondary infection and re-test in 3-5 days if respiratory signs worsen."
        )
        return interp, action

    # MHAE alone (no viral co-detection)
    if flags.mhae and n_viral == 0:
        interp = (
            "Mannheimia haemolytica (lktA) detected without a concurrent viral trigger. "
            "Bacterial BRD without detectable virus may reflect: (1) the viral stage "
            "has passed and virus is no longer shedding at detectable levels, (2) stress-induced "
            "M. haemolytica overgrowth from the nasopharyngeal reservoir without a viral predisposing "
            "event, or (3) a missed viral cause (re-test with a fresh swab if clinical signs are acute)."
        )
        action = (
            "Initiate appropriate antimicrobial therapy per veterinary prescription and "
            "susceptibility patterns (florfenicol, tulathromycin, enrofloxacin, or "
            "ceftiofur as first-line options). Notify site veterinarian. Re-test for "
            "viral co-pathogens if clinical response is incomplete within 48 h."
        )
        return interp, action

    # Generic multi-target fallback
    interp = (
        f"Bovine respiratory LAMP panel: {n_total} target(s) detected -- {target_str}. "
        f"Alert level: {alert_level.value}."
    )
    action = (
        "Consult the site veterinarian. Follow treatment protocols appropriate for each "
        "detected pathogen. Prioritise antimicrobial therapy if M. haemolytica is "
        "co-detected. Monitor herd for spread and re-test index animals in 3-5 days."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_bov_risk(flags: BovPanelFlags) -> BovRiskAssessment:
    """Compute a structured bovine respiratory alert from five-pathogen panel flags.

    Weights each pathogen by its BRD impact (viral immunosuppression capacity,
    bacterial virulence), sums to a 0-100 score, and maps it to one of five
    categorical alert levels.  Sets ``bacterial_coinfection`` when *M. haemolytica*
    (lktA) is positive alongside any viral target -- the highest-mortality pattern.

    Args:
        flags: Positivity flags for all five bovine respiratory panel targets.

    Returns:
        :class:`BovRiskAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.bov_risk import BovPanelFlags, assess_bov_risk

        flags = BovPanelFlags(brsv=True, bcov=False, bvdv=False, ibr=False, mhae=True)
        assessment = assess_bov_risk(flags)
        print(assessment.alert_level)          # BovAlertLevel.CRITICAL
        print(assessment.bacterial_coinfection)  # True
    """
    score = 0
    if flags.brsv:
        score += _WEIGHT_BRSV
    if flags.bcov:
        score += _WEIGHT_BCOV
    if flags.bvdv:
        score += _WEIGHT_BVDV
    if flags.ibr:
        score += _WEIGHT_IBR
    if flags.mhae:
        score += _WEIGHT_MHAE
    if flags.mbovis:
        score += _WEIGHT_MBOVIS

    score = max(0, min(100, score))

    target_order = [
        ("BRSV", flags.brsv),
        ("BCOV", flags.bcov),
        ("BVDV", flags.bvdv),
        ("IBR", flags.ibr),
        ("MHAE", flags.mhae),
        ("MBOVIS", flags.mbovis),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    any_viral = any([flags.brsv, flags.bcov, flags.bvdv, flags.ibr])
    bacterial_coinfection = flags.mhae and any_viral
    ibr_mandatory_notification = flags.ibr
    mbovis_pan_resistant = flags.mbovis

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(
        flags, active_targets, alert_level, bacterial_coinfection
    )

    return BovRiskAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        ibr_mandatory_notification=ibr_mandatory_notification,
        bacterial_coinfection=bacterial_coinfection,
        mbovis_pan_resistant=mbovis_pan_resistant,
        interpretation=interpretation,
        recommended_action=recommended_action,
    )


def flags_from_dict(data: dict[str, object]) -> BovPanelFlags:
    """Construct :class:`BovPanelFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``brsv``, ``bcov``, ``bvdv``, ``ibr``, ``mhae`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`BovPanelFlags` instance.
    """
    return BovPanelFlags(
        brsv=bool(data.get("brsv", False)),
        bcov=bool(data.get("bcov", False)),
        bvdv=bool(data.get("bvdv", False)),
        ibr=bool(data.get("ibr", False)),
        mhae=bool(data.get("mhae", False)),
        mbovis=bool(data.get("mbovis", False)),
    )


def write_assessment_json(assessment: BovRiskAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`BovRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: BovRiskAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`BovRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        d = assessment.to_dict()
        for key, val in d.items():
            if key == "active_targets" and isinstance(val, list):
                writer.writerow([key, "; ".join(str(v) for v in val)])
            elif key == "panel_flags" and isinstance(val, dict):
                for pk, pv in val.items():
                    writer.writerow([f"target_{pk}", "positive" if pv else "negative"])
            else:
                writer.writerow([key, val])
