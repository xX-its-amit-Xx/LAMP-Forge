"""Human POC respiratory differential LAMP panel risk assessment engine.

Interprets the positivity flags from a three-target portable respiratory
LAMP panel (hRSV + Influenza A + SARS-CoV-2) as a structured clinical alert
that directs the clinician to the correct antiviral -- or confirms purely
supportive management -- before the patient leaves the room.

Clinical context
----------------
Lower respiratory tract illness in infants, children, elderly, and
immunocompromised adults presents with overlapping symptoms across three
viral aetiologies that now each have distinct treatment pathways:

    hRSV   N gene   Human respiratory syncytial virus    nirsevimab (paediatric)
    IAV    M gene   Influenza A (pan-subtype)             oseltamivir within 48 h
    SARS2  N gene   SARS-CoV-2                            nirmatrelvir/remdesivir

A portable RT-LAMP panel resolves the differential at the bedside in
under 60 minutes, directing therapy within treatment windows -- exactly
the use-case for BioVind's POC vertical (Recipe 34, cookbook).

Scoring model
-------------
    SARS2 : 25  (treatment window 5 days; high-risk groups at severe-disease risk)
    IAV   : 20  (oseltamivir most effective within 48 h; influenza isolation)
    hRSV  : 10  (supportive for adults; paediatric nirsevimab for infants)

Alert level boundaries
----------------------
    Score >= 30  ->  HIGH      (any co-infection -- dual or triple viral burden)
    Score >= 15  ->  MODERATE  (SARS-CoV-2 or Influenza A detected alone)
    Score >=  5  ->  LOW       (hRSV detected alone; supportive care)
    Score ==  0  ->  NEGATIVE  (all targets negative)

No CRITICAL level is defined: unlike the general POC panel (poc_risk),
none of these three targets are mandatorily notifiable in most jurisdictions.
Influenza A subtypes associated with zoonotic risk (H5N1, H7N9) may warrant
public health notification in some jurisdictions -- but the pan-IAV M-gene
assay cannot subtype, so that determination remains with the clinician and
local public health authority.

References:
    Shi T et al. (2017) RSV ALRI burden. Lancet 390(10098):946-958.
        doi:10.1016/S0140-6736(17)30938-8
    Hamid S et al. (2022) Nirsevimab for prevention of RSV in infants.
        N Engl J Med 386:837-846. doi:10.1056/NEJMoa2110275
    WHO (2022) Pharmacological treatment of influenza A guidelines.
    Hammond J et al. (2022) Nirmatrelvir for high-risk COVID-19.
        N Engl J Med 386:1397-1408. doi:10.1056/NEJMoa2118542
    Nie K et al. (2017) Rapid hRSV detection by RT-LAMP.
        Arch Virol 162(4):1077-1083. doi:10.1007/s00705-016-3182-3
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


class PocRespAlertLevel(StrEnum):
    """Respiratory differential panel clinical alert severity.

    Attributes:
        HIGH: Two or more respiratory viruses co-detected; dual or triple
            antiviral consideration, enhanced isolation precautions.
        MODERATE: SARS-CoV-2 or Influenza A detected alone; treatment-window
            antivirals should be considered for eligible patients.
        LOW: hRSV detected alone; supportive care for most patients;
            nirsevimab consideration for infants and high-risk children.
        NEGATIVE: All three targets negative.
    """

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class PocRespFlags:
    """Three-target positivity flags for one respiratory differential LAMP result.

    Each flag is True if the corresponding RT-LAMP assay returned a positive
    amplification signal for that pathogen.

    Attributes:
        hrsv: Human respiratory syncytial virus (N gene, pan-RSV-A/B) positive.
        iav: Influenza A virus (M gene, pan-IAV) positive.
        sars2: SARS-CoV-2 (N gene) positive.
    """

    hrsv: bool
    iav: bool
    sars2: bool


@dataclass(frozen=True, slots=True)
class PocRespAssessment:
    """Structured clinical alert from a respiratory differential LAMP panel.

    Attributes:
        flags: Pathogen positivity flags used to compute this assessment.
        alert_level: Categorical alert level (HIGH -> NEGATIVE).
        alert_score: Numeric 0-100 summary; higher = more urgent.
        detected_targets: Labels of positive targets in canonical order.
        antiviral_indicated: True when Influenza A or SARS-CoV-2 positive;
            treatment-window antivirals (oseltamivir, nirmatrelvir) apply.
        rsv_indicated: True when hRSV positive; nirsevimab guidance for
            infants and high-risk children; supportive care for adults.
        dual_infection: True when two or more targets positive.
        interpretation: Concise clinical interpretation of the panel result.
        recommended_action: Actionable clinical recommendation for this
            detection pattern.
    """

    flags: PocRespFlags
    alert_level: PocRespAlertLevel
    alert_score: int
    detected_targets: tuple[str, ...]
    antiviral_indicated: bool
    rsv_indicated: bool
    dual_infection: bool
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
            "antiviral_indicated": self.antiviral_indicated,
            "rsv_indicated": self.rsv_indicated,
            "dual_infection": self.dual_infection,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "panel_flags": {
                "hrsv": self.flags.hrsv,
                "iav": self.flags.iav,
                "sars2": self.flags.sars2,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_SARS2 = 25
_WEIGHT_IAV = 20
_WEIGHT_HRSV = 10

_HIGH_THRESHOLD = 30
_MODERATE_THRESHOLD = 15
_LOW_THRESHOLD = 5


def _alert_level_from_score(score: int) -> PocRespAlertLevel:
    if score >= _HIGH_THRESHOLD:
        return PocRespAlertLevel.HIGH
    if score >= _MODERATE_THRESHOLD:
        return PocRespAlertLevel.MODERATE
    if score >= _LOW_THRESHOLD:
        return PocRespAlertLevel.LOW
    return PocRespAlertLevel.NEGATIVE


# ---------------------------------------------------------------------------
# Internal text builder
# ---------------------------------------------------------------------------


def _build_text(
    flags: PocRespFlags,
    detected_targets: tuple[str, ...],
    alert_level: PocRespAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        detected_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Tuple of (interpretation, recommended_action) strings.
    """
    hrsv, iav, sars2 = flags.hrsv, flags.iav, flags.sars2
    n_pos = sum([hrsv, iav, sars2])

    if n_pos == 0:
        return (
            "All three respiratory differential targets negative (hRSV, IAV, SARS-CoV-2). "
            "Respiratory viral illness due to these three pathogens is unlikely.",
            "Confirm sample adequacy with a positive extraction control before reporting "
            "a clean result. If clinical suspicion persists, consider alternative "
            "pathogens (human metapneumovirus, parainfluenza virus, rhinovirus) or "
            "a bacterial cause. Standard supportive care applies.",
        )

    # Triple co-infection
    if n_pos == 3:
        return (
            "hRSV, Influenza A, and SARS-CoV-2 all detected simultaneously. "
            "Triple respiratory viral co-infection is rare but documented in "
            "heavily immunocompromised patients and represents a high-severity "
            "clinical scenario. Antiviral treatment is indicated for both IAV "
            "(oseltamivir within 48 h) and SARS-CoV-2 (nirmatrelvir/ritonavir or "
            "remdesivir within 5 days for eligible patients). The hRSV co-infection "
            "compounds respiratory burden; nirsevimab is approved only for infant "
            "prophylaxis, not adult treatment.",
            "Enhanced respiratory and contact isolation. Consult infectious diseases "
            "immediately. Initiate oseltamivir for IAV and assess eligibility for "
            "COVID-19 antivirals per local protocol. Monitor closely for rapid "
            "deterioration -- triple co-infection carries a high risk of progression "
            "to ARDS, especially in the elderly and immunocompromised. Paediatric "
            "patients: discuss nirsevimab eligibility with paediatric ID consultant.",
        )

    # Dual: SARS2 + IAV
    if sars2 and iav:
        return (
            "SARS-CoV-2 (N gene) and Influenza A (M gene) co-detected. "
            "'Flurona' dual infection is clinically significant: both antivirals "
            "may be warranted and respiratory isolation applies for both pathogens "
            "simultaneously. hRSV is negative.",
            "Respiratory isolation (droplet + contact precautions). Consider "
            "antivirals for both: oseltamivir for IAV (within 48 h of symptom onset) "
            "and nirmatrelvir/ritonavir or remdesivir for SARS-CoV-2 (within 5 days "
            "for high-risk patients). Report to local surveillance per applicable "
            "protocol. Review drug interactions before prescribing both antivirals "
            "concurrently.",
        )

    # Dual: SARS2 + hRSV
    if sars2 and hrsv:
        return (
            "SARS-CoV-2 (N gene) and hRSV (N gene) co-detected. "
            "Co-infection with these two respiratory viruses increases the risk "
            "of severe lower respiratory tract disease, particularly in infants, "
            "the elderly, and immunocompromised adults. IAV is negative.",
            "Respiratory and contact isolation for SARS-CoV-2. Assess eligibility "
            "for COVID-19 antivirals (nirmatrelvir/ritonavir or remdesivir within "
            "5 days for high-risk patients). Supportive care for hRSV (hydration, "
            "supplemental oxygen as needed). For infants who have not received "
            "nirsevimab prophylaxis, consult paediatric guidelines on eligibility. "
            "Monitor oxygen saturation closely -- dual viral pneumonitis can "
            "progress rapidly.",
        )

    # Dual: IAV + hRSV
    if iav and hrsv:
        return (
            "Influenza A (M gene) and hRSV (N gene) co-detected. "
            "IAV + RSV co-infection is associated with increased severity in "
            "infants and elderly patients. Antiviral treatment for IAV is indicated "
            "within 48 h of symptom onset. SARS-CoV-2 is negative.",
            "Respiratory isolation. Initiate oseltamivir for Influenza A within "
            "48 h of symptom onset -- delay beyond this window substantially reduces "
            "efficacy. Supportive care for hRSV component (antipyretics, "
            "supplemental oxygen, hydration). For infants: assess nirsevimab "
            "eligibility per paediatric guidelines. Monitor for secondary bacterial "
            "pneumonia -- IAV disrupts mucosal barriers and predisposes to "
            "Staphylococcus aureus and Streptococcus pneumoniae superinfection.",
        )

    # Single: SARS2 only
    if sars2:
        return (
            "SARS-CoV-2 (N gene RT-LAMP) detected. Influenza A and hRSV are "
            "negative. Respiratory isolation is required. Treatment-window antivirals "
            "should be considered for high-risk patients.",
            "Respiratory and contact isolation. Assess eligibility for outpatient "
            "antivirals per local COVID-19 treatment protocol (nirmatrelvir/ritonavir "
            "within 5 days of symptom onset; remdesivir for those unable to take "
            "oral medications). Advise patient on isolation duration per current "
            "public-health guidance. Report to local surveillance if variant-of-"
            "concern guidance applies.",
        )

    # Single: IAV only
    if iav:
        return (
            "Influenza A (M gene RT-LAMP, pan-IAV) detected. hRSV and SARS-CoV-2 "
            "are negative. The M-gene assay detects all IAV subtypes (H1N1, H3N2, "
            "and potentially H5N1/H7N9 -- subtyping requires a follow-up assay). "
            "Oseltamivir is most effective within 48 h of symptom onset.",
            "Respiratory isolation (droplet precautions). Prescribe oseltamivir for "
            "eligible patients within 48 h. Antiviral benefit persists up to 5 days "
            "for hospitalised or severely ill patients. In healthcare or care-home "
            "settings, consider post-exposure prophylaxis for high-risk contacts. "
            "If zoonotic IAV exposure is suspected (H5N1 poultry contact, H7N9 "
            "market exposure), notify the local public health authority promptly -- "
            "the M-gene assay cannot exclude H5N1.",
        )

    # Single: hRSV only
    return (
        "Human RSV (N gene RT-LAMP, pan-RSV-A/B) detected. Influenza A and "
        "SARS-CoV-2 are negative. hRSV is the leading cause of acute lower "
        "respiratory illness in infants and causes significant morbidity in "
        "older adults and immunocompromised patients.",
        "Supportive care: supplemental oxygen to maintain SpO2 >= 94%, careful "
        "hydration, and antipyretics as needed. Nebulised bronchodilators may "
        "provide symptomatic relief in wheezing children but do not alter the "
        "course of infection. For infants who have not previously received "
        "nirsevimab prophylaxis, assess eligibility per current ACIP/JCVI "
        "guidance. Ribavirin is not recommended outside severely immunocompromised "
        "patients in specialised settings. Contact precautions are appropriate "
        "for hospitalised patients.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_poc_resp_risk(flags: PocRespFlags) -> PocRespAssessment:
    """Compute a structured clinical alert for a respiratory differential LAMP result.

    Args:
        flags: Positivity flags from the three-target respiratory panel.

    Returns:
        A PocRespAssessment with alert level, clinical flags, and text guidance.
    """
    score = (
        (_WEIGHT_SARS2 if flags.sars2 else 0)
        + (_WEIGHT_IAV if flags.iav else 0)
        + (_WEIGHT_HRSV if flags.hrsv else 0)
    )

    detected: list[str] = []
    if flags.hrsv:
        detected.append("hRSV")
    if flags.iav:
        detected.append("IAV")
    if flags.sars2:
        detected.append("SARS2")
    detected_targets = tuple(detected)

    alert_level = _alert_level_from_score(score)
    antiviral_indicated = flags.iav or flags.sars2
    rsv_indicated = flags.hrsv
    dual_infection = sum([flags.hrsv, flags.iav, flags.sars2]) >= 2

    interpretation, recommended_action = _build_text(flags, detected_targets, alert_level)

    return PocRespAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=min(score, 100),
        detected_targets=detected_targets,
        antiviral_indicated=antiviral_indicated,
        rsv_indicated=rsv_indicated,
        dual_infection=dual_infection,
        interpretation=interpretation,
        recommended_action=recommended_action,
    )


def flags_from_dict(d: dict[str, object]) -> PocRespFlags:
    """Parse a flags dict (e.g. from JSON) into a PocRespFlags instance.

    Missing keys default to False so that partial JSON payloads (e.g. from a
    device that only ran a subset of assays) are handled gracefully.

    Args:
        d: Dict with optional bool-valued keys ``hrsv``, ``iav``, ``sars2``.

    Returns:
        PocRespFlags instance.
    """
    return PocRespFlags(
        hrsv=bool(d.get("hrsv", False)),
        iav=bool(d.get("iav", False)),
        sars2=bool(d.get("sars2", False)),
    )


def write_assessment_json(assessment: PocRespAssessment, path: Path) -> None:
    """Write a PocRespAssessment to a JSON file.

    Args:
        assessment: The assessment to serialise.
        path: Destination file path (parent directory must exist).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(assessment.to_dict(), indent=2), encoding="utf-8")


def write_assessment_csv(assessment: PocRespAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of a PocRespAssessment.

    Args:
        assessment: The assessment to serialise.
        path: Destination file path (parent directory must exist).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, object]] = [
        ("alert_level", assessment.alert_level.value),
        ("alert_score", assessment.alert_score),
        ("detected_targets", ";".join(assessment.detected_targets)),
        ("antiviral_indicated", assessment.antiviral_indicated),
        ("rsv_indicated", assessment.rsv_indicated),
        ("dual_infection", assessment.dual_infection),
        ("hrsv", assessment.flags.hrsv),
        ("iav", assessment.flags.iav),
        ("sars2", assessment.flags.sars2),
        ("interpretation", assessment.interpretation),
        ("recommended_action", assessment.recommended_action),
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["key", "value"])
        writer.writerows(rows)
