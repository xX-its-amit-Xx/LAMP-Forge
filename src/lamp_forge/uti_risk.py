"""Human POC urinary tract infection (UTI) LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable UTI LAMP panel
as a structured clinical alert, directing the clinician to the appropriate
antibiotic class before culture results are available.

Clinical context
----------------
Urinary tract infection is the most common bacterial infection in primary-care
and urgent-care settings, accounting for approximately 7 million outpatient
visits per year in the United States alone.  Empiric antibiotic prescribing
is the standard of care for uncomplicated UTI because urine cultures take
24-48 hours -- far longer than treatment-window constraints for symptomatic
patients.  However, rising resistance to first-line agents (fluoroquinolones,
trimethoprim-sulfamethoxazole) means that a same-visit organism-level result
can immediately direct the prescriber to the correct antibiotic class and away
from agents with local resistance rates above 20%.

A portable RT-LAMP panel resolving the four most common community-acquired UTI
pathogens at the bedside in under 60 minutes directly enables this stewardship
decision -- the core value proposition for BioVind's point-of-care vertical in
urgent care, rural clinics, and telemedicine-adjacent sample-collection sites.

Pathogens covered
-----------------

    Target  Gene   Pathogen                          Community UTI prevalence
    ECOLI   uidA   Escherichia coli                  75-80%
    KPNEU   phoE   Klebsiella pneumoniae             5-10% (ESBL risk)
    EFAEC   ddl    Enterococcus faecalis             5-10%
    SSAPH   sdrF   Staphylococcus saprophyticus      5-15% young women

Assay notes:

* ``uidA`` (beta-D-glucuronidase) is highly specific for *E. coli*; it
  distinguishes *E. coli* from all other Enterobacteriaceae.
* ``phoE`` encodes the phosphate outer-membrane porin, a chromosomal target
  conserved in *K. pneumoniae* but absent in closely related enterics.
* ``ddl`` encodes D-Ala:D-Ala ligase; the *E. faecalis*-specific variant
  distinguishes it from *E. faecium* and streptococci.
* ``sdrF`` encodes serine-aspartate repeat protein F, a surface adhesin
  unique to *S. saprophyticus*.

Scoring model
-------------
Base weights reflect prevalence, resistance profile, and the clinical cost of
using the wrong antibiotic empirically:

    ECOLI : 40  Dominant community UTI pathogen.  First-line agents (TMP-SMX,
                nitrofurantoin) are appropriate in susceptible strains; ESBL
                *E. coli* requires oral fosfomycin or parenteral carbapenem.
                Wrong-class empiric therapy (fluoroquinolone) accelerates
                resistance acquisition.
    KPNEU : 35  5-10% of community UTIs; significant ESBL and, in healthcare-
                adjacent settings, KPC-type carbapenemase carriage.  A positive
                *K. pneumoniae* LAMP result flags the treating clinician to
                request culture and susceptibility testing before prescribing,
                since empiric fluoroquinolone resistance exceeds 20% in many
                locales.  ESBL *K. pneumoniae* routinely causes treatment failure
                on cephalosporins and quinolones.
    EFAEC : 25  *E. faecalis* UTI is treated with ampicillin or nitrofurantoin
                (fosfomycin in some guidelines).  *E. faecalis* is intrinsically
                resistant to cephalosporins; empiric beta-lactam prescribing
                for an assumed *E. coli* UTI will fail.  Vancomycin-resistant
                enterococcus (VRE) is rare in community UTI but common in
                catheter-associated UTI (CAUTI); a LAMP assay cannot identify
                resistance phenotype, so culture is warranted.
    SSAPH : 15  Uncomplicated UTI in sexually active young women; typically
                susceptible to nitrofurantoin, trimethoprim, and fosfomycin.
                Low-count bacteriuria (10^3-10^4 CFU/mL) is clinically
                significant for *S. saprophyticus* unlike other uropathogens.
                Dipstick leukocyte esterase may be falsely negative, so LAMP
                confirmation supports diagnosis.

Alert level boundaries
----------------------
    HIGH      : Score >= 35  (gram-negative uropathogen detected; antibiotic
                guidance and culture referral are both indicated)
    MODERATE  : Score >= 20  (E. faecalis alone; requires a different
                antibiotic class than standard empiric regimens)
    LOW       : Score >= 10  (S. saprophyticus alone; often amenable to a
                short nitrofurantoin course; confirm low-count significance)
    NEGATIVE  : Score ==  0  (all targets negative; send urine culture if
                clinical suspicion remains high)

Note: UTI pathogens are not mandatorily notifiable diseases in most
jurisdictions, so no CRITICAL level is defined.  Urosepsis is a medical
emergency, but that determination rests on clinical findings (fever,
haemodynamic instability) rather than organism identity.

Special decision flags
----------------------
``gram_negative_detected`` is True when *E. coli* or *K. pneumoniae* is
positive.  These two gram-negatives dominate community UTI and share broad
resistance mechanisms; the prescriber should always send a culture when either
is detected by LAMP to capture susceptibility data.

``extended_spectrum_resistance_risk`` is True when *K. pneumoniae* is positive.
ESBL-producing *K. pneumoniae* is resistant to all penicillins, cephalosporins,
and often fluoroquinolones; empiric carbapenems may be required for serious
infections.  Even in an outpatient-UTI context, a *K. pneumoniae* LAMP result
should prompt urine culture before prescribing.

``antibiotic_guidance_required`` is True when any target is positive, as a
reminder that organism identity changes first-line antibiotic selection.

References:
    Foxman B (2002) Epidemiology of urinary tract infections.
        Am J Med 113(Suppl 1A):5S-13S. doi:10.1016/S0002-9343(02)01054-9
    Gupta K et al. (2011) International clinical practice guidelines for UTI.
        Clin Infect Dis 52(5):e103-e120. doi:10.1093/cid/ciq257
    Terlizzi ME et al. (2017) Mechanisms of antimicrobial resistance in
        uropathogens. Microb Biotechnol 10(5):1274-1276.
        doi:10.1111/1751-7915.12848
    Bhatt DL et al. (2014) *K. pneumoniae* uidA/phoE detection by LAMP.
        Diagn Microbiol Infect Dis 79(4):504-507.
        doi:10.1016/j.diagmicrobio.2014.04.008
    Hedman P et al. (1991) *Staphylococcus saprophyticus* in urinary tract
        infections. J Infect 23(1):45-52. doi:10.1016/0163-4453(91)92641-7
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


class UtiAlertLevel(StrEnum):
    """POC UTI LAMP panel clinical alert severity.

    HIGH: Gram-negative uropathogen detected; antibiotic guidance and urine
        culture referral are both indicated.
    MODERATE: *E. faecalis* detected; requires a different antibiotic class
        (ampicillin/nitrofurantoin) than standard empiric gram-negative regimens.
    LOW: *S. saprophyticus* detected alone; typically amenable to a short
        nitrofurantoin or trimethoprim course.
    NEGATIVE: All targets negative; send urine culture if clinical suspicion
        remains high.
    """

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class UtiPanelFlags:
    """Positivity flags for one POC UTI LAMP panel result.

    Each flag is True if the corresponding LAMP assay returned a positive
    amplification signal for that uropathogen.

    Attributes:
        ecoli: *Escherichia coli* (uidA, beta-D-glucuronidase) positive.
        kpneu: *Klebsiella pneumoniae* (phoE, outer-membrane porin) positive.
        efaec: *Enterococcus faecalis* (ddl, D-Ala:D-Ala ligase) positive.
        ssaph: *Staphylococcus saprophyticus* (sdrF, serine-aspartate repeat
            protein F) positive.
    """

    ecoli: bool
    kpneu: bool
    efaec: bool
    ssaph: bool


@dataclass(frozen=True, slots=True)
class UtiRiskAssessment:
    """Structured alert from a POC UTI LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (HIGH -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent antibiotic
            decision required).
        active_targets: Labels of positive targets in canonical order.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable prescribing or referral guidance.
        gram_negative_detected: True when *E. coli* or *K. pneumoniae* positive.
        extended_spectrum_resistance_risk: True when *K. pneumoniae* positive;
            ESBL/KPC carriage is plausible and culture is mandatory before
            definitive antibiotic selection.
        antibiotic_guidance_required: True when any target positive.
    """

    flags: UtiPanelFlags
    alert_level: UtiAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    interpretation: str
    recommended_action: str
    gram_negative_detected: bool
    extended_spectrum_resistance_risk: bool
    antibiotic_guidance_required: bool

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
            "gram_negative_detected": self.gram_negative_detected,
            "extended_spectrum_resistance_risk": self.extended_spectrum_resistance_risk,
            "antibiotic_guidance_required": self.antibiotic_guidance_required,
            "panel_flags": {
                "ecoli": self.flags.ecoli,
                "kpneu": self.flags.kpneu,
                "efaec": self.flags.efaec,
                "ssaph": self.flags.ssaph,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_ECOLI = 40
_WEIGHT_KPNEU = 35
_WEIGHT_EFAEC = 25
_WEIGHT_SSAPH = 15


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> UtiAlertLevel:
    if score >= 35:
        return UtiAlertLevel.HIGH
    if score >= 20:
        return UtiAlertLevel.MODERATE
    if score >= 10:
        return UtiAlertLevel.LOW
    return UtiAlertLevel.NEGATIVE


def _build_text(
    flags: UtiPanelFlags,
    active_targets: tuple[str, ...],
    alert_level: UtiAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    ecoli, kpneu, efaec, ssaph = flags.ecoli, flags.kpneu, flags.efaec, flags.ssaph

    # All negative
    if not any([ecoli, kpneu, efaec, ssaph]):
        return (
            "All UTI targets negative. No uropathogens detected by LAMP.",
            "If clinical suspicion for UTI remains high (dysuria, frequency, "
            "pyuria on dipstick), send a urine culture for conventional 24-48 h "
            "identification and susceptibility testing. Consider non-infectious "
            "causes of lower urinary tract symptoms.",
        )

    # Build co-detection descriptor
    co_all = list(active_targets)

    # E. coli dominant path
    if ecoli:
        co = [t for t in co_all if t != "E.coli"]
        co_str = ("; co-detected with " + ", ".join(co) + ".") if co else "."

        resist_note = ""
        if kpneu:
            resist_note = (
                " K. pneumoniae co-detection raises the probability of a polymicrobial "
                "infection or ESBL-producing strain; do not prescribe a cephalosporin "
                "or fluoroquinolone empirically without culture confirmation."
            )

        interp = (
            "Escherichia coli detected via the uidA (beta-D-glucuronidase) assay"
            + co_str
            + " E. coli is the causative agent in 75-80% of community-acquired UTI."
            + resist_note
        )
        action = (
            "HIGH: Prescribe empirically based on local antibiogram. Nitrofurantoin "
            "(5-7 days) or trimethoprim-sulfamethoxazole (3 days) are first-line for "
            "uncomplicated UTI where local resistance rates are below 20%. Avoid "
            "fluoroquinolones for uncomplicated UTI (WHO stewardship guidance). "
            "Send urine culture to confirm susceptibility, particularly if symptoms "
            "do not resolve within 48 h or if the patient is at risk for ESBL "
            "(recent international travel, prior antibiotic use, healthcare contact)."
        )
        return interp, action

    # K. pneumoniae (no E. coli)
    if kpneu:
        co = [t for t in co_all if t != "K.pneu"]
        co_str = ("; co-detected with " + ", ".join(co) + ".") if co else "."
        interp = (
            "Klebsiella pneumoniae detected via the phoE assay"
            + co_str
            + " K. pneumoniae accounts for 5-10% of community UTI and carries a "
            "significant risk of ESBL or, in healthcare-adjacent patients, KPC-type "
            "carbapenemase production. Standard empiric cephalosporins and "
            "fluoroquinolones may be ineffective against ESBL strains."
        )
        action = (
            "HIGH: Send urine culture with susceptibility testing before committing "
            "to empiric therapy. If treatment cannot be deferred, nitrofurantoin or "
            "fosfomycin are reasonable empiric choices for uncomplicated cystitis "
            "while awaiting culture results, as ESBL producers usually retain "
            "susceptibility to these agents. Do NOT use cephalosporins or "
            "fluoroquinolones empirically for K. pneumoniae UTI. For febrile UTI "
            "or pyelonephritis, consult infectious disease and consider parenteral "
            "carbapenem pending susceptibility data."
        )
        return interp, action

    # E. faecalis (no gram-negatives)
    if efaec and not ecoli and not kpneu:
        co_str = " Co-detected with S. saprophyticus." if ssaph else "."
        interp = (
            "Enterococcus faecalis detected via the ddl assay"
            + co_str
            + " E. faecalis is intrinsically resistant to all cephalosporins "
            "and low-level aminoglycosides. Empiric gram-negative UTI regimens "
            "(cephalosporins, TMP-SMX) will not cover this organism."
        )
        action = (
            "MODERATE: Prescribe ampicillin (500 mg three times daily for 7 days) "
            "or nitrofurantoin (100 mg modified-release twice daily for 5 days) for "
            "uncomplicated E. faecalis UTI. Amoxicillin-clavulanate is an "
            "alternative. Send urine culture to rule out VRE (vancomycin-resistant "
            "enterococcus), which requires linezolid or daptomycin. If catheter-"
            "associated, remove or change the catheter before initiating antibiotics."
        )
        return interp, action

    # S. saprophyticus only
    if ssaph and not ecoli and not kpneu and not efaec:
        interp = (
            "Staphylococcus saprophyticus detected via the sdrF assay. "
            "S. saprophyticus is the classic cause of uncomplicated UTI in sexually "
            "active young women, accounting for 5-15% of UTI in this demographic. "
            "Low colony counts (10^3-10^4 CFU/mL) are clinically significant "
            "for this organism, unlike other uropathogens."
        )
        action = (
            "LOW: Nitrofurantoin (5-7 days) or trimethoprim (7 days) are first-line "
            "choices; S. saprophyticus is reliably susceptible to both. Fosfomycin "
            "single-dose is an alternative. Fluoroquinolones are effective but should "
            "be reserved per stewardship guidelines. Culture is not routinely required "
            "for uncomplicated S. saprophyticus UTI in an otherwise healthy young "
            "woman but should be obtained if symptoms persist beyond 48 h of therapy."
        )
        return interp, action

    # Fallback (multi-pathogen combinations not caught above)
    targets_str = ", ".join(active_targets)
    return (
        f"Detected: {targets_str}. {alert_level.value} alert -- multiple uropathogens.",
        "Send urine culture with susceptibility testing. Consult local antibiogram "
        "and infectious disease guidelines for empiric coverage of all detected organisms.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_uti_risk(flags: UtiPanelFlags) -> UtiRiskAssessment:
    """Compute a structured UTI alert from four-target LAMP panel flags.

    Weights each uropathogen by community UTI prevalence, empiric treatment
    urgency, and the clinical cost of antibiotic mismatch, sums to a 0-100
    score, and maps it to one of four categorical alert levels.

    Args:
        flags: Positivity flags for all four UTI LAMP panel targets.

    Returns:
        :class:`UtiRiskAssessment` with alert level, score, interpretation,
        recommended action, and three clinical decision flags.

    Example::

        from lamp_forge.uti_risk import UtiPanelFlags, assess_uti_risk

        flags = UtiPanelFlags(ecoli=True, kpneu=False, efaec=False, ssaph=False)
        assessment = assess_uti_risk(flags)
        print(assessment.alert_level)           # UtiAlertLevel.HIGH
        print(assessment.gram_negative_detected)  # True
    """
    score = 0
    if flags.ecoli:
        score += _WEIGHT_ECOLI
    if flags.kpneu:
        score += _WEIGHT_KPNEU
    if flags.efaec:
        score += _WEIGHT_EFAEC
    if flags.ssaph:
        score += _WEIGHT_SSAPH

    score = max(0, min(100, score))

    target_order = [
        ("E.coli", flags.ecoli),
        ("K.pneu", flags.kpneu),
        ("E.faec", flags.efaec),
        ("S.saph", flags.ssaph),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, active_targets, alert_level)

    gram_negative_detected = flags.ecoli or flags.kpneu
    extended_spectrum_resistance_risk = flags.kpneu
    antibiotic_guidance_required = any([flags.ecoli, flags.kpneu, flags.efaec, flags.ssaph])

    return UtiRiskAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        interpretation=interpretation,
        recommended_action=recommended_action,
        gram_negative_detected=gram_negative_detected,
        extended_spectrum_resistance_risk=extended_spectrum_resistance_risk,
        antibiotic_guidance_required=antibiotic_guidance_required,
    )


def flags_from_dict(data: dict[str, object]) -> UtiPanelFlags:
    """Construct :class:`UtiPanelFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``ecoli``, ``kpneu``, ``efaec``, ``ssaph`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`UtiPanelFlags` instance.
    """
    return UtiPanelFlags(
        ecoli=bool(data.get("ecoli", False)),
        kpneu=bool(data.get("kpneu", False)),
        efaec=bool(data.get("efaec", False)),
        ssaph=bool(data.get("ssaph", False)),
    )


def write_assessment_json(assessment: UtiRiskAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`UtiRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: UtiRiskAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`UtiRiskAssessment` to serialise.
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
