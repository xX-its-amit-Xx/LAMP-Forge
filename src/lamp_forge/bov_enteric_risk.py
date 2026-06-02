"""Bovine neonatal calf diarrhea LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable bovine neonatal
enteric LAMP panel as a structured alert, prioritising pathogens by their
time-to-death in untreated neonatal calves, systemic-disease risk, and the
urgency of antimicrobial intervention.

Neonatal calf diarrhea (NCD / "scours") is the leading cause of calf mortality
worldwide, responsible for more deaths in beef and dairy calves than all other
diseases combined.  Unlike the adult bovine respiratory panel (:mod:`lamp_forge
.bov_risk`), which covers BRSV/BCoV/BVDV/IBR/M. haemolytica/M. bovis, the
enteric panel covers the five pathogens implicated in neonatal scours.  The
critical clinical question is not just "what is present?" but **"does this
calf need intravenous antibiotics within the next two hours?"**  A portable
LAMP result distinguishes bacterial (ETEC, Salmonella) from viral (BCoV,
RotaA) and parasitic (Crypto) disease at the pen-side, allowing the
farm veterinarian or calf-care technician to make that decision without
waiting for a culture result.

Pathogens covered
-----------------

    Target      Gene     Pathogen                            CFR untreated
    ETEC_K99    faeG     ETEC K99/F5 (enterotoxigenic)      >90% neonates <4d
    SALMONELLA  invA     Salmonella sp. (Dublin/Typhimurium) 40-70% neonates
    CRYPTO      18S/hsp  Cryptosporidium parvum              10-50% neonates
    BCOV        N gene   Bovine coronavirus (enteric)        5-20% neonates
    ROTA_A      VP6      Bovine rotavirus A                  5-15% neonates

Scoring model
-------------
Weights reflect time-to-death without intervention and antimicrobial urgency:

    ETEC_K99   : 50  Near-certain mortality in calves <4 days old without IV
                     antibiotics; kills within 24-48 h via hypersecretory
                     diarrhoea and metabolic acidosis.  Prompt IV ampicillin
                     or amoxicillin-clavulanate is the standard of care.
    SALMONELLA : 35  Bacteremia and septicemia risk in calves <2 weeks old;
                     40-70% mortality without treatment; zoonotic (S. Dublin
                     is a human food-safety pathogen).  Antibiotic therapy
                     required.
    CRYPTO     : 20  10-50% neonatal mortality; unique treatment (halofuginone
                     prophylaxis / treatment); zoonotic (C. parvum causes
                     human cryptosporidiosis).
    BCOV       : 15  Supportive care only; 5-20% mortality; neonates respond
                     to oral rehydration therapy.  Frequently co-infects with
                     rotavirus.
    ROTA_A     : 10  Most common cause of neonatal calf scours worldwide;
                     self-limiting with oral rehydration; 5-15% mortality in
                     well-managed herds.

Alert level boundaries
----------------------
    Score >= 50  -> CRITICAL   (ETEC K99: IV antibiotic therapy within 2 h)
    Score >= 35  -> HIGH       (Salmonella: antibiotic + systemic workup)
    Score >= 15  -> MODERATE   (Crypto or BCoV: specific or supportive care)
    Score >=  1  -> LOW        (Rotavirus A: oral rehydration; monitor)
    Score ==  0  -> NEGATIVE   (all negative; confirm extraction control positive)

Special decision flags
----------------------
``antibiotic_indicated``
    True when ETEC K99 or Salmonella is positive.  These are the only
    pathogens in the panel for which antimicrobial therapy is the
    standard of care; flagging them prevents unnecessary antibiotic use
    for viral/parasitic co-infections.

``zoonotic_risk``
    True when Salmonella or Cryptosporidium is positive.  Both are
    notifiable zoonotic pathogens in most jurisdictions.  Personnel
    handling affected calves should wear gloves, wash hands, and avoid
    touching mucous membranes until hands are washed; calves should be
    isolated from human food-handling areas.

References:
    Cho Y & Yoon K (2014) An overview of calf diarrhea -- infectious etiology,
        diagnosis, and intervention. J Vet Sci 15(1):1-17.
        doi:10.4142/jvs.2014.15.1.1
    Gulliksen SM et al. (2009) Enteropathogens and risk factors for diarrhea
        in Norwegian dairy calves. J Dairy Sci 92(10):5057-5066.
        doi:10.3168/jds.2009-2080
    Bartels CJM et al. (2010) Prevalence of bovine respiratory disease in
        veal calves. Vet Rec 167(25):960-963. doi:10.1136/vr.c5129
    de Graaf MI et al. (1999) Bovine coronavirus N protein gene expression
        in the enteric tract.  J Virol 73(10):8485-8490.
        doi:10.1128/JVI.73.10.8485-8490.1999
    Smith DR (2012) Field diseases of beef calves in the first 45 days of
        life: enteric and respiratory diseases. Vet Clin Food Anim
        28(1):57-82. doi:10.1016/j.cvfa.2012.01.001
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


class BovEntericAlertLevel(StrEnum):
    """Bovine neonatal calf diarrhea LAMP panel alert severity.

    CRITICAL: ETEC K99 detected; near-certain mortality without IV antibiotic
        therapy within 2 hours; calves <4 days old are at highest risk.
    HIGH: Salmonella detected without ETEC K99; bacteremia and septicemia risk;
        IV antibiotic therapy and systemic workup required; zoonotic.
    MODERATE: Cryptosporidium or Bovine coronavirus detected; specific or
        supportive treatment required; zoonotic risk for Crypto.
    LOW: Rotavirus A only detected; oral rehydration therapy; self-limiting
        with good supportive care; monitor for secondary bacterial infection.
    NEGATIVE: All targets negative.  Confirm extraction control (16S rRNA)
        is positive before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class BovEntericFlags:
    """Positivity flags for one bovine neonatal enteric LAMP panel result.

    Each flag is True if the corresponding LAMP or RT-LAMP assay returned a
    positive amplification signal for that pathogen.

    Attributes:
        etec_k99: ETEC K99 / F5 fimbria (faeG gene) positive.
        salmonella: Salmonella sp. (invA gene) positive.
        crypto: Cryptosporidium parvum (18S rRNA or HSP70 gene) positive.
        bcov: Bovine coronavirus enteric pathotype (N gene) positive.
        rota_a: Bovine rotavirus A (VP6 gene) positive.
    """

    etec_k99: bool
    salmonella: bool
    crypto: bool
    bcov: bool
    rota_a: bool


@dataclass(frozen=True, slots=True)
class BovEntericAssessment:
    """Structured alert from a bovine neonatal calf diarrhea LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation.
        antibiotic_indicated: True when ETEC K99 or Salmonella is positive.
            These are the only pathogens in this panel for which antimicrobial
            therapy is the standard of care.  Prevents antibiotic overuse in
            pure viral or parasitic diarrhea cases.
        zoonotic_risk: True when Salmonella or Cryptosporidium is positive.
            Both are notifiable food-safety / zoonotic pathogens; enhanced
            hygiene precautions apply for all personnel handling affected
            calves.
    """

    flags: BovEntericFlags
    alert_level: BovEntericAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    interpretation: str
    recommended_action: str
    antibiotic_indicated: bool
    zoonotic_risk: bool

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
            "antibiotic_indicated": self.antibiotic_indicated,
            "zoonotic_risk": self.zoonotic_risk,
            "panel_flags": {
                "etec_k99": self.flags.etec_k99,
                "salmonella": self.flags.salmonella,
                "crypto": self.flags.crypto,
                "bcov": self.flags.bcov,
                "rota_a": self.flags.rota_a,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_ETEC_K99 = 50
_WEIGHT_SALMONELLA = 35
_WEIGHT_CRYPTO = 20
_WEIGHT_BCOV = 15
_WEIGHT_ROTA_A = 10


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> BovEntericAlertLevel:
    if score >= 50:
        return BovEntericAlertLevel.CRITICAL
    if score >= 35:
        return BovEntericAlertLevel.HIGH
    if score >= 15:
        return BovEntericAlertLevel.MODERATE
    if score >= 1:
        return BovEntericAlertLevel.LOW
    return BovEntericAlertLevel.NEGATIVE


def _build_text(
    flags: BovEntericFlags,
    active_targets: tuple[str, ...],
    alert_level: BovEntericAlertLevel,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    etec = flags.etec_k99
    salm = flags.salmonella
    crypto = flags.crypto
    bcov = flags.bcov
    rota = flags.rota_a

    # All negative
    if not any([etec, salm, crypto, bcov, rota]):
        return (
            "All bovine neonatal enteric targets negative. No enteric pathogens detected.",
            "Confirm an internal extraction control (e.g. 16S rRNA) is positive "
            "before reporting a clean result. Re-test if watery diarrhoea, depression, "
            "or dehydration in neonatal calves is observed within 24 h.",
        )

    # ETEC K99 positive (highest priority; dominant regardless of co-detection)
    if etec:
        co: list[str] = []
        if salm:
            co.append("Salmonella")
        if crypto:
            co.append("Cryptosporidium")
        if bcov:
            co.append("BCoV")
        if rota:
            co.append("Rotavirus A")
        co_str = "; co-detected with " + ", ".join(co) + "." if co else "."
        interp = (
            "Enterotoxigenic E. coli K99/F5 (ETEC K99) detected via the faeG assay"
            + co_str
            + " ETEC K99 causes hypersecretory diarrhoea in calves under 4 days old "
            "through heat-stable enterotoxin (STa/STb) and K99/F5 fimbrial adhesion to "
            "small-intestinal epithelium. Without intravenous antibiotic therapy and "
            "fluid replacement, mortality exceeds 90% within 24-48 h of onset."
        )
        action = (
            "CRITICAL: Start IV or oral electrolyte rehydration therapy immediately. "
            "Administer IV antibiotics (ampicillin or amoxicillin-clavulanate) within "
            "2 h of diagnosis -- do NOT wait for culture sensitivity. Isolate affected "
            "calves from unaffected littermates. Notify the farm veterinarian. Review "
            "colostrum management (ETEC K99 protection depends on passive transfer of "
            "K99-specific immunoglobulins from dam colostrum or a K99 vaccine programme). "
            "Sanitise pens with a disinfectant active against gram-negative bacteria "
            "before reuse."
        )
        return interp, action

    # Salmonella positive (no ETEC)
    if salm:
        co_str_parts = []
        if crypto:
            co_str_parts.append("Cryptosporidium")
        if bcov:
            co_str_parts.append("BCoV")
        if rota:
            co_str_parts.append("Rotavirus A")
        co_str = " Co-detected with " + ", ".join(co_str_parts) + "." if co_str_parts else ""
        interp = (
            "Salmonella sp. detected via the invA assay."
            + co_str
            + " Salmonella causes bacteraemia and septicaemia in neonatal calves with "
            "40-70% mortality without treatment. S. Dublin (cattle-adapted) and "
            "S. Typhimurium DT104 are the most common serovars in neonatal bovine "
            "scours and are notifiable zoonotic pathogens -- personnel handling "
            "affected calves are at risk of infection. Confirmatory culture and "
            "antibiotic sensitivity testing are required before selecting the "
            "antimicrobial protocol."
        )
        action = (
            "HIGH: Initiate fluid therapy immediately (IV if severe dehydration). "
            "Notify the farm veterinarian for antibiotic selection pending sensitivity "
            "results; avoid fluoroquinolones as first-line therapy in food animals "
            "without veterinary authorisation. Isolate affected calves. Report "
            "Salmonella detection to the farm herd-health programme and relevant "
            "veterinary or public-health authority (zoonotic). Enforce strict hygiene: "
            "gloves, boot covers, and handwashing for all personnel handling scouring "
            "calves. Collect fresh faecal samples for culture and sensitivity testing "
            "at a diagnostic laboratory."
        )
        return interp, action

    # Cryptosporidium positive (no ETEC, no Salmonella)
    if crypto:
        co_str_parts2 = []
        if bcov:
            co_str_parts2.append("BCoV")
        if rota:
            co_str_parts2.append("Rotavirus A")
        co_str2 = " Co-detected with " + ", ".join(co_str_parts2) + "." if co_str_parts2 else ""
        interp = (
            "Cryptosporidium parvum detected."
            + co_str2
            + " C. parvum is a zoonotic protozoan parasite; neonatal calves are the "
            "primary environmental reservoir for human cryptosporidiosis. In calves, "
            "it causes 10-50% mortality in untreated neonates and self-resolves in "
            "adequately nourished older calves. Halofuginone (lactate) is the only "
            "licensed specific treatment / prophylactic in most jurisdictions; antibiotics "
            "are ineffective."
        )
        action = (
            "Consult the farm veterinarian for halofuginone lactate administration "
            "protocol (prophylactic: 100 mcg/kg PO once daily from birth for 7 days; "
            "therapeutic: same dose for 7 days from onset of diarrhoea). Provide oral "
            "rehydration therapy. Enforce strict zoonotic hygiene -- C. parvum oocysts "
            "are chlorine-resistant; use steam cleaning or approved oocysticide. "
            "Avoid feeding waste milk from affected calves to healthy calves without "
            "heat treatment. Notify the herd veterinarian; Cryptosporidiosis is a "
            "notifiable zoonosis in some jurisdictions."
        )
        return interp, action

    # BCoV positive (no ETEC, Salmonella, or Crypto)
    if bcov:
        co_str3 = " Co-detected with Rotavirus A." if rota else ""
        interp = (
            "Bovine coronavirus (enteric pathotype) detected via the N-gene assay."
            + co_str3
            + " BCoV causes neonatal calf diarrhea and winter dysentery (adult cattle). "
            "Neonatal diarrhea from BCoV is managed with supportive care; mortality is "
            "5-20% in untreated neonates but substantially lower with early oral "
            "rehydration therapy. Antibiotics are NOT indicated for primary BCoV "
            "diarrhea."
        )
        action = (
            "Initiate oral electrolyte rehydration therapy (ORT) immediately (4-6 L/day "
            "for a 40 kg calf). Continue milk feeding alongside ORT if the calf is "
            "bright and able to stand. Do not administer antibiotics unless secondary "
            "bacterial infection is confirmed or strongly suspected clinically. Monitor "
            "for dehydration: skin-tent test and sunken eyes. Isolate affected calves "
            "to prevent pen contamination. Consult the herd veterinarian if mortality "
            "occurs or diarrhoea persists beyond 5 days."
        )
        return interp, action

    # Rotavirus A only
    if rota and not etec and not salm and not crypto and not bcov:
        interp = (
            "Bovine rotavirus A detected via the VP6 assay. Rotavirus A is the most "
            "common cause of neonatal calf scours worldwide, with 5-15% mortality in "
            "well-managed herds and higher mortality when colostrum intake was inadequate. "
            "The disease is self-limiting in calves that receive adequate colostrum and "
            "rehydration support. Antibiotics are NOT indicated for primary rotavirus "
            "diarrhea."
        )
        action = (
            "Initiate oral electrolyte rehydration therapy (ORT). Continue milk feeding "
            "alongside ORT if the calf is alert. Monitor for dehydration and secondary "
            "bacterial infection (which would increase the urgency of re-testing). "
            "Review colostrum management programme and dam vaccination history for "
            "rotavirus coverage. Most rotavirus outbreaks resolve within 5-7 days "
            "with good supportive care."
        )
        return interp, action

    # Fallback (mixed combinations not covered above)
    targets_str = ", ".join(active_targets)
    return (
        f"Detected: {targets_str}. {alert_level.value} alert.",
        "Review with farm veterinarian. Initiate supportive care and follow "
        "standard neonatal calf diarrhea protocol.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_bov_enteric_risk(flags: BovEntericFlags) -> BovEntericAssessment:
    """Compute a structured neonatal calf diarrhea alert from five-target panel flags.

    Weights each pathogen by time-to-death in untreated neonates and antibiotic
    urgency, sums to a 0-100 score, and maps it to one of five categorical alert
    levels.  Sets ``antibiotic_indicated`` when ETEC K99 or Salmonella is positive
    and ``zoonotic_risk`` when Salmonella or Cryptosporidium is positive.

    Args:
        flags: Positivity flags for all five bovine neonatal enteric targets.

    Returns:
        :class:`BovEntericAssessment` with alert level, score, interpretation,
        and recommended action.

    Example::

        from lamp_forge.bov_enteric_risk import BovEntericFlags, assess_bov_enteric_risk

        flags = BovEntericFlags(
            etec_k99=True, salmonella=False, crypto=False, bcov=False, rota_a=False
        )
        assessment = assess_bov_enteric_risk(flags)
        print(assessment.alert_level)         # BovEntericAlertLevel.CRITICAL
        print(assessment.antibiotic_indicated)  # True
    """
    score = 0
    if flags.etec_k99:
        score += _WEIGHT_ETEC_K99
    if flags.salmonella:
        score += _WEIGHT_SALMONELLA
    if flags.crypto:
        score += _WEIGHT_CRYPTO
    if flags.bcov:
        score += _WEIGHT_BCOV
    if flags.rota_a:
        score += _WEIGHT_ROTA_A

    score = max(0, min(100, score))

    target_order = [
        ("ETEC_K99", flags.etec_k99),
        ("Salmonella", flags.salmonella),
        ("Crypto", flags.crypto),
        ("BCoV", flags.bcov),
        ("RotaA", flags.rota_a),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(flags, active_targets, alert_level)

    return BovEntericAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        interpretation=interpretation,
        recommended_action=recommended_action,
        antibiotic_indicated=flags.etec_k99 or flags.salmonella,
        zoonotic_risk=flags.salmonella or flags.crypto,
    )


def flags_from_dict(data: dict[str, object]) -> BovEntericFlags:
    """Construct :class:`BovEntericFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``etec_k99``, ``salmonella``, ``crypto``, ``bcov``,
    ``rota_a`` (all bool).  Missing keys default to False; unknown keys
    are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`BovEntericFlags` instance.
    """
    return BovEntericFlags(
        etec_k99=bool(data.get("etec_k99", False)),
        salmonella=bool(data.get("salmonella", False)),
        crypto=bool(data.get("crypto", False)),
        bcov=bool(data.get("bcov", False)),
        rota_a=bool(data.get("rota_a", False)),
    )


def write_assessment_json(assessment: BovEntericAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`BovEntericAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: BovEntericAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`BovEntericAssessment` to serialise.
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
