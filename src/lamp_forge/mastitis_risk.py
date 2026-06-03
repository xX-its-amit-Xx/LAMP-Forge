"""Bovine mastitis LAMP panel risk assessment engine.

Interprets the positivity flags from a BioVind-style portable bovine mastitis
LAMP panel as a structured alert, guiding on-farm treatment and herd-management
decisions before milk culture results are available.

Bovine mastitis is the most economically significant disease in dairy cattle,
causing approximately USD 2 billion in annual losses in the United States alone
(Hogeveen et al. 2011).  Differentiating contagious pathogens
(*Staphylococcus aureus*, *Streptococcus agalactiae*) from environmental ones
(*Streptococcus uberis*, *Escherichia coli*) at the time of first clinical
signs is the single most important decision-shaping factor in mastitis
management:

* Contagious pathogens spread cow-to-cow through the milking process and require
  segregation, targeted dry-cow therapy, or culling.  Early identification
  prevents herd-level spread.
* Environmental pathogens arise from bedding, manure, and water; they indicate
  husbandry-level deficiencies rather than a milking-hygiene failure.
* *E. coli* mastitis is frequently self-limiting and produces endotoxin-driven
  acute illness; antimicrobials are often withheld in favour of supportive care.

A portable multiplex LAMP panel resolving all four major mastitis agents
from a milking-time milk sample in under 60 minutes enables same-visit treatment
decisions -- the core value proposition for BioVind's farm diagnostic vertical.

Pathogens covered
-----------------

    Target  Gene  Pathogen                  Transmission  NA type
    SAUR    nuc   Staphylococcus aureus     Contagious    DNA
    SAGA    cfb   Streptococcus agalactiae  Contagious    DNA
    SUBE    sue   Streptococcus uberis      Environmental DNA
    ECOLI   uidA  Escherichia coli          Environmental DNA

Assay notes:

* ``nuc`` encodes the thermostable nuclease (thermonuclease) of *S. aureus*.
  The nuc target is chromosomal, present in all *S. aureus* strains including
  methicillin-resistant (MRSA), and absent from coagulase-negative staphylococci
  (CoNS).  It is the gold-standard PCR/LAMP target for bovine *S. aureus*
  mastitis (Brakstad et al. 1992).
* ``cfb`` encodes the CAMP (Christie-Atkins-Munch-Petersen) reaction factor,
  which is unique to *S. agalactiae* among streptococcal species, making it
  exquisitely specific for *Streptococcus agalactiae* in milk samples.
* ``sue`` is the *S. uberis* 16S-rRNA-adjacent repeat locus used in LAMP and
  PCR assays for environmental streptococcal mastitis (Phuektes et al. 2001).
* ``uidA`` encodes beta-D-glucuronidase and is the same chromosomal marker
  used in the LAMP-Forge UTI panel for *E. coli*.  It is absent from all
  other coliforms (Salmonella, Klebsiella, Proteus).

Scoring model
-------------
Base weights reflect contagiousness, treatment-success rate, and the clinical
cost of a delayed or wrong empiric decision:

    SAUR : 40  S. aureus is the most important contagious mastitis pathogen.
               Intramammary treatment success for subclinical S. aureus mastitis
               is only 20-30%; biofilm-forming strains are untreatable with
               intramammary antibiotics; culling is often the recommended
               action.  Every case of S. aureus mastitis detected late risks
               becoming a herd reservoir (carrier cows shed intermittently at
               milking; teat-liners become fomites).
    SAGA : 35  S. agalactiae is an obligate udder pathogen (does not persist
               in the environment); it is, paradoxically, one of the most
               CURABLE causes of mastitis if treated promptly with intramammary
               penicillin or amoxicillin.  However, a single missed carrier cow
               in the milking herd causes continuous cow-to-cow transmission;
               eradication programmes (identify by LAMP, treat or cull) are
               effective only when all carriers are identified simultaneously.
    SUBE : 20  S. uberis is the leading environmental mastitis pathogen in
               housed dairy herds (Europe/UK) and a major cause of dry-period
               mastitis.  Treatment response is moderate (40-60%);
               recurrence within the same lactation is common.  Unlike
               contagious mastitis, husbandry changes (deep-litter management,
               cubicle hygiene, teat-spray protocol) are the primary control
               levers.
    ECOLI: 15  E. coli causes the majority of acute peracute mastitis
               episodes (sudden onset, systemic toxaemia, reduced milk yield).
               The pathology is endotoxin-driven; most cases resolve
               spontaneously in 2-5 days without antimicrobials; NSAIDs +
               fluid therapy are first-line supportive care.  Antimicrobials
               are indicated only in severe cases with fever and systemic signs.

Alert level boundaries
----------------------
    CRITICAL  : score >= 70  (contagious co-infection or polymicrobial outbreak)
    HIGH      : score >= 35  (single contagious pathogen detected)
    MODERATE  : score >= 20  (environmental streptococcus, treatable)
    LOW       : score >= 10  (E. coli alone; often self-limiting)
    NEGATIVE  : score ==  0  (all targets negative)

Special decision flags
----------------------
``contagious_pathogen_detected`` is True when *S. aureus* (SAUR) or
*S. agalactiae* (SAGA) is positive.  Contagious mastitis pathogens require
immediate herd-management action (segregation of affected cow during milking,
review of milking-machine hygiene, consideration of blanket dry-cow therapy)
beyond just treating the individual case.

``contagious_coinfection`` is True when *both* SAUR and SAGA are positive.
This pattern is unusual in a single cow but can occur in a herd with both
endemic *S. aureus* and circulating *S. agalactiae*.  Combined scoring reaches
75, triggering CRITICAL.

``s_aureus_biofilm_risk`` is True when SAUR is positive.  The majority of
bovine *S. aureus* mastitis isolates form biofilms on the teat-canal mucosa
and within the alveolar epithelium, which dramatically reduces the penetration
and efficacy of intramammary antibiotics.  Culling of chronic or recurring
*S. aureus* cases is the economically dominant strategy for the majority of
herds with SCC > 500,000 cells/mL in the affected quarter.

References:
----------
Hogeveen H et al. (2011) Economic aspects of mastitis: new developments.
    N Z Vet J 59(1):16-23. doi:10.1080/00480169.2011.552nucleotide541
Brakstad OG et al. (1992) Detection of Staphylococcus aureus by polymerase
    chain reaction amplification of the nuc gene.
    J Clin Microbiol 30(7):1654-1660. doi:10.1128/jcm.30.7.1654-1660.1992
Phuektes P et al. (2001) Multiplex polymerase chain reaction as a
    mastitis screening test for Staphylococcus aureus, Streptococcus agalactiae,
    Streptococcus dysgalactiae, and Streptococcus uberis.
    J Dairy Sci 84(5):1140-1148. doi:10.3168/jds.S0022-0302(01)74573-6
Bradley AJ & Green MJ (2001) Adaptation and survival of bovine somatic cells
    in the presence of Streptococcus uberis.
    J Dairy Sci 84(5):1099-1108. doi:10.3168/jds.S0022-0302(01)74569-4
Erskine RJ et al. (2002) Intramammary antibiotic therapy.
    J Am Vet Med Assoc 221(10):1395-1402. doi:10.2460/javma.2002.221.1395
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


class MastitisAlertLevel(StrEnum):
    """Bovine mastitis LAMP panel alert severity.

    CRITICAL: Multiple contagious pathogens or polymicrobial co-infection;
        immediate veterinary assessment and herd-level management action.
    HIGH: Single major contagious pathogen (S. aureus or S. agalactiae);
        cow-level and milking-hygiene intervention required.
    MODERATE: Environmental streptococcus (S. uberis); targeted intramammary
        therapy and bedding/hygiene review.
    LOW: E. coli alone; supportive care, antimicrobials contingent on severity.
    NEGATIVE: All targets negative; confirm sample adequacy with a positive
        16S rRNA internal control before reporting a clean result.
    """

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"


@dataclass(frozen=True, slots=True)
class MastitisPanelFlags:
    """Positivity flags for one bovine mastitis LAMP panel result.

    Each flag is True if the corresponding LAMP assay returned a positive
    amplification signal for that pathogen in the milk sample.

    Attributes:
        saur: Staphylococcus aureus (nuc thermostable nuclease gene) positive.
        saga: Streptococcus agalactiae (cfb CAMP factor gene) positive.
        sube: Streptococcus uberis (sue repeat locus) positive.
        ecoli: Escherichia coli (uidA beta-D-glucuronidase gene) positive.
    """

    saur: bool
    saga: bool
    sube: bool
    ecoli: bool


@dataclass(frozen=True, slots=True)
class MastitisRiskAssessment:
    """Structured alert from a bovine mastitis LAMP panel result.

    Attributes:
        flags: Pathogen-positivity flags used to compute this assessment.
        alert_level: Categorical alert level (CRITICAL -> NEGATIVE).
        alert_score: Numeric 0-100 summary (higher = more urgent).
        active_targets: Labels of positive targets in canonical order.
        contagious_pathogen_detected: True when S. aureus or S. agalactiae
            is positive; triggers herd-level management beyond individual-cow
            treatment.
        contagious_coinfection: True when both S. aureus and S. agalactiae
            are positive simultaneously; unusual in one cow but indicates
            herd-level co-circulation of two contagious pathogens.
        s_aureus_biofilm_risk: True when S. aureus is positive; biofilm-forming
            strains have low intramammary antibiotic penetration; culling should
            be evaluated.
        interpretation: Concise interpretation of the detection pattern.
        recommended_action: Field-actionable recommendation for the farmer
            or attending veterinarian.
    """

    flags: MastitisPanelFlags
    alert_level: MastitisAlertLevel
    alert_score: int
    active_targets: tuple[str, ...]
    contagious_pathogen_detected: bool
    contagious_coinfection: bool
    s_aureus_biofilm_risk: bool
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
            "contagious_pathogen_detected": self.contagious_pathogen_detected,
            "contagious_coinfection": self.contagious_coinfection,
            "s_aureus_biofilm_risk": self.s_aureus_biofilm_risk,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "panel_flags": {
                "saur": self.flags.saur,
                "saga": self.flags.saga,
                "sube": self.flags.sube,
                "ecoli": self.flags.ecoli,
            },
        }


# ---------------------------------------------------------------------------
# Scoring constants (see module docstring for rationale)
# ---------------------------------------------------------------------------

_WEIGHT_SAUR = 40
_WEIGHT_SAGA = 35
_WEIGHT_SUBE = 20
_WEIGHT_ECOLI = 15

_CRITICAL_THRESHOLD = 70
_HIGH_THRESHOLD = 35
_MODERATE_THRESHOLD = 20
_LOW_THRESHOLD = 10

_TARGET_DESCRIPTION: dict[str, str] = {
    "SAUR": "Staphylococcus aureus (nuc) -- contagious; biofilm-forming; low treatment success",
    "SAGA": "Streptococcus agalactiae (cfb) -- contagious; eradicable with targeted penicillin",
    "SUBE": "Streptococcus uberis (sue) -- environmental; most common mastitis in housed herds",
    "ECOLI": "Escherichia coli (uidA) -- environmental; acute peracute; usually self-limiting",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _alert_level_from_score(score: int) -> MastitisAlertLevel:
    if score >= _CRITICAL_THRESHOLD:
        return MastitisAlertLevel.CRITICAL
    if score >= _HIGH_THRESHOLD:
        return MastitisAlertLevel.HIGH
    if score >= _MODERATE_THRESHOLD:
        return MastitisAlertLevel.MODERATE
    if score >= _LOW_THRESHOLD:
        return MastitisAlertLevel.LOW
    return MastitisAlertLevel.NEGATIVE


def _build_text(
    flags: MastitisPanelFlags,
    active_targets: tuple[str, ...],
    alert_level: MastitisAlertLevel,
    contagious_coinfection: bool,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detection pattern.

    Args:
        flags: Pathogen positivity flags.
        active_targets: Positive target labels in canonical order.
        alert_level: Pre-computed alert level.
        contagious_coinfection: True when both SAUR and SAGA are positive.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    any_positive = any([flags.saur, flags.saga, flags.sube, flags.ecoli])

    if not any_positive:
        return (
            "All bovine mastitis targets negative. No major mastitis pathogens detected "
            "in this milk sample.",
            "Confirm a 16S rRNA internal control is positive before reporting a "
            "clean result. If clinical signs (abnormal milk, swollen quarter, pain) "
            "persist with a negative LAMP result, consider mycoplasma or prototheca "
            "as alternative causes, or repeat sampling from a fresh strip.",
        )

    target_str = ", ".join(active_targets)
    n_total = len(active_targets)

    # Both contagious pathogens together
    if contagious_coinfection:
        sube_note = (
            " S. uberis also detected: environmental co-contamination likely; "
            "review cubicle hygiene independently."
            if flags.sube
            else ""
        )
        ecoli_note = (
            " E. coli also detected: monitor affected cow for systemic signs "
            "(fever, depression, reduced feed intake) and initiate NSAIDs if required."
            if flags.ecoli
            else ""
        )
        interp = (
            "Both S. aureus (nuc) and S. agalactiae (cfb) detected simultaneously. "
            "This pattern indicates co-circulation of two contagious mastitis pathogens "
            "in the milking herd. S. aureus and S. agalactiae both spread cow-to-cow "
            "during milking via contaminated teat cups. Combined detection in a single "
            "sample may represent a polymicrobial infection of one quarter, or represent "
            "carry-over from two positive cows in the same sampling session. Herd-level "
            f"investigation is required to identify all carriers.{sube_note}{ecoli_note}"
        )
        action = (
            "IMMEDIATE: Segregate the LAMP-positive cow to the end of the milking order "
            "or to a separate milking unit. Conduct a whole-herd LAMP screen (bulk-tank "
            "LAMP or individual cow composite samples) to identify all S. aureus and "
            "S. agalactiae carriers. Notify the herd veterinarian. For S. agalactiae: "
            "intramammary penicillin G achieves >90% bacteriological cure if administered "
            "promptly. For S. aureus: evaluate individual cow history (SCC, lactation "
            "number, quarter chronicity) -- culling is often the most cost-effective "
            "decision for cows with SCC > 500,000 cells/mL or repeat infections. "
            "Disinfect and back-flush teat cups after each positive cow."
        )
        return interp, action

    # S. aureus alone or dominant
    if flags.saur and not flags.saga:
        env_note = ""
        if flags.sube:
            env_note += (
                " S. uberis also detected: secondary environmental mastitis concurrent "
                "with the contagious infection; treat the S. uberis component with "
                "intramammary amoxicillin or pirlimycin per veterinary formulary."
            )
        if flags.ecoli:
            env_note += (
                " E. coli also detected: monitor for endotoxaemia (fever, reduced milk, "
                "depression); supportive NSAIDs + fluid therapy if systemic signs develop."
            )
        interp = (
            "S. aureus (nuc) detected. S. aureus is the most economically important "
            "contagious mastitis pathogen in dairy herds. It colonises the teat-canal "
            "mucosa and alveolar epithelium, where biofilm formation protects it from "
            "both host defences and intramammary antibiotic penetration. Treatment success "
            "for subclinical S. aureus mastitis averages only 20-30% with current "
            "intramammary products; the probability of cure falls further in older cows, "
            "high-SCC quarters, and chronic infections (> 60 days somatic cell count "
            f"elevation). Herd-level spread occurs via contaminated teat cups.{env_note}"
        )
        action = (
            "Segregate the cow to the end of the milking order. Evaluate cure probability: "
            "cows with SCC < 200,000 cells/mL in the affected quarter, first or second "
            "lactation, and < 3 weeks of elevated SCC have the best prognosis for "
            "intramammary therapy. Consult the herd veterinarian to decide between "
            "intramammary treatment (dry-off at end of lactation preferred over in-lactation "
            "therapy for subclinical cases) and culling. Conduct a whole-herd LAMP screen "
            "to identify all S. aureus-positive cows. Enforce pre- and post-milking teat "
            "disinfection and check teat-cup liner integrity and replacement schedule."
        )
        return interp, action

    # S. agalactiae alone or dominant
    if flags.saga and not flags.saur:
        env_note = ""
        if flags.sube or flags.ecoli:
            env_note = (
                " Additional environmental pathogen(s) detected: manage S. agalactiae "
                "eradication and the environmental mastitis programme in parallel."
            )
        interp = (
            "S. agalactiae (cfb) detected. S. agalactiae is an obligate udder pathogen "
            "that does not persist in the environment; cow-to-cow transmission during "
            "milking is the sole route of spread. This makes eradication feasible when "
            "all carriers are identified and treated simultaneously. S. agalactiae is "
            "exquisitely sensitive to penicillin G; bacteriological cure rates exceed "
            "90% when intramammary penicillin G is administered at the correct dose and "
            "duration. However, if only some cows in the herd are treated while carriers "
            f"remain, the herd-level cycle of transmission continues.{env_note}"
        )
        action = (
            "Initiate intramammary penicillin G therapy per veterinary prescription "
            "(typically 3-5 consecutive milkings). Segregate the positive cow during "
            "the treatment course. Conduct a whole-herd LAMP screen to identify ALL "
            "S. agalactiae carriers -- only simultaneous treatment of all positives "
            "achieves eradication. Review teat-cup liner hygiene and post-milking teat "
            "disinfection. S. agalactiae eradication is a realistic goal for herds "
            "with a robust identification-and-treat programme; coordinate with the "
            "herd veterinarian to design the programme."
        )
        return interp, action

    # S. uberis alone (no contagious pathogen)
    if flags.sube and not flags.saur and not flags.saga:
        ecoli_note = (
            " E. coli also detected: two environmental pathogens co-present; "
            "review bedding management, cubicle hygiene, and teat-end condition urgently."
            if flags.ecoli
            else ""
        )
        interp = (
            "S. uberis (sue) detected. S. uberis is the leading cause of environmental "
            "mastitis in housed dairy herds and the most common cause of dry-period "
            "new infections in UK and European herds. It colonises the teat canal from "
            "bedding (particularly organic matter: straw, sawdust, sand), teat-end "
            "contamination, and standing water. Treatment response is moderate "
            "(40-60% bacteriological cure at 14 days), and reinfection within the same "
            f"lactation is common.{ecoli_note}"
        )
        action = (
            "Initiate intramammary therapy per veterinary prescription (amoxicillin, "
            "pirlimycin, or cephalosporin-based intramammary products commonly used; "
            "duration per product SPC). Assess teat-end condition for hyperkeratosis "
            "or teat-end roughness that increases colonisation risk. Review cubicle "
            "management: replace or treat contaminated bedding; ensure correct lying "
            "surface depth; clean and disinfect before re-filling. If S. uberis "
            "mastitis is a recurring herd-level problem, assess the dry-period "
            "management protocol (dry-cow therapy + teat sealant) and post-calving "
            "management period."
        )
        return interp, action

    # E. coli alone
    if flags.ecoli and not flags.saur and not flags.saga and not flags.sube:
        interp = (
            "E. coli (uidA) detected. E. coli mastitis is an environmental disease "
            "caused by faecal contamination of the teat end; it typically presents as "
            "an acute peracute mastitis (sudden onset, swollen hot quarter, watery "
            "milk, systemic signs). The pathology is driven by lipopolysaccharide (LPS) "
            "endotoxin from lysed E. coli; the majority of uncomplicated clinical "
            "cases resolve without antimicrobial therapy within 2-5 days, provided "
            "that frequent stripping of the affected quarter is performed to remove "
            "toxin-laden milk."
        )
        action = (
            "Strip the affected quarter frequently (every 2-4 h) to remove endotoxin. "
            "Administer NSAIDs (meloxicam or flunixin meglumine) for analgesia and "
            "anti-endotoxaemia per veterinary prescription. Provide IV fluid therapy "
            "if systemic signs are severe (recumbency, heart rate > 100 bpm, "
            "temperature > 39.5 degC or < 38.0 degC). Antimicrobials are indicated "
            "only in severe or bacteraemic cases -- consult veterinarian before "
            "prescribing. Improve hygiene in the calving and lying areas to reduce "
            "teat-end contamination with coliform bacteria."
        )
        return interp, action

    # Generic multi-target fallback (mixed patterns not covered above)
    interp = (
        f"Bovine mastitis LAMP panel: {n_total} target(s) detected -- {target_str}. "
        f"Alert level: {alert_level.value}."
    )
    action = (
        "Consult the herd veterinarian. Apply treatment appropriate to each detected "
        "pathogen (contagious pathogens require herd-level investigation; environmental "
        "pathogens require husbandry review). Strip affected quarters and monitor "
        "systemic signs closely."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def assess_mastitis_risk(flags: MastitisPanelFlags) -> MastitisRiskAssessment:
    """Compute a structured bovine mastitis alert from four-pathogen panel flags.

    Weights each pathogen by its contagiousness, treatment-success rate, and
    economic impact; sums to a 0-100 score; maps to one of five categorical
    alert levels.  Sets ``contagious_coinfection`` when both *S. aureus*
    and *S. agalactiae* are positive simultaneously.

    Args:
        flags: Positivity flags for all four bovine mastitis panel targets.

    Returns:
        :class:`MastitisRiskAssessment` with alert level, score,
        interpretation, and recommended action.

    Example::

        from lamp_forge.mastitis_risk import MastitisPanelFlags, assess_mastitis_risk

        flags = MastitisPanelFlags(saur=True, saga=False, sube=False, ecoli=False)
        assessment = assess_mastitis_risk(flags)
        print(assessment.alert_level)                  # MastitisAlertLevel.HIGH
        print(assessment.s_aureus_biofilm_risk)        # True
        print(assessment.contagious_pathogen_detected) # True
    """
    score = 0
    if flags.saur:
        score += _WEIGHT_SAUR
    if flags.saga:
        score += _WEIGHT_SAGA
    if flags.sube:
        score += _WEIGHT_SUBE
    if flags.ecoli:
        score += _WEIGHT_ECOLI

    score = max(0, min(100, score))

    target_order = [
        ("SAUR", flags.saur),
        ("SAGA", flags.saga),
        ("SUBE", flags.sube),
        ("ECOLI", flags.ecoli),
    ]
    active_targets = tuple(label for label, pos in target_order if pos)

    contagious_pathogen_detected = flags.saur or flags.saga
    contagious_coinfection = flags.saur and flags.saga
    s_aureus_biofilm_risk = flags.saur

    alert_level = _alert_level_from_score(score)
    interpretation, recommended_action = _build_text(
        flags, active_targets, alert_level, contagious_coinfection
    )

    return MastitisRiskAssessment(
        flags=flags,
        alert_level=alert_level,
        alert_score=score,
        active_targets=active_targets,
        contagious_pathogen_detected=contagious_pathogen_detected,
        contagious_coinfection=contagious_coinfection,
        s_aureus_biofilm_risk=s_aureus_biofilm_risk,
        interpretation=interpretation,
        recommended_action=recommended_action,
    )


def flags_from_dict(data: dict[str, object]) -> MastitisPanelFlags:
    """Construct :class:`MastitisPanelFlags` from a plain dict (e.g. parsed JSON).

    Expected keys: ``saur``, ``saga``, ``sube``, ``ecoli`` (all bool).
    Missing keys default to False; unknown keys are ignored.

    Args:
        data: Dict with optional bool values for each target.

    Returns:
        :class:`MastitisPanelFlags` instance.
    """
    return MastitisPanelFlags(
        saur=bool(data.get("saur", False)),
        saga=bool(data.get("saga", False)),
        sube=bool(data.get("sube", False)),
        ecoli=bool(data.get("ecoli", False)),
    )


def write_assessment_json(assessment: MastitisRiskAssessment, path: Path) -> None:
    """Write the assessment to a JSON file.

    Args:
        assessment: :class:`MastitisRiskAssessment` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(assessment.to_dict(), fh, indent=2)
        fh.write("\n")


def write_assessment_csv(assessment: MastitisRiskAssessment, path: Path) -> None:
    """Write a flat key-value CSV summary of the assessment.

    Args:
        assessment: :class:`MastitisRiskAssessment` to serialise.
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
