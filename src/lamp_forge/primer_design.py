"""LAMP primer set design.

LAMP geometry (Notomi 2000):

::

    5' ──F3── ──F2── ──F1c (rev-comp of F1)── ... ──B1c (rev-comp of B1)── ──B2── ──B3── 3'
                ↑                                                       ↑
                F2-B2 amplicon span: 120-160 bp (inner amplicon)
        Outer flanks (F3-F2 gap, B2-B3 gap): typically 0-20 bp

Six primers per set:

- **F3** anneals to F3c on the template (outer forward).
- **B3** anneals to B3c on the template (outer backward).
- **FIP** = F1c + F2 — chimeric inner forward.
- **BIP** = B1c + B2 — chimeric inner backward.
- **LF** sits between F2 and F1 on the (-) strand (optional, ~halves reaction time).
- **LB** sits between B1 and B2 on the (+) strand (optional, ~halves reaction time).

Constraints we enforce:

- Tm 60-65 °C, spread ≤ tm_match_tolerance across the set.
- GC 40-65%.
- Hairpin and homodimer ΔG above (less negative than) configured thresholds.
- F2-B2 inner amplicon between ``f2_b2_min`` and ``f2_b2_max``.
- No ambiguous bases (N) in any primer — those columns came from MSA tie-breaks
  and would produce ambiguous synthesis orders.

We use the primer3 Python bindings (``primer3-py``) for Tm + thermodynamic
calculations. primer3 is the canonical engine for primer design (Untergasser
et al. 2012) and its salt-corrected Tm is what most synthesis houses use.

This module deliberately does NOT call primer3's `designPrimers` orchestrator
— that's PCR-oriented and doesn't know about LAMP geometry. We generate
candidates ourselves and score them with primer3's thermodynamic primitives.
"""

from __future__ import annotations

import itertools
import logging
from collections.abc import Iterable
from dataclasses import dataclass

import primer3

from lamp_forge.types import ConservedRegion, LampPrimerSet, PipelineConfig, Primer, PrimerRole

logger = logging.getLogger(__name__)

# LAMP primer length bands. Notomi 2000 + Eiken's "A Guide to LAMP primer
# design" recommendations. F1c/B1c are typically 20-22, F2/B2 18-22,
# F3/B3 18-22, loop primers 18-22.
F3_LEN_RANGE = (18, 22)
B3_LEN_RANGE = (18, 22)
F2_LEN_RANGE = (18, 22)
B2_LEN_RANGE = (18, 22)
F1C_LEN_RANGE = (20, 22)
B1C_LEN_RANGE = (20, 22)
LOOP_LEN_RANGE = (18, 22)

# Geometric gaps in bp.
F3_F2_GAP = (0, 20)   # F3 ends 0-20 bp before F2 starts (5' side)
B2_B3_GAP = (0, 20)   # B3 starts 0-20 bp after B2 ends (3' side)
F2_F1_GAP = (40, 60)  # F2 ends 40-60 bp before F1 starts (inner spacing)
B1_B2_GAP = (40, 60)  # B1 ends 40-60 bp before B2 starts (inner spacing)


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (preserves case)."""
    return seq.translate(COMPLEMENT)[::-1]


def gc_percent(seq: str) -> float:
    """GC content (%) of a sequence. Returns 0.0 for empty input."""
    if not seq:
        return 0.0
    seq_u = seq.upper()
    gc = seq_u.count("G") + seq_u.count("C")
    return 100.0 * gc / len(seq_u)


def tm(seq: str) -> float:
    """Salt-corrected Tm (°C) via primer3.

    Defaults match LAMP buffer conditions: 50 mM Na+, 8 mM Mg2+, 0.8 mM dNTPs,
    primer at 0.2 µM. Adjust at your own risk — LAMP buffers vary by vendor.
    """
    return float(
        primer3.calc_tm(
            seq.upper(),
            mv_conc=50.0,
            dv_conc=8.0,
            dntp_conc=0.8,
            dna_conc=200.0,
        )
    )


def hairpin_dg(seq: str) -> float:
    """Hairpin ΔG (kcal/mol). More negative = more stable hairpin = worse."""
    res = primer3.calc_hairpin(seq.upper())
    return float(res.dg) / 1000.0  # primer3 returns cal/mol; convert to kcal/mol


def homodimer_dg(seq: str) -> float:
    """Homodimer ΔG (kcal/mol). More negative = more stable dimer = worse."""
    res = primer3.calc_homodimer(seq.upper())
    return float(res.dg) / 1000.0


def heterodimer_dg(seq_a: str, seq_b: str) -> float:
    """Heterodimer ΔG (kcal/mol) between two primers."""
    res = primer3.calc_heterodimer(seq_a.upper(), seq_b.upper())
    return float(res.dg) / 1000.0


# ----------------------------------------------------------------------------
# Candidate primer generation
# ----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _Candidate:
    """A scored single-primer candidate, before set assembly."""

    sequence: str
    start: int
    end: int
    strand: str  # "+" or "-"
    tm: float
    gc: float
    hairpin_dg: float
    homodimer_dg: float


def _slide_candidates(
    sequence: str,
    region_offset: int,
    length_range: tuple[int, int],
    strand: str,
    cfg: PipelineConfig,
) -> list[_Candidate]:
    """Slide every (length, position) window over ``sequence`` and filter.

    Returns candidates that pass single-primer constraints (Tm, GC, hairpin,
    homodimer, no Ns). Strand is used to set Primer.strand later — the
    sequence returned is already in 5'→3' orientation (rev-comp applied for
    minus strand).

    Args:
        sequence: The consensus subsequence to scan (already extracted from
            the conserved region; passing the whole consensus would be
            wasteful).
        region_offset: Absolute start of ``sequence`` within the consensus,
            so candidate ``start``/``end`` are in consensus coordinates.
        length_range: (min, max) inclusive primer length.
        strand: "+" or "-". For "-" we rev-comp the slice before evaluation.
        cfg: Pipeline config (Tm, GC, ΔG thresholds).

    Returns:
        List of passing :class:`_Candidate` objects.
    """
    min_len, max_len = length_range
    candidates: list[_Candidate] = []
    for length in range(min_len, max_len + 1):
        for i in range(0, len(sequence) - length + 1):
            window = sequence[i : i + length]
            if "N" in window.upper():
                continue
            primer_seq = window if strand == "+" else reverse_complement(window)
            primer_tm = tm(primer_seq)
            if not (cfg.tm_min <= primer_tm <= cfg.tm_max):
                continue
            primer_gc = gc_percent(primer_seq)
            if not (cfg.gc_min <= primer_gc <= cfg.gc_max):
                continue
            h_dg = hairpin_dg(primer_seq)
            if h_dg < cfg.hairpin_dg_threshold:
                continue
            d_dg = homodimer_dg(primer_seq)
            if d_dg < cfg.dimer_dg_threshold:
                continue
            candidates.append(
                _Candidate(
                    sequence=primer_seq,
                    start=region_offset + i,
                    end=region_offset + i + length,
                    strand=strand,
                    tm=primer_tm,
                    gc=primer_gc,
                    hairpin_dg=h_dg,
                    homodimer_dg=d_dg,
                )
            )
    return candidates


def _make_primer(c: _Candidate, role: PrimerRole) -> Primer:
    """Wrap a `_Candidate` into a typed `Primer` for the output set."""
    return Primer(
        role=role,
        sequence=c.sequence,
        start=c.start,
        end=c.end,
        strand=c.strand,  # type: ignore[arg-type]
        tm=c.tm,
        gc_percent=c.gc,
        hairpin_dg=c.hairpin_dg,
        homodimer_dg=c.homodimer_dg,
    )


def _make_chimeric_primer(
    role: PrimerRole,
    inner_subseq: str,
    outer_subseq: str,
    cfg: PipelineConfig,
) -> Primer | None:
    """Build a FIP (F1c+F2) or BIP (B1c+B2) chimeric primer.

    Returns None if either arm fails the per-arm thermodynamic checks (Tm in
    band, no hairpin/dimer issues) — both arms must work independently or
    the assembled primer is unlikely to behave.

    The assembled chimera's ``start``/``end`` are set to the outer arm's
    coordinates on the consensus, because that's where it anneals first; the
    inner arm only matters once it's flipped on the loop-back. This is a
    convenience for plotting; the report shows both arms explicitly.
    """
    assembled = inner_subseq + outer_subseq

    # Check each arm: Tm should be in band for both. The full chimera Tm is
    # not directly relevant — LAMP works because each arm primes
    # independently in distinct phases of the reaction.
    inner_tm = tm(inner_subseq)
    outer_tm = tm(outer_subseq)
    if not (cfg.tm_min <= inner_tm <= cfg.tm_max):
        return None
    if not (cfg.tm_min <= outer_tm <= cfg.tm_max):
        return None

    # Hairpin/dimer check the full assembled primer — this is what the
    # synthesised oligo will actually see in solution.
    h_dg = hairpin_dg(assembled)
    if h_dg < cfg.hairpin_dg_threshold:
        return None
    d_dg = homodimer_dg(assembled)
    if d_dg < cfg.dimer_dg_threshold:
        return None

    # We treat the chimera's "primer tm" as the mean of the two arms; that's
    # what enters the set-level Tm-spread check.
    return Primer(
        role=role,
        sequence=assembled,
        start=0,  # filled by caller
        end=0,    # filled by caller
        strand="+" if role == "FIP" else "-",
        tm=(inner_tm + outer_tm) / 2,
        gc_percent=gc_percent(assembled),
        hairpin_dg=h_dg,
        homodimer_dg=d_dg,
        inner_subseq=inner_subseq,
        outer_subseq=outer_subseq,
    )


# ----------------------------------------------------------------------------
# Set assembly
# ----------------------------------------------------------------------------


def _tm_spread(values: Iterable[float]) -> float:
    """Max - min over an iterable of Tm values."""
    vs = list(values)
    return max(vs) - min(vs)


def _thermodynamic_score(primer_set: LampPrimerSet, cfg: PipelineConfig) -> float:
    """Score 0-1 reflecting how well-balanced the set is thermodynamically.

    Components, each in [0, 1]:
      - Tm tightness: 1 when spread = 0, linear → 0 at ``tm_match_tolerance``.
      - Tm centering: 1 when mean Tm sits at the midpoint of [tm_min, tm_max].
      - ΔG margin: how comfortably the worst hairpin/dimer clears the
        threshold (1 = far above, 0 = right at the threshold).
    """
    tms = [p.tm for p in primer_set.primers]
    spread = _tm_spread(tms)
    tightness = max(0.0, 1.0 - spread / cfg.tm_match_tolerance)

    centre = (cfg.tm_min + cfg.tm_max) / 2
    half_width = (cfg.tm_max - cfg.tm_min) / 2
    mean_tm = sum(tms) / len(tms)
    centering = max(0.0, 1.0 - abs(mean_tm - centre) / half_width)

    worst_hairpin = min(p.hairpin_dg for p in primer_set.primers)
    worst_dimer = min(p.homodimer_dg for p in primer_set.primers)
    hp_margin = max(0.0, min(1.0, (worst_hairpin - cfg.hairpin_dg_threshold) / 5.0))
    dm_margin = max(0.0, min(1.0, (worst_dimer - cfg.dimer_dg_threshold) / 5.0))

    return (tightness + centering + hp_margin + dm_margin) / 4


def _composite_score(primer_set: LampPrimerSet) -> float:
    """Combine conservation, specificity, and thermodynamic balance.

    Weights chosen for the "single-target microbial diagnostic" use case:
    specificity (0.45) dominates because false positives are the failure mode
    that kills clinical utility. Conservation (0.35) ensures the primers
    work on real isolates. Thermo (0.20) is necessary but the constraint
    filter already enforces the hard floor, so weight is lower.
    """
    return (
        0.35 * primer_set.conservation_score
        + 0.45 * primer_set.specificity_score
        + 0.20 * primer_set.thermodynamic_score
    )


def design_sets_for_region(
    region: ConservedRegion,
    cfg: PipelineConfig,
    *,
    max_sets: int = 50,
) -> list[LampPrimerSet]:
    """Generate candidate LAMP primer sets within one conserved region.

    Approach: pick an F2 position; given F2, the F1c position is constrained
    by the F2-F1 gap; F3 sits 5' of F2 within the F3-F2 gap. Mirror for the
    B-side. Then test every (F-side, B-side) pairing whose F2-B2 distance
    falls in [f2_b2_min, f2_b2_max].

    We cap candidates per role to keep combinatorics tractable. With
    ``max_sets=50`` and ~10-20 candidates per role this returns within
    seconds even on a slow laptop.

    Args:
        region: Conserved region to design within.
        cfg: Pipeline config.
        max_sets: Hard cap on output sets per region (top sets by composite
            score are kept).

    Returns:
        List of candidate :class:`LampPrimerSet`, sorted by descending
        composite score (specificity placeholder = 1.0 here; the
        specificity stage rescales it).
    """
    seq = region.consensus
    if len(seq) < cfg.f2_b2_max + max(F3_F2_GAP[1], B2_B3_GAP[1]) * 2:
        logger.debug("Region %s too short for primer design", region.region_id)
        return []

    # Generate candidates for each role across the whole region.
    f3_cands = _slide_candidates(seq, region.start, F3_LEN_RANGE, "+", cfg)
    b3_cands = _slide_candidates(seq, region.start, B3_LEN_RANGE, "-", cfg)
    f2_cands = _slide_candidates(seq, region.start, F2_LEN_RANGE, "+", cfg)
    b2_cands = _slide_candidates(seq, region.start, B2_LEN_RANGE, "-", cfg)
    f1c_cands = _slide_candidates(seq, region.start, F1C_LEN_RANGE, "-", cfg)
    b1c_cands = _slide_candidates(seq, region.start, B1C_LEN_RANGE, "+", cfg)

    logger.info(
        "Region %s candidates: F3=%d B3=%d F2=%d B2=%d F1c=%d B1c=%d",
        region.region_id,
        len(f3_cands),
        len(b3_cands),
        len(f2_cands),
        len(b2_cands),
        len(f1c_cands),
        len(b1c_cands),
    )

    sets: list[LampPrimerSet] = []
    set_idx = 0

    # F-side combinations: F3 → F2 → F1c with geometric gaps.
    f_side_combos = []
    for f3, f2 in itertools.product(f3_cands, f2_cands):
        gap = f2.start - f3.end
        if not (F3_F2_GAP[0] <= gap <= F3_F2_GAP[1]):
            continue
        for f1c in f1c_cands:
            f1_start = f1c.start  # F1c is rev-comp; F1 sits at the same span
            f2_f1_gap = f1_start - f2.end
            if not (F2_F1_GAP[0] <= f2_f1_gap <= F2_F1_GAP[1]):
                continue
            f_side_combos.append((f3, f2, f1c))
        if len(f_side_combos) >= max_sets * 4:
            break

    b_side_combos = []
    for b3, b2 in itertools.product(b3_cands, b2_cands):
        gap = b3.start - b2.end
        if not (B2_B3_GAP[0] <= gap <= B2_B3_GAP[1]):
            continue
        for b1c in b1c_cands:
            b1_end = b1c.end
            b1_b2_gap = b2.start - b1_end
            if not (B1_B2_GAP[0] <= b1_b2_gap <= B1_B2_GAP[1]):
                continue
            b_side_combos.append((b3, b2, b1c))
        if len(b_side_combos) >= max_sets * 4:
            break

    logger.info(
        "Region %s F-side combos: %d, B-side combos: %d",
        region.region_id,
        len(f_side_combos),
        len(b_side_combos),
    )

    for (f3, f2, f1c), (b3, b2, b1c) in itertools.product(f_side_combos, b_side_combos):
        # F2-B2 inner amplicon = b2.end - f2.start (b2 is on minus strand;
        # its end is the 5' end on plus strand).
        inner_amplicon = b2.end - f2.start
        if not (cfg.f2_b2_min <= inner_amplicon <= cfg.f2_b2_max):
            continue

        # Build FIP and BIP chimeras. FIP = F1c (rev-comp F1) + F2.
        # We take the F1c candidate's sequence (already rev-comp) directly.
        fip = _make_chimeric_primer(
            "FIP", inner_subseq=f1c.sequence, outer_subseq=f2.sequence, cfg=cfg
        )
        if fip is None:
            continue
        bip = _make_chimeric_primer(
            "BIP", inner_subseq=b1c.sequence, outer_subseq=b2.sequence, cfg=cfg
        )
        if bip is None:
            continue

        f3_primer = _make_primer(f3, "F3")
        b3_primer = _make_primer(b3, "B3")

        # Set Tm spread across F3, B3, FIP, BIP.
        all_tms = [f3.tm, b3.tm, fip.tm, bip.tm]
        spread = _tm_spread(all_tms)
        if spread > cfg.tm_match_tolerance:
            continue

        # Loop primers: scan the inter-arm space for LF (between F1 and F2 on
        # minus strand) and LB (between B1 and B2 on plus strand).
        lf_region_start = f2.end - region.start
        lf_region_end = f1c.start - region.start
        lf_seq = seq[lf_region_start:lf_region_end] if lf_region_end > lf_region_start else ""
        lf_cands = (
            _slide_candidates(lf_seq, f2.end, LOOP_LEN_RANGE, "-", cfg) if lf_seq else []
        )

        lb_region_start = b1c.end - region.start
        lb_region_end = b2.start - region.start
        lb_seq = seq[lb_region_start:lb_region_end] if lb_region_end > lb_region_start else ""
        lb_cands = (
            _slide_candidates(lb_seq, b1c.end, LOOP_LEN_RANGE, "+", cfg) if lb_seq else []
        )

        lf = _make_primer(lf_cands[0], "LF") if lf_cands else None
        lb = _make_primer(lb_cands[0], "LB") if lb_cands else None

        primer_set = LampPrimerSet(
            set_id=f"{region.region_id}_S{set_idx:03d}",
            region_id=region.region_id,
            f3=f3_primer,
            b3=b3_primer,
            fip=fip,
            bip=bip,
            lf=lf,
            lb=lb,
            f2_b2_distance=inner_amplicon,
            tm_spread=spread,
            mean_tm=sum(all_tms) / len(all_tms),
        )
        primer_set.thermodynamic_score = _thermodynamic_score(primer_set, cfg)
        primer_set.composite_score = _composite_score(primer_set)
        sets.append(primer_set)
        set_idx += 1

    sets.sort(key=lambda s: s.composite_score, reverse=True)
    return sets[:max_sets]


def design_all(
    regions: list[ConservedRegion],
    cfg: PipelineConfig,
    *,
    max_sets_per_region: int = 50,
) -> list[LampPrimerSet]:
    """Design primer sets across every conserved region and return them all.

    Conservation score is filled here (specificity is filled later by
    :mod:`lamp_forge.specificity`).
    """
    from lamp_forge.conserve import conservation_score

    all_sets: list[LampPrimerSet] = []
    for region in regions:
        c_score = conservation_score(region)
        sets = design_sets_for_region(region, cfg, max_sets=max_sets_per_region)
        for s in sets:
            s.conservation_score = c_score
            s.composite_score = _composite_score(s)
        logger.info(
            "Region %s: designed %d primer sets (mean entropy %.3f bits)",
            region.region_id,
            len(sets),
            region.mean_entropy,
        )
        all_sets.extend(sets)
    return all_sets
