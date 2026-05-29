"""Conservation analysis: per-position Shannon entropy + region detection.

Algorithm:

1. For each column of the alignment, compute Shannon entropy over {A, C, G, T}.
   Gaps are ignored (we don't penalise indels here — that's a separate
   coverage track and most LAMP-relevant gaps reflect alignment artefacts at
   the edges, not real biological variation).
2. Smooth with a sliding window of width ``window_size`` (default 30, ~LAMP
   primer length). The smoothed entropy at position *i* is the mean of
   columns [i, i+window_size).
3. Find contiguous runs of smoothed entropy ≤ ``entropy_threshold`` that are
   at least ``min_region_length`` long. These are the candidate target windows
   for primer design.
4. For each region, compute the consensus sequence (majority base, ties → IUPAC
   N to be conservative — primer design will reject N-containing primers).

Why entropy and not pairwise identity: entropy generalises to >2 sequences
cleanly, treats partial conservation correctly (a column that's 80%A/20%G is
less conserved than a column that's 99%A/1%G), and aligns with how
conservation is reported in marker-gene literature (Schneider & Stephens 1990).

A note on units: we report entropy in bits (log base 2). Maximum entropy at a
column with four equal bases is 2.0 bits. Common "highly conserved" thresholds
in the rRNA / marker-gene literature land around 0.2-0.4 bits.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from math import log2

import numpy as np
from Bio.Align import MultipleSeqAlignment

from lamp_forge.types import ConservedRegion

logger = logging.getLogger(__name__)

NUCLEOTIDES = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class ConservationTrack:
    """The per-position outputs of conservation analysis.

    ``raw_entropy`` and ``smoothed_entropy`` are 1-D numpy arrays of length
    equal to the alignment width. ``consensus`` is the per-column majority
    base (with N for ties or all-gap columns) as a string of the same length.
    """

    raw_entropy: np.ndarray
    smoothed_entropy: np.ndarray
    consensus: str
    coverage: np.ndarray  # fraction of non-gap rows per column


def _column_entropy(column: str) -> float:
    """Shannon entropy (bits) of a single MSA column over {A, C, G, T}.

    Gaps and ambiguous bases are excluded from the calculation. If a column
    has no canonical bases (all-gap or all-N), entropy is defined as 0 by
    convention — those columns will be excluded from primer design by the
    coverage track anyway.
    """
    counts = Counter(b.upper() for b in column if b.upper() in NUCLEOTIDES)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for n in NUCLEOTIDES:
        c = counts.get(n, 0)
        if c == 0:
            continue
        p = c / total
        entropy -= p * log2(p)
    return entropy


def _column_consensus(column: str) -> str:
    """Majority base of a column. Returns 'N' for ties or no-canonical columns."""
    counts = Counter(b.upper() for b in column if b.upper() in NUCLEOTIDES)
    if not counts:
        return "N"
    top_count = max(counts.values())
    top_bases = [b for b, c in counts.items() if c == top_count]
    if len(top_bases) > 1:
        return "N"  # tie — conservative N rather than IUPAC ambiguity
    return top_bases[0]


def _column_coverage(column: str) -> float:
    """Fraction of rows in this column that are canonical bases (not gap/N)."""
    total = len(column)
    if total == 0:
        return 0.0
    canon = sum(1 for b in column if b.upper() in NUCLEOTIDES)
    return canon / total


def compute_track(alignment: MultipleSeqAlignment, window_size: int) -> ConservationTrack:
    """Compute per-column entropy, consensus, and coverage for the alignment.

    Args:
        alignment: A multi-sequence alignment.
        window_size: Width of the smoothing window in columns. Should
            approximately match the LAMP primer length (~30 nt) so the
            smoothed track reflects "is a primer-sized stretch around here
            uniformly conserved?".

    Returns:
        A :class:`ConservationTrack` with raw + smoothed entropy, consensus
        sequence, and per-column coverage.
    """
    if window_size < 1:
        raise ValueError(f"window_size must be ≥1, got {window_size}")

    width = alignment.get_alignment_length()
    raw = np.zeros(width, dtype=np.float64)
    coverage = np.zeros(width, dtype=np.float64)
    consensus_chars: list[str] = []

    # Materialise columns once instead of N indexed lookups.
    seqs = [str(rec.seq) for rec in alignment]
    for i in range(width):
        column = "".join(s[i] for s in seqs)
        raw[i] = _column_entropy(column)
        consensus_chars.append(_column_consensus(column))
        coverage[i] = _column_coverage(column)

    # Smoothing: mean of [i, i+window_size). For the last (window_size - 1)
    # positions we shrink the window so we don't introduce phantom conservation
    # past the alignment end.
    smoothed = np.zeros(width, dtype=np.float64)
    cumsum = np.concatenate(([0.0], np.cumsum(raw)))
    for i in range(width):
        end = min(i + window_size, width)
        n = end - i
        smoothed[i] = (cumsum[end] - cumsum[i]) / n if n > 0 else 0.0

    return ConservationTrack(
        raw_entropy=raw,
        smoothed_entropy=smoothed,
        consensus="".join(consensus_chars),
        coverage=coverage,
    )


def find_conserved_regions(
    track: ConservationTrack,
    *,
    entropy_threshold: float,
    min_region_length: int,
    min_coverage: float = 0.8,
) -> list[ConservedRegion]:
    """Identify stretches of low smoothed entropy suitable for primer design.

    A position qualifies if its smoothed entropy ≤ ``entropy_threshold`` AND
    its coverage ≥ ``min_coverage``. Contiguous qualifying positions are
    merged into regions; regions shorter than ``min_region_length`` are
    discarded.

    Args:
        track: Output of :func:`compute_track`.
        entropy_threshold: Per-position entropy ceiling (bits). Lower =
            stricter. 0.2 is a sensible default for marker genes; tighten to
            0.1 for highly divergent clades, loosen to 0.4 for highly
            conserved housekeeping genes.
        min_region_length: Minimum length in bp. LAMP needs ~200 bp between
            F3 and B3 inclusive, so 200 is the floor.
        min_coverage: Minimum per-column non-gap coverage. Defaults to 0.8 —
            columns missing in >20% of sequences are excluded so we don't
            design primers against alignment artefacts.

    Returns:
        List of :class:`ConservedRegion`, sorted by start coordinate.
    """
    qualifies = (track.smoothed_entropy <= entropy_threshold) & (
        track.coverage >= min_coverage
    )

    regions: list[ConservedRegion] = []
    width = len(qualifies)
    i = 0
    region_idx = 0
    while i < width:
        if not qualifies[i]:
            i += 1
            continue
        start = i
        while i < width and qualifies[i]:
            i += 1
        end = i
        length = end - start
        if length >= min_region_length:
            mean_entropy = float(track.smoothed_entropy[start:end].mean())
            max_entropy = float(track.smoothed_entropy[start:end].max())
            consensus = track.consensus[start:end]
            regions.append(
                ConservedRegion(
                    region_id=f"R{region_idx:03d}",
                    start=start,
                    end=end,
                    mean_entropy=mean_entropy,
                    max_entropy=max_entropy,
                    consensus=consensus,
                )
            )
            region_idx += 1

    logger.info(
        "Found %d conserved regions ≥%dbp with smoothed entropy ≤%.3f bits",
        len(regions),
        min_region_length,
        entropy_threshold,
    )
    return regions


def conservation_score(region: ConservedRegion) -> float:
    """Map a region's mean entropy onto a 0-1 conservation score.

    Score = 1 - (mean_entropy / 2.0), clamped to [0, 1]. 2.0 bits is the
    theoretical maximum for {A,C,G,T} so this scales linearly. We expose this
    as its own function so the composite ranker in primer_design.py can
    re-use it.
    """
    score = 1.0 - (region.mean_entropy / 2.0)
    return float(max(0.0, min(1.0, score)))
