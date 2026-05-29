"""Shared dataclasses for LAMP-Forge.

Every module passes these around instead of dicts so the pipeline contract is
explicit. Dataclasses are frozen where the value should not mutate after
construction (Primer, SpecificityHit) and mutable where downstream stages
annotate (LampPrimerSet — specificity scoring fills fields after design).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

PrimerRole = Literal["F3", "B3", "FIP", "BIP", "LF", "LB"]
"""The six primers that make up a complete LAMP set, plus the two optional loop primers.

LF and LB are technically optional — a 4-primer LAMP set (F3, B3, FIP, BIP) will
amplify, but loop primers cut reaction time roughly in half (Nagamine 2002).
LAMP-Forge designs all six by default and reports loop-primer-included sets first.
"""


@dataclass(frozen=True, slots=True)
class Primer:
    """A single oligonucleotide that participates in a LAMP set.

    Coordinates are 0-based, half-open against the consensus sequence of the
    conserved region the set was designed in. ``strand`` reflects which strand
    the primer anneals to: F* primers are sense (+), B* primers are antisense (-).

    FIP and BIP are chimeric primers — ``sequence`` is the assembled chimera
    (F1c+F2 for FIP, B1c+B2 for BIP), and ``inner_subseq`` / ``outer_subseq``
    store the two halves separately so we can compute Tm for the F1c/F2 arms
    independently (LAMP convention: both arms should be in the 60-65°C band).
    """

    role: PrimerRole
    sequence: str
    start: int
    end: int
    strand: Literal["+", "-"]
    tm: float
    gc_percent: float
    hairpin_dg: float
    homodimer_dg: float
    inner_subseq: str | None = None
    outer_subseq: str | None = None

    @property
    def length(self) -> int:
        """Primer length in nt (including both arms for chimeric FIP/BIP)."""
        return len(self.sequence)


@dataclass(slots=True)
class LampPrimerSet:
    """A full LAMP primer set (F3, B3, FIP, BIP, +/- LF, LB) for one region.

    Mutable because scoring stages (specificity, composite score) populate
    fields after the design stage creates the set.
    """

    set_id: str
    region_id: str
    f3: Primer
    b3: Primer
    fip: Primer
    bip: Primer
    lf: Primer | None = None
    lb: Primer | None = None
    f2_b2_distance: int = 0
    tm_spread: float = 0.0
    mean_tm: float = 0.0
    conservation_score: float = 0.0
    specificity_score: float = 1.0
    thermodynamic_score: float = 0.0
    composite_score: float = 0.0
    cross_reactivity_hits: list[SpecificityHit] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def primers(self) -> list[Primer]:
        """All non-None primers in canonical LAMP order."""
        ordered: list[Primer | None] = [self.f3, self.b3, self.fip, self.bip, self.lf, self.lb]
        return [p for p in ordered if p is not None]

    @property
    def has_loop_primers(self) -> bool:
        """True if both LF and LB were successfully designed."""
        return self.lf is not None and self.lb is not None


@dataclass(frozen=True, slots=True)
class ConservedRegion:
    """A stretch of the MSA flagged as suitable for LAMP target placement.

    ``mean_entropy`` is averaged across the window-smoothed per-position
    Shannon entropy (bits) over [start, end). Lower means more conserved.
    The consensus sequence is what primer design actually operates on — it
    is the IUPAC-collapsed majority sequence of the underlying MSA columns.
    """

    region_id: str
    start: int
    end: int
    mean_entropy: float
    max_entropy: float
    consensus: str

    @property
    def length(self) -> int:
        """Region length in bp."""
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class SpecificityHit:
    """A BLAST hit linking one primer to one off-target sequence.

    Flagged hits (``is_flagged=True``) exceed the user-configured identity
    and coverage thresholds and contribute to the set-level cross-reactivity
    penalty. Unflagged hits are kept for the report so reviewers can see the
    full picture.
    """

    primer_role: PrimerRole
    primer_sequence: str
    off_target_name: str
    percent_identity: float
    alignment_length: int
    primer_length: int
    coverage: float
    e_value: float
    bitscore: float
    is_flagged: bool


@dataclass(slots=True)
class PipelineConfig:
    """Parsed and validated YAML config.

    See ``config/example_config.yaml`` for the schema. We keep this as a
    plain dataclass rather than pulling in Pydantic — the surface is small
    and the YAML is already validated structurally on load.
    """

    target_name: str
    taxon_id: int | None
    accessions: list[str]
    gene: str | None
    max_sequences: int
    email: str

    off_target_dir: Path
    min_identity_threshold: float
    min_coverage_threshold: float

    window_size: int
    entropy_threshold: float
    min_region_length: int

    tm_min: float
    tm_max: float
    tm_match_tolerance: float
    gc_min: float
    gc_max: float
    hairpin_dg_threshold: float
    dimer_dg_threshold: float
    f2_b2_min: int
    f2_b2_max: int

    output_dir: Path
    top_n: int
    generate_html: bool
    generate_csv: bool

    cache_dir: Path = Path("cache")
