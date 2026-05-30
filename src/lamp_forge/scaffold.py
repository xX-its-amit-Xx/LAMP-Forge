r"""Config scaffold generator for new LAMP assay targets.

Generates a LAMP-Forge YAML config pre-populated with vertical-specific
defaults so users can go from target identification to a validated config
in one command.  Three BioVind deployment verticals are supported:

* ``oil-gas`` — industrial MIC/souring: functional genes (dsrB, mcrA, omcA)
  under strong catalytic constraint but high synonymous diversity.
* ``farm``    — on-farm biosecurity: DNA or RNA livestock virus targets
  (ASFV, PRRSV, FMDV, AIV); RNA mode enables RT-LAMP defaults.
* ``poc``     — human point-of-care: tight conservation requirements for
  clinical sensitivity across isolates.
* ``generic`` — neutral defaults; no vertical-specific tuning.

Usage::

    lamp-forge scaffold \\
      --target-name prrsv_orf7 \\
      --vertical farm \\
      --na-type rna \\
      --taxon-id 28344 \\
      --gene ORF7 \\
      --out config/prrsv_orf7.yaml

References:
    NEB WarmStart(R) RTLAMP Kit manual (E1700): 63-65 degC one-step protocol.
    Li et al. (2018) Multiplex LAMP protocol for isothermal amplification.
    Notomi T et al. (2000) Nucleic Acids Res 28(12):e63.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from textwrap import dedent


class Vertical(StrEnum):
    """BioVind deployment vertical — governs default parameter choices.

    Attributes:
        OIL_GAS: Industrial asset integrity (MIC/souring); functional marker
            genes span diverse genera and tolerate a wider entropy window.
        FARM: On-farm animal biosecurity; DNA or RNA livestock virus targets.
        POC: Human point-of-care; tight conservation for clinical sensitivity.
        GENERIC: Vertical-neutral defaults; safe starting point for any target.
    """

    OIL_GAS = "oil-gas"
    FARM = "farm"
    POC = "poc"
    GENERIC = "generic"


class NaType(StrEnum):
    """Whether the diagnostic target molecule is DNA or RNA.

    RNA targets (viruses, mRNA markers) require a reverse-transcription step
    before the isothermal amplification; one-step RT-LAMP runs at 63-65 degC
    (higher than the 60 degC floor for pure DNA-LAMP) to maintain
    reverse-transcriptase co-activity.

    Attributes:
        DNA: Double-stranded DNA target; standard LAMP at 60-65 degC.
        RNA: RNA target; one-step RT-LAMP Tm floor raised to 63 degC.
    """

    DNA = "dna"
    RNA = "rna"


@dataclass(frozen=True, slots=True)
class _VerticalDefaults:
    """Internal defaults bundle for one (Vertical, NaType) combination."""

    entropy_threshold: float
    gc_min: float
    gc_max: float
    tm_min: float
    tm_max: float
    max_sequences: int
    id_threshold: float
    cov_threshold: float
    min_region_length: int
    f2_b2_min: int
    f2_b2_max: int
    entropy_note: str
    tm_note: str


# Keyed by (Vertical, NaType).  Every combination is explicitly covered so
# that a KeyError is a programming error (caught at import time by tests).
_DEFAULTS: dict[tuple[Vertical, NaType], _VerticalDefaults] = {
    (Vertical.OIL_GAS, NaType.DNA): _VerticalDefaults(
        entropy_threshold=0.35,
        gc_min=35.0,
        gc_max=65.0,
        tm_min=60.0,
        tm_max=65.0,
        max_sequences=40,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=220,
        f2_b2_min=120,
        f2_b2_max=160,
        entropy_note=(
            "functional genes (dsrB/mcrA/omcA) are sequence-diverse "
            "despite catalytic constraint; loosen to find usable windows"
        ),
        tm_note="standard LAMP; single incubation temperature 60-65 degC",
    ),
    (Vertical.OIL_GAS, NaType.RNA): _VerticalDefaults(
        entropy_threshold=0.35,
        gc_min=35.0,
        gc_max=65.0,
        tm_min=63.0,
        tm_max=65.0,
        max_sequences=40,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=220,
        f2_b2_min=120,
        f2_b2_max=160,
        entropy_note=("functional genes; loosen entropy to find windows in diverse genera"),
        tm_note=("one-step RT-LAMP floor 63 degC for co-activity with NEB RTx / AMV-RT"),
    ),
    (Vertical.FARM, NaType.DNA): _VerticalDefaults(
        entropy_threshold=0.25,
        gc_min=35.0,
        gc_max=65.0,
        tm_min=60.0,
        tm_max=65.0,
        max_sequences=25,
        id_threshold=0.85,
        cov_threshold=0.85,
        min_region_length=220,
        f2_b2_min=120,
        f2_b2_max=150,
        entropy_note=("DNA virus target (e.g. ASFV); tight for well-conserved capsid genes"),
        tm_note="standard LAMP; single incubation temperature 60-65 degC",
    ),
    (Vertical.FARM, NaType.RNA): _VerticalDefaults(
        entropy_threshold=0.30,
        gc_min=35.0,
        gc_max=60.0,
        tm_min=63.0,
        tm_max=65.0,
        max_sequences=30,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=200,
        f2_b2_min=120,
        f2_b2_max=150,
        entropy_note=(
            "RNA virus (e.g. PRRSV/AIV/FMDV); loosen slightly for synonymous drift across genotypes"
        ),
        tm_note=(
            "one-step RT-LAMP floor 63 degC (PRRSV / avian influenza / FMDV one-step protocol)"
        ),
    ),
    (Vertical.POC, NaType.DNA): _VerticalDefaults(
        entropy_threshold=0.20,
        gc_min=40.0,
        gc_max=65.0,
        tm_min=60.0,
        tm_max=65.0,
        max_sequences=20,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=220,
        f2_b2_min=120,
        f2_b2_max=160,
        entropy_note=("clinical target; tight conservation for cross-isolate sensitivity"),
        tm_note="standard LAMP; single incubation temperature 60-65 degC",
    ),
    (Vertical.POC, NaType.RNA): _VerticalDefaults(
        entropy_threshold=0.25,
        gc_min=35.0,
        gc_max=60.0,
        tm_min=63.0,
        tm_max=65.0,
        max_sequences=25,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=200,
        f2_b2_min=120,
        f2_b2_max=150,
        entropy_note=(
            "RNA virus (e.g. SARS-CoV-2 / influenza); slightly looser "
            "for cross-lineage conservation"
        ),
        tm_note="one-step RT-LAMP floor 63 degC",
    ),
    (Vertical.GENERIC, NaType.DNA): _VerticalDefaults(
        entropy_threshold=0.25,
        gc_min=40.0,
        gc_max=65.0,
        tm_min=60.0,
        tm_max=65.0,
        max_sequences=20,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=200,
        f2_b2_min=120,
        f2_b2_max=160,
        entropy_note="general default; adjust based on target sequence diversity",
        tm_note="standard LAMP; single incubation temperature 60-65 degC",
    ),
    (Vertical.GENERIC, NaType.RNA): _VerticalDefaults(
        entropy_threshold=0.30,
        gc_min=35.0,
        gc_max=60.0,
        tm_min=63.0,
        tm_max=65.0,
        max_sequences=20,
        id_threshold=0.80,
        cov_threshold=0.80,
        min_region_length=200,
        f2_b2_min=120,
        f2_b2_max=150,
        entropy_note="RNA virus general default; loosen if too few conserved regions found",
        tm_note="one-step RT-LAMP floor 63 degC",
    ),
}


@dataclass(frozen=True, slots=True)
class ScaffoldParams:
    """Parameters for generating a LAMP-Forge scaffold config.

    Attributes:
        target_name: Short identifier used in ``target.name`` and output paths
            (e.g. ``prrsv_orf7``, ``srb_dsrB``).
        vertical: Deployment context controlling default thresholds.
        na_type: DNA or RNA target; RNA enables RT-LAMP temperature defaults.
        taxon_id: NCBI taxonomy ID for the target organism.  When ``None``,
            the generated YAML contains a ``null`` placeholder that the user
            must replace with a taxon ID or an explicit accessions list.
        gene: Target gene name as annotated in NCBI records (e.g. ``dsrB``,
            ``ORF7``).  ``None`` leaves the field as ``null``.
        email: NCBI courtesy email address written into the config.  Defaults
            to a placeholder that must be replaced before running the pipeline.
        max_sequences: Override the vertical default maximum sequence count.
            ``None`` uses the vertical default.
    """

    target_name: str
    vertical: Vertical
    na_type: NaType
    taxon_id: int | None = None
    gene: str | None = None
    email: str = "your-email@institution.edu"
    max_sequences: int | None = None


def scaffold_yaml(params: ScaffoldParams) -> str:
    """Return a YAML string for a LAMP-Forge config with vertical-appropriate defaults.

    The generated YAML is a valid ``yaml.safe_load``-parseable document.  It
    includes inline comments explaining each parameter choice so a wet-lab
    user can understand and adjust the values without reading the full docs.

    Note: if ``params.taxon_id`` is ``None``, the generated config cannot be
    validated by ``lamp-forge validate`` until the user sets either
    ``target.taxon_id`` or an explicit ``target.accessions`` list.

    Args:
        params: Scaffold parameters specifying target, vertical, and NA type.

    Returns:
        YAML text suitable for writing to a ``.yaml`` file.
    """
    key = (params.vertical, params.na_type)
    defs = _DEFAULTS[key]

    taxon_str = str(params.taxon_id) if params.taxon_id is not None else "null"
    gene_str = f'"{params.gene}"' if params.gene is not None else "null"
    max_seqs = params.max_sequences if params.max_sequences is not None else defs.max_sequences
    na_label = "RNA" if params.na_type is NaType.RNA else "DNA"

    rt_tip = ""
    if params.na_type is NaType.RNA:
        rt_tip = (
            f"\n# After design, verify RT-LAMP readiness:\n"
            f"#   lamp-forge rt-check "
            f"--input results/{params.target_name}/primer_sets.json --na-type rna"
        )

    return dedent(
        f"""\
        # LAMP-Forge config: {params.target_name}
        # Vertical: {params.vertical.value}  |  Target NA type: {na_label}
        # Generated by: lamp-forge scaffold
        #
        # Next steps:
        #   1. Replace target.email with your NCBI courtesy email.
        #   2. Set target.taxon_id (or add explicit accessions under target.accessions).
        #   3. Drop off-target FASTAs into input/off_targets/
        #      (see docs/cookbook.md for recommended panels per vertical).
        #   4. Validate: lamp-forge validate --config <this_file>
        #   5. Run:      lamp-forge run    --config <this_file>{rt_tip}

        target:
          name: {params.target_name}
          taxon_id: {taxon_str}
          gene: {gene_str}
          max_sequences: {max_seqs}
          email: {params.email}

        off_targets:
          fasta_dir: input/off_targets
          min_identity_threshold: {defs.id_threshold:.2f}
          min_coverage_threshold: {defs.cov_threshold:.2f}

        conservation:
          window_size: 30
          entropy_threshold: {defs.entropy_threshold:.2f}  # {defs.entropy_note}
          min_region_length: {defs.min_region_length}

        primer:
          tm_min: {defs.tm_min:.1f}  # {defs.tm_note}
          tm_max: {defs.tm_max:.1f}
          tm_match_tolerance: 2.0
          gc_min: {defs.gc_min:.1f}
          gc_max: {defs.gc_max:.1f}
          hairpin_dg_threshold: -2.0
          dimer_dg_threshold: -5.0
          amplicon_size:
            f2_b2_min: {defs.f2_b2_min}
            f2_b2_max: {defs.f2_b2_max}

        output:
          dir: results/{params.target_name}
          top_n: 10
          generate_html: true
          generate_csv: true
        """
    )


def write_scaffold(params: ScaffoldParams, path: Path) -> None:
    """Write scaffold YAML to ``path``, creating parent directories as needed.

    Args:
        params: Scaffold parameters; forwarded to :func:`scaffold_yaml`.
        path: Destination file path.  The file is written in UTF-8.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(scaffold_yaml(params), encoding="utf-8")
