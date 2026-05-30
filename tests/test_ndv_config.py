"""Validate the NDV M-gene on-farm biosecurity config (config/ndv_M_gene.yaml).

Semantic checks specific to Newcastle Disease Virus (NDV) via the M gene:
- Avian orthoavulavirus 1 anchor (NCBI taxon 11234)
- Entropy threshold appropriate for a moderately conserved RNA-virus gene
- GC window covers NDV M-gene GC content (~47%) and inter-genotype diversity
- RT-LAMP Tm floor (tm_min >= 63 degC) required for one-step RNA detection
- Tight off-target identity threshold (>= 0.83) for avian paramyxovirus family

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_NDV_CONFIG = CONFIG_DIR / "ndv_M_gene.yaml"


@pytest.fixture
def ndv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the NDV M-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_NDV_CONFIG)


def test_ndv_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_NDV_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestNdvMGeneConfig:
    """Semantic checks for the NDV M-gene on-farm biosecurity config.

    Newcastle Disease Virus (NDV; Avian orthoavulavirus 1) is a WOAH-listed
    poultry pathogen with near-100% mortality in velogenic outbreaks.  The M
    (matrix protein) gene is the most conserved coding region across all Class I
    and Class II genotypes and the basis of OIE-endorsed molecular assays.
    NDV is a negative-sense ssRNA paramyxovirus, so detection requires RT-LAMP
    with a one-step reaction at 63-65 degC.
    """

    def test_target_name(self, ndv_cfg: PipelineConfig) -> None:
        assert ndv_cfg.target_name == "ndv_M_gene"

    def test_gene_is_M(self, ndv_cfg: PipelineConfig) -> None:
        """M (matrix protein) must be the target gene.

        The M gene is more conserved across all NDV pathotypes (velogenic,
        mesogenic, lentogenic) and Class II genotypes than the F gene, which
        harbours the cleavage-site motif used for pathotype determination.
        M-gene-based assays detect any NDV regardless of virulence, which
        is the correct scope for a first-pass diagnostic screen.
        """
        assert ndv_cfg.gene == "M"

    def test_taxon_is_newcastle_disease_virus(self, ndv_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Newcastle disease virus (NCBI taxon 11234).

        NCBI taxon 11234 covers Avian orthoavulavirus 1 (the current ICTV
        name) under its historical Newcastle disease virus accession set.
        Sequences span Class I (single genotype) and Class II (genotypes
        I-XXI, including the virulent genotypes responsible for active global
        outbreaks) -- the full diversity needed for a pan-NDV assay.
        """
        assert ndv_cfg.taxon_id == 11234

    def test_entropy_threshold_appropriate_for_rna_virus(
        self, ndv_cfg: PipelineConfig
    ) -> None:
        """M-gene entropy threshold must reflect RNA-virus inter-genotype variability.

        The NDV M gene is under structural constraint (matrix-layer assembly),
        so it is more conserved than the F gene across genotypes, but RNA-virus
        synonymous drift is still higher than for bacterial housekeeping genes.
        Thresholds below 0.20 bits are too strict and will yield no windows
        across the 20+ genotypes; above 0.40 bits admits positions too variable
        to support specific primers across the full genotype range.
        """
        assert ndv_cfg.entropy_threshold >= 0.20, (
            "NDV M gene diverges across Class II genotypes (I-XXI); "
            "entropy_threshold < 0.20 is too strict and will yield no usable windows"
        )
        assert ndv_cfg.entropy_threshold <= 0.40, (
            "NDV M gene is under structural assembly constraint; "
            "entropy_threshold > 0.40 admits overly variable positions"
        )

    def test_gc_range_covers_ndv_diversity(self, ndv_cfg: PipelineConfig) -> None:
        """GC window must span NDV M-gene GC content (~47%) across genotype diversity.

        The NDV M gene is ~47% GC in the OIE reference strain (LaSota, Class II
        genotype II).  Inter-genotype variability shifts GC composition +/- 5
        percentage points.  The design window must cover 40-55% GC to produce
        primers that work across velogenic outbreak strains (genotypes VII, XXI)
        and classical vaccine strains (genotypes II).
        """
        assert ndv_cfg.gc_min <= 42.0, (
            "gc_min must be <= 42% to accommodate NDV genotypes with lower GC content"
        )
        assert ndv_cfg.gc_max >= 55.0, (
            "gc_max must be >= 55% to accommodate NDV genotypes with higher GC content"
        )

    def test_gc_window_wide_enough_for_genotype_diversity(
        self, ndv_cfg: PipelineConfig
    ) -> None:
        """GC window must span NDV inter-genotype GC diversity."""
        assert ndv_cfg.gc_max - ndv_cfg.gc_min >= 15.0, (
            "NDV genotypes span significant GC diversity; "
            "gc window must be at least 15 units wide to support multi-genotype primer design"
        )

    def test_max_sequences_spans_genotype_diversity(self, ndv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover Class I and the major Class II genotypes."""
        assert ndv_cfg.max_sequences >= 30, (
            "NDV has 20+ Class II genotypes; at least 30 sequences needed for "
            "cross-genotype conservation analysis"
        )

    def test_rt_lamp_tm_floor_is_set(self, ndv_cfg: PipelineConfig) -> None:
        """NDV is an RNA virus; tm_min must meet the RT-LAMP one-step co-activity floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx at 63-65 degC) requires
        that core primers have Tm >= 63 degC to maintain reverse-transcriptase
        co-activity.  This is the same floor used for PRRSV (Recipe 11), AIV
        (Recipe 12), and FMDV (Recipe 14).
        """
        assert ndv_cfg.tm_min >= 63.0, (
            "NDV is a negative-sense ssRNA paramyxovirus: tm_min must be >= 63 degC "
            "for one-step RT-LAMP co-activity with NEB RTx / Bst 2.0 WarmStart"
        )

    def test_tm_max_within_rt_lamp_window(self, ndv_cfg: PipelineConfig) -> None:
        """Tm max must stay within the RT-LAMP optimal window (63-65 degC)."""
        assert ndv_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_tight_specificity_threshold_for_paramyxovirus_family(
        self, ndv_cfg: PipelineConfig
    ) -> None:
        """Off-target identity threshold must be tight for avian paramyxovirus family.

        APMV-2 (Yucaipa virus) and APMV-3 (turkey paramyxovirus) are the
        closest non-target paramyxoviruses and share conserved M-gene motifs.
        A threshold >= 0.83 is needed to catch potential cross-reactivity from
        these relatives while remaining strict enough to require actual NDV
        sequence for a positive call.
        """
        assert ndv_cfg.min_identity_threshold >= 0.83, (
            "APMV-2/3 share paramyxovirus M-gene structural motifs with NDV; "
            "identity threshold must be >= 0.83 to flag potential cross-reactivity"
        )

    def test_amplicon_size_within_rt_lamp_window(self, ndv_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= ndv_cfg.f2_b2_min < ndv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, ndv_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking the F2-B2 span."""
        assert ndv_cfg.min_region_length >= ndv_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_ndv(self, ndv_cfg: PipelineConfig) -> None:
        assert "ndv" in str(ndv_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, ndv_cfg: PipelineConfig) -> None:
        assert 0.70 <= ndv_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, ndv_cfg: PipelineConfig) -> None:
        assert 0.70 <= ndv_cfg.min_coverage_threshold <= 0.95
