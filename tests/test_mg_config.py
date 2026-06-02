"""Validate the Mycoplasma gallisepticum mgc2 poultry biosecurity config (config/mg_mgc2.yaml).

Semantic checks specific to Mycoplasma gallisepticum detection via the mgc2 gene
(cytadherence-related protein MGC2 -- the leading published LAMP target for this
pathogen, first validated by Goto et al. 2003):

- MG taxon anchor (NCBI taxon 2097)
- mgc2 entropy threshold appropriate for a conserved single-copy gene
- GC window correctly lowered for the extremely AT-rich MG genome (~30.9% GC)
- Standard LAMP Tm window, not RT-LAMP (mgc2 is a chromosomal DNA target)
- Tight off-target thresholds to exclude avian Mycoplasma cross-reactivity
  (especially Mycoplasma synoviae, the #1 clinical confounder)

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_MG_CONFIG = CONFIG_DIR / "mg_mgc2.yaml"

_RT_LAMP_FLOOR = 63.0  # degC; one-step RT-LAMP co-activity floor
_DNA_LAMP_FLOOR = 60.0  # degC; standard DNA-LAMP Bst 2.0 floor
_MG_GC = 30.9  # approximate genomic GC content of M. gallisepticum


@pytest.fixture
def mg_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the MG mgc2 config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_MG_CONFIG)


def test_mg_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_MG_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestMgMgc2Config:
    """Semantic checks for the Mycoplasma gallisepticum mgc2 poultry biosecurity config.

    M. gallisepticum is the leading mycoplasmal pathogen of commercial poultry
    worldwide, causing chronic respiratory disease in chickens and infectious
    sinusitis in turkeys.  Its low genomic GC content (~30.9%) is the primary
    primer-design challenge; the mgc2 cytadherence gene is the published LAMP
    target of choice (Goto et al. 2003).  This is a DNA target -- no reverse
    transcription is required.
    """

    def test_target_name(self, mg_cfg: PipelineConfig) -> None:
        assert mg_cfg.target_name == "mg_mgc2"

    def test_gene_is_mgc2(self, mg_cfg: PipelineConfig) -> None:
        """mgc2 (cytadherence-related protein 2) is the published MG LAMP marker."""
        assert mg_cfg.gene == "mgc2"

    def test_taxon_is_mycoplasma_gallisepticum(self, mg_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Mycoplasma gallisepticum (NCBI TaxID 2097).

        TaxID 2097 retrieves sequences from the major MG genetic clusters
        (Rlow reference, S6, F-strain vaccine, AP3AS field strains) needed
        for a conserved-region analysis that reflects clinical diversity.
        """
        assert mg_cfg.taxon_id == 2097

    def test_not_rt_lamp(self, mg_cfg: PipelineConfig) -> None:
        """M. gallisepticum is a DNA pathogen; tm_min must NOT be raised to RT-LAMP floor.

        Using the 63 degC RT-LAMP Tm floor (required for AIV, NDV, IBDV, IBV RNA
        targets in the poultry panel) is incorrect for a chromosomal DNA target.
        Raising it unnecessarily would bias primer design toward higher-Tm
        sequences, further reducing the already-scarce primer candidates in
        this extremely AT-rich genome.
        """
        assert mg_cfg.tm_min < _RT_LAMP_FLOOR, (
            "M. gallisepticum mgc2 is a chromosomal DNA target: tm_min should be at "
            f"the standard DNA-LAMP {_DNA_LAMP_FLOOR} degC floor, not the "
            f"RT-LAMP {_RT_LAMP_FLOOR} degC floor"
        )

    def test_standard_dna_lamp_tm_window(self, mg_cfg: PipelineConfig) -> None:
        """Standard DNA-LAMP Tm window: 60-65 degC."""
        assert mg_cfg.tm_min >= _DNA_LAMP_FLOOR
        assert mg_cfg.tm_max <= 65.0

    def test_gc_min_low_enough_for_at_rich_genome(self, mg_cfg: PipelineConfig) -> None:
        """gc_min must be well below 40% to find primers in MG (~30.9% GC).

        M. gallisepticum has one of the lowest genomic GC contents of any avian
        pathogen (~30.9%).  The mgc2 gene has a similar AT-rich codon composition.
        Standard LAMP-Forge configs use gc_min >= 35-40%, which would yield zero
        primer candidates for this organism.  gc_min <= 30% is required to find
        primers in the conserved attachment-domain windows of mgc2.
        """
        assert mg_cfg.gc_min <= 30.0, (
            f"M. gallisepticum genomic GC is ~{_MG_GC}%; gc_min must be <= 30% to find "
            "primer candidates in the AT-rich mgc2 CDS (standard 35-40% floor fails here)"
        )

    def test_gc_max_caps_below_off_target_gc(self, mg_cfg: PipelineConfig) -> None:
        """gc_max must cap primers below the GC-richer avian Mycoplasma off-targets.

        Avian respiratory flora and M. synoviae (~30.6% GC) occupy the same
        oropharyngeal/tracheal niche.  A gc_max < 50% reduces the risk of
        designing primers that also match slightly GC-richer commensal organisms.
        """
        assert mg_cfg.gc_max <= 50.0, (
            "gc_max must be <= 50% to cap primers below GC-richer avian respiratory "
            "off-targets (avian respiratory flora typically 40-65% GC)"
        )

    def test_gc_window_is_wide_enough_for_at_rich_target(self, mg_cfg: PipelineConfig) -> None:
        """GC window must be wide enough to find primer candidates in an AT-rich gene."""
        assert mg_cfg.gc_max - mg_cfg.gc_min >= 15.0, (
            "AT-rich targets need a GC window of at least 15 units; "
            "MG mgc2 spans ~25-45% GC across the CDS"
        )

    def test_entropy_threshold_reflects_mgc2_conservation(self, mg_cfg: PipelineConfig) -> None:
        """mgc2 is under functional constraint; entropy threshold must be in housekeeping range.

        Within-MG mgc2 identity is > 95% across all sequenced strains (Rlow,
        S6, F-strain, AP3AS).  The threshold must be tight enough to exploit the
        invariant cytadherence-domain blocks but not so strict that synonymous
        drift at codon wobble sites is excluded.
        """
        assert mg_cfg.entropy_threshold >= 0.10, (
            "mgc2 has some synonymous inter-strain drift; "
            "entropy_threshold < 0.10 is too strict for a multi-strain analysis"
        )
        assert mg_cfg.entropy_threshold <= 0.30, (
            "mgc2 cytadherence residues are near-invariant within MG; "
            "entropy_threshold > 0.30 admits poorly-conserved positions"
        )

    def test_tight_specificity_thresholds_for_avian_mycoplasma(
        self, mg_cfg: PipelineConfig
    ) -> None:
        """Off-target thresholds must be tight to exclude avian Mycoplasma cross-reactivity.

        Mycoplasma synoviae (MS) shares the same avian respiratory niche and has
        similar genomic GC (~30.6%).  The 0.85 / 0.85 threshold is required to
        surface residual mgc2 cross-reactive primer candidates in silico before
        ordering.
        """
        assert mg_cfg.min_identity_threshold >= 0.85, (
            "Avian Mycoplasma spp. require >= 0.85 identity threshold "
            "to detect cross-reactive primer candidates in silico"
        )
        assert mg_cfg.min_coverage_threshold >= 0.85, (
            "Coverage threshold must be >= 0.85 for avian Mycoplasma specificity checking"
        )

    def test_max_sequences_spans_strain_diversity(self, mg_cfg: PipelineConfig) -> None:
        """Need enough sequences to span key MG genetic clusters for conserved-region analysis."""
        assert mg_cfg.max_sequences >= 20

    def test_amplicon_size_within_lamp_window(self, mg_cfg: PipelineConfig) -> None:
        assert 80 <= mg_cfg.f2_b2_min < mg_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, mg_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert mg_cfg.min_region_length >= mg_cfg.f2_b2_max + 40

    def test_output_dir_references_mg(self, mg_cfg: PipelineConfig) -> None:
        assert "mg" in str(mg_cfg.output_dir).lower()

    def test_identity_threshold_is_in_valid_range(self, mg_cfg: PipelineConfig) -> None:
        assert 0.70 <= mg_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_in_valid_range(self, mg_cfg: PipelineConfig) -> None:
        assert 0.70 <= mg_cfg.min_coverage_threshold <= 0.95

    def test_source_specified(self, mg_cfg: PipelineConfig) -> None:
        """Config must have a taxon_id or non-empty accessions list."""
        assert mg_cfg.taxon_id is not None or len(mg_cfg.accessions) > 0
