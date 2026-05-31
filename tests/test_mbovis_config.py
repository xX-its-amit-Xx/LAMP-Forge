"""Validate the Mycoplasma bovis uvrC BRD supplemental config (config/mbovis_uvrC.yaml).

Semantic checks specific to Mycoplasma bovis detection via the uvrC gene
(DNA repair subunit C -- the leading published LAMP target for this pathogen):

- M. bovis taxon anchor (NCBI taxon 28903)
- uvrC entropy threshold appropriate for a functionally constrained single-species gene
- GC window correctly lowered for the extremely AT-rich M. bovis genome (~29.5% GC)
- Standard LAMP Tm window, not RT-LAMP (uvrC is a chromosomal DNA target)
- Tight off-target thresholds to exclude cross-reactive ruminant Mycoplasma spp.

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_MBOVIS_CONFIG = CONFIG_DIR / "mbovis_uvrC.yaml"

_RT_LAMP_FLOOR = 63.0  # degC; one-step RT-LAMP co-activity floor
_DNA_LAMP_FLOOR = 60.0  # degC; standard DNA-LAMP Bst 2.0 floor
_MBOVIS_GC = 29.5  # approximate genomic GC content of M. bovis


@pytest.fixture
def mbovis_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the M. bovis uvrC config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_MBOVIS_CONFIG)


def test_mbovis_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_MBOVIS_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestMbovIsUvrCConfig:
    """Semantic checks for the Mycoplasma bovis uvrC BRD supplemental config.

    M. bovis is the leading mycoplasmal cause of bovine respiratory disease,
    mastitis, and arthritis.  Its extremely low genomic GC content (~29.5%)
    is the primary design challenge; the uvrC DNA repair gene is the published
    LAMP target of choice (Chen et al. 2017).  This is a DNA target -- no
    reverse transcription is required.
    """

    def test_target_name(self, mbovis_cfg: PipelineConfig) -> None:
        assert mbovis_cfg.target_name == "mbovis_uvrC"

    def test_gene_is_uvrc(self, mbovis_cfg: PipelineConfig) -> None:
        """uvrC (DNA repair subunit C) is the published M. bovis LAMP marker."""
        assert mbovis_cfg.gene == "uvrC"

    def test_taxon_is_mycoplasma_bovis(self, mbovis_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Mycoplasma bovis (NCBI TaxID 28903).

        TaxID 28903 retrieves sequences from diverse M. bovis field strains
        (PG45 reference, HB0801, NP151, Ningxia-1, 08M) needed for a
        conserved-region analysis that reflects clinical diversity.
        """
        assert mbovis_cfg.taxon_id == 28903

    def test_not_rt_lamp(self, mbovis_cfg: PipelineConfig) -> None:
        """M. bovis is a DNA pathogen; tm_min must NOT be raised to RT-LAMP floor.

        Using the 63 degC RT-LAMP Tm floor (required for BRSV, BCoV, BVDV RNA
        targets in the BRD panel) is incorrect for a chromosomal DNA target.
        Raising it unnecessarily would bias primer design toward higher-Tm
        sequences, further reducing the already-scarce primer candidates in
        this extremely AT-rich genome.
        """
        assert mbovis_cfg.tm_min < _RT_LAMP_FLOOR, (
            "M. bovis uvrC is a chromosomal DNA target: tm_min should be at the "
            f"standard DNA-LAMP {_DNA_LAMP_FLOOR} degC floor, not the "
            f"RT-LAMP {_RT_LAMP_FLOOR} degC floor"
        )

    def test_standard_dna_lamp_tm_window(self, mbovis_cfg: PipelineConfig) -> None:
        """Standard DNA-LAMP Tm window: 60-65 degC."""
        assert mbovis_cfg.tm_min >= _DNA_LAMP_FLOOR
        assert mbovis_cfg.tm_max <= 65.0

    def test_gc_min_low_enough_for_at_rich_genome(self, mbovis_cfg: PipelineConfig) -> None:
        """gc_min must be well below 40% to find primers in M. bovis (~29.5% GC).

        M. bovis is one of the most AT-rich livestock pathogens (genomic GC
        ~29.5%).  The uvrC gene itself has a similar AT-rich codon composition.
        Standard LAMP-Forge configs use gc_min >= 35-40%, which would yield
        zero primer candidates for this organism.  gc_min <= 28% is required
        to find primers in the conserved catalytic-domain windows of uvrC.
        """
        assert mbovis_cfg.gc_min <= 28.0, (
            f"M. bovis genomic GC is ~{_MBOVIS_GC}%; gc_min must be <= 28% to find primer "
            "candidates in the AT-rich uvrC CDS (standard 35-40% floor fails here)"
        )

    def test_gc_max_caps_below_off_target_gc(self, mbovis_cfg: PipelineConfig) -> None:
        """gc_max must cap primers below the GC-richer ruminant Mycoplasma off-targets.

        Mycoplasma agalactiae (~33% GC), M. bovirhinis (~31% GC), and M. dispar
        (~28% GC) all occupy the same bovine respiratory niche.  A gc_max < 50%
        reduces the risk of designing primers that also match off-target
        Mycoplasma species with slightly higher GC.
        """
        assert mbovis_cfg.gc_max <= 50.0, (
            "gc_max must be <= 50% to cap primers below the GC-richer ruminant "
            "Mycoplasma off-targets (M. agalactiae ~33%, M. bovirhinis ~31%)"
        )

    def test_gc_window_is_wide_enough_for_at_rich_target(self, mbovis_cfg: PipelineConfig) -> None:
        """GC window must be wide enough to find primer candidates in an AT-rich gene."""
        assert mbovis_cfg.gc_max - mbovis_cfg.gc_min >= 15.0, (
            "AT-rich targets need a GC window of at least 15 units; "
            "M. bovis uvrC spans ~25-45% GC across the CDS"
        )

    def test_entropy_threshold_reflects_uvrc_conservation(self, mbovis_cfg: PipelineConfig) -> None:
        """uvrC is under functional constraint; entropy threshold must be in housekeeping range.

        Within-M. bovis uvrC identity is > 97% across all sequenced strains.
        The threshold must be tight enough to exploit the invariant ATP-binding
        and DNA-binding domain blocks (< 0.25 bits) but not so strict that the
        few synonymous-drift positions at codon wobble sites are excluded
        (>= 0.10 bits is the practical lower bound for any multi-strain analysis).
        """
        assert mbovis_cfg.entropy_threshold >= 0.10, (
            "uvrC has some synonymous inter-strain drift; "
            "entropy_threshold < 0.10 is too strict for a multi-strain analysis"
        )
        assert mbovis_cfg.entropy_threshold <= 0.25, (
            "uvrC catalytic residues are near-invariant within M. bovis; "
            "entropy_threshold > 0.25 admits poorly-conserved positions"
        )

    def test_tight_specificity_thresholds_for_ruminant_mycoplasma(
        self, mbovis_cfg: PipelineConfig
    ) -> None:
        """Off-target thresholds must be tight to exclude ruminant Mycoplasma cross-reactivity.

        Mycoplasma agalactiae shares > 80% 16S rRNA identity with M. bovis.
        The 0.85 / 0.85 threshold is required to surface residual uvrC
        cross-reactive primer candidates in silico before ordering.
        """
        assert mbovis_cfg.min_identity_threshold >= 0.85, (
            "Ruminant Mycoplasma spp. require >= 0.85 identity threshold "
            "to detect cross-reactive primer candidates in silico"
        )
        assert mbovis_cfg.min_coverage_threshold >= 0.85, (
            "Coverage threshold must be >= 0.85 for ruminant Mycoplasma specificity checking"
        )

    def test_max_sequences_spans_strain_diversity(self, mbovis_cfg: PipelineConfig) -> None:
        """Need enough sequences to span key M. bovis strains for conserved-region analysis."""
        assert mbovis_cfg.max_sequences >= 20

    def test_amplicon_size_within_lamp_window(self, mbovis_cfg: PipelineConfig) -> None:
        assert 80 <= mbovis_cfg.f2_b2_min < mbovis_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, mbovis_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert mbovis_cfg.min_region_length >= mbovis_cfg.f2_b2_max + 40

    def test_output_dir_references_mbovis(self, mbovis_cfg: PipelineConfig) -> None:
        assert "mbovis" in str(mbovis_cfg.output_dir).lower()

    def test_identity_threshold_is_in_valid_range(self, mbovis_cfg: PipelineConfig) -> None:
        assert 0.70 <= mbovis_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_in_valid_range(self, mbovis_cfg: PipelineConfig) -> None:
        assert 0.70 <= mbovis_cfg.min_coverage_threshold <= 0.95

    def test_source_specified(self, mbovis_cfg: PipelineConfig) -> None:
        """Config must have a taxon_id or non-empty accessions list."""
        assert mbovis_cfg.taxon_id is not None or len(mbovis_cfg.accessions) > 0
