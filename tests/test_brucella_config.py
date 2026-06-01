"""Validate the Brucella IS711 farm/POC config (config/brucella_is711.yaml).

Semantic checks specific to Brucella genus detection via IS711:

- Brucella genus taxon anchor (NCBI taxon 234) for pan-species coverage
- IS711 entropy threshold appropriate for a conserved insertion sequence
- GC window set for GC-rich Brucella genome (~57% GC)
- Standard LAMP Tm window, not RT-LAMP (IS711 is a chromosomal DNA target)
- Tight off-target thresholds to exclude Ochrobactrum and other alpha-proteobacteria
- Multi-copy nature of IS711 (5-40 copies per genome) makes it ideal for
  high-sensitivity LAMP detection of WOAH-notifiable brucellosis

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_BRUCELLA_CONFIG = CONFIG_DIR / "brucella_is711.yaml"


@pytest.fixture
def brucella_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the Brucella IS711 config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_BRUCELLA_CONFIG)


def test_brucella_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_BRUCELLA_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestBrucellaIS711Config:
    """Semantic checks for the Brucella IS711 farm/point-of-care config.

    IS711 (also annotated IS6501/ISBm2) is a Brucella-specific insertion sequence
    present in 5-40 copies per genome across all Brucella species.  Its multi-copy
    nature gives exceptional LAMP sensitivity (10-100x lower LOD than single-copy
    targets), making it the gold-standard target for brucellosis diagnostics
    (Abd El Wahed et al. 2015 PLoS NTD).

    Brucella abortus (bovine), B. melitensis (caprine/ovine), and B. suis (porcine)
    are WOAH-notifiable pathogens that affect BioVind's farm-animal biosecurity
    vertical; B. melitensis and B. abortus also cause human brucellosis, covering
    the point-of-care vertical.
    """

    def test_target_name(self, brucella_cfg: PipelineConfig) -> None:
        assert brucella_cfg.target_name == "brucella_is711"

    def test_gene_is_is711(self, brucella_cfg: PipelineConfig) -> None:
        assert brucella_cfg.gene == "IS711"

    def test_taxon_is_brucella_genus(self, brucella_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Brucella genus (NCBI taxon 234).

        Taxon 234 (Brucella) retrieves sequences spanning all nine species and their
        biovars.  This is required for pan-species IS711 detection: B. abortus
        biovars 1-7/9, B. melitensis biovars 1-3, B. suis biovars 1-5, B. ovis,
        B. canis, and the newly described marine species share the IS711 element
        with high sequence identity (>98%) across a conserved internal region.
        """
        assert brucella_cfg.taxon_id == 234

    def test_entropy_threshold_reflects_is711_conservation(
        self, brucella_cfg: PipelineConfig
    ) -> None:
        """IS711 is highly conserved within Brucella; threshold must be in a tight range.

        The IS711 transposase open reading frame is under strong purifying selection
        across all Brucella species.  The internal conserved region used in published
        LAMP assays shows <2% nucleotide divergence across species, justifying a tight
        entropy threshold.  However, the threshold must not be so strict that it
        excludes windows when minor insertion polymorphisms are present.
        """
        assert brucella_cfg.entropy_threshold >= 0.10, (
            "IS711 has minor inter-species drift; entropy_threshold < 0.10 is too strict"
        )
        assert brucella_cfg.entropy_threshold <= 0.30, (
            "IS711 transposase region is near-invariant across Brucella; "
            "entropy_threshold > 0.30 admits poorly-conserved positions"
        )

    def test_gc_range_covers_gc_rich_brucella(self, brucella_cfg: PipelineConfig) -> None:
        """GC window must accommodate Brucella's GC-rich genome (~57% GC).

        Brucella has one of the highest GC contents among the alpha-proteobacteria
        (~57% genomic GC).  IS711 itself is approximately 55-60% GC.  The gc_min
        must be high enough to be realistic for primers in this GC-rich target
        (gc_min >= 45.0) while gc_max allows the upper GC range of the element.
        """
        assert brucella_cfg.gc_min >= 45.0, (
            "gc_min must be >= 45% to reflect IS711 GC content in Brucella (~57% genomic GC); "
            "lower values would admit primers with unrealistically low GC for this target"
        )
        assert brucella_cfg.gc_max <= 72.0, (
            "gc_max must be <= 72% to cap primers at the upper range of IS711 GC content"
        )

    def test_gc_window_is_wide_enough_for_gc_rich_target(
        self, brucella_cfg: PipelineConfig
    ) -> None:
        """GC window must be wide enough to find primer candidates in a GC-rich gene."""
        assert brucella_cfg.gc_max - brucella_cfg.gc_min >= 15.0, (
            "GC-rich targets still vary locally; a window of at least 15 units is needed "
            "to find all valid primer candidates across IS711"
        )

    def test_max_sequences_spans_species_diversity(self, brucella_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover major Brucella species for conserved-region analysis."""
        assert brucella_cfg.max_sequences >= 20

    def test_standard_lamp_tm_window(self, brucella_cfg: PipelineConfig) -> None:
        """IS711 is a chromosomal DNA target; standard 60-65 degC LAMP window is correct."""
        assert brucella_cfg.tm_min >= 60.0
        assert brucella_cfg.tm_max <= 65.0

    def test_not_rt_lamp(self, brucella_cfg: PipelineConfig) -> None:
        """IS711 is a chromosomal DNA target; tm_min must not be raised to RT-LAMP floor.

        The RT-LAMP one-step 63 degC floor (used for PRRSV, FMDV, AIV, NDV RNA
        targets) is not needed here.  IS711 is a DNA insertion sequence; raising
        tm_min unnecessarily would bias primer design toward higher-Tm sequences
        and reduce candidates in the GC-rich target.
        """
        assert brucella_cfg.tm_min < 63.5, (
            "Brucella / IS711 is a DNA target: tm_min should be at the standard "
            "LAMP 60 degC floor, not the RT-LAMP one-step co-activity floor"
        )

    def test_tight_specificity_thresholds(self, brucella_cfg: PipelineConfig) -> None:
        """Off-target thresholds must be tight to exclude alpha-proteobacterial relatives.

        Ochrobactrum anthropi (now Brucella anthropi) and other Rhizobiaceae share
        significant genomic similarity with Brucella spp.  The 0.85 x 0.85 threshold
        is required to catch false-positive primer matches against these off-targets
        before wet-lab validation.
        """
        assert brucella_cfg.min_identity_threshold >= 0.85, (
            "Alpha-proteobacterial relatives require >= 0.85 identity threshold "
            "to detect cross-reactive primer candidates in silico"
        )
        assert brucella_cfg.min_coverage_threshold >= 0.85, (
            "Coverage threshold must be >= 0.85 for specificity checking "
            "against Ochrobactrum and other Rhizobiaceae off-targets"
        )

    def test_amplicon_size_within_lamp_window(self, brucella_cfg: PipelineConfig) -> None:
        assert 80 <= brucella_cfg.f2_b2_min < brucella_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, brucella_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert brucella_cfg.min_region_length >= brucella_cfg.f2_b2_max + 40

    def test_output_dir_references_brucella(self, brucella_cfg: PipelineConfig) -> None:
        assert "brucella" in str(brucella_cfg.output_dir).lower()

    def test_identity_threshold_is_in_recommended_range(self, brucella_cfg: PipelineConfig) -> None:
        assert 0.70 <= brucella_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_in_recommended_range(self, brucella_cfg: PipelineConfig) -> None:
        assert 0.70 <= brucella_cfg.min_coverage_threshold <= 0.95
