"""Validate the C. diff tcdB hospital point-of-care config (config/cdiff_tcdB.yaml).

Semantic checks specific to Clostridioides difficile detection via the tcdB gene
(toxin B large clostridial glucosyltransferase):

- C. difficile taxon anchor (NCBI taxon 1496)
- tcdB entropy threshold appropriate for a virulence gene with inter-ribotype drift
- GC window set for the extremely AT-rich C. diff genome (~28.6% GC)
- Standard LAMP Tm window, not RT-LAMP (tcdB is a chromosomal DNA target)
- Tight off-target thresholds to exclude Clostridioides sordellii lethal toxin
- GC cap below GC-richer gut-flora off-targets (Bacteroides, E. coli)

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_CDIFF_CONFIG = CONFIG_DIR / "cdiff_tcdB.yaml"


@pytest.fixture
def cdiff_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the C. diff tcdB config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_CDIFF_CONFIG)


def test_cdiff_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_CDIFF_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestCdiffTcdBConfig:
    """Semantic checks for the C. diff tcdB hospital point-of-care config.

    tcdB (toxin B glucosyltransferase) is the IDSA/SHEA-recommended molecular
    target for Clostridioides difficile infection (CDI) diagnosis.  Its presence
    in stool drives the cohort isolation and antibiotic-prescribing decision.
    A portable 30-minute LAMP assay replaces the 1-4 hour laboratory NAAT and
    enables same-visit decisions in urgent-care and hospital settings.
    """

    def test_target_name(self, cdiff_cfg: PipelineConfig) -> None:
        assert cdiff_cfg.target_name == "cdiff_tcdB_poc"

    def test_gene_is_tcdB(self, cdiff_cfg: PipelineConfig) -> None:
        assert cdiff_cfg.gene == "tcdB"

    def test_taxon_is_clostridioides_difficile(self, cdiff_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Clostridioides difficile (NCBI taxon 1496).

        Taxon 1496 retrieves sequences across the clinically dominant ribotypes
        (RT001, RT027/NAP1, RT106, RT078, RT014, RT017) needed for
        conserved-region analysis that spans current epidemic lineages.
        """
        assert cdiff_cfg.taxon_id == 1496

    def test_entropy_threshold_reflects_tcdB_conservation(self, cdiff_cfg: PipelineConfig) -> None:
        """tcdB catalytic core is conserved but with inter-ribotype drift.

        The N-terminal glucosyltransferase domain of tcdB is under strong
        purifying selection at the catalytic DXD motif and substrate-binding
        loops.  However, synonymous drift between major ribotypes (RT001 vs.
        RT027 vs. RT078) is greater than in chromosomal housekeeping genes
        (rpoB, speB).  The entropy threshold should be slightly relaxed
        relative to housekeeping genes (>= 0.20) but not so loose as to admit
        hypervariable inter-ribotype positions (< 0.35 bits).
        """
        assert cdiff_cfg.entropy_threshold >= 0.20, (
            "tcdB catalytic core is more conserved than 0.20 bits; "
            "entropy_threshold < 0.20 would be unusably strict for this gene"
        )
        assert cdiff_cfg.entropy_threshold <= 0.35, (
            "entropy_threshold > 0.35 admits too much inter-ribotype drift "
            "and may generate primers outside the conserved catalytic domain"
        )

    def test_gc_min_accommodates_at_rich_genome(self, cdiff_cfg: PipelineConfig) -> None:
        """gc_min must be low enough to find primers in C. diff's AT-rich genome.

        C. difficile has one of the most AT-rich genomes among Gram-positive
        pathogens (~28.6% GC genome-wide; Sebaihia et al. 2006).  The tcdB
        gene and the surrounding PaLoc are ~26-30% GC.  Standard LAMP primer
        GC floors of 40% would yield no candidate windows; gc_min must be
        set below 30% to allow primer design in this organism.
        """
        assert cdiff_cfg.gc_min <= 30.0, (
            "C. difficile genome is ~28.6% GC (extremely AT-rich); "
            "gc_min must be <= 30% to find primer candidates in the tcdB gene"
        )

    def test_gc_max_stays_below_gut_flora_off_targets(self, cdiff_cfg: PipelineConfig) -> None:
        """gc_max must be capped below GC-richer gut-flora off-targets.

        Major stool background organisms have higher GC content than C. diff:
        Bacteroides fragilis ~43-48% GC, Escherichia coli K-12 ~51% GC.
        Capping gc_max at or below 52% reduces the chance of designing primers
        that cross-amplify these abundant gut commensals.
        """
        assert cdiff_cfg.gc_max <= 52.0, (
            "gc_max should be <= 52% to avoid cross-amplification of "
            "GC-richer gut-flora off-targets (Bacteroides ~48%, E. coli ~51%)"
        )

    def test_gc_window_is_wide_enough_for_at_rich_target(self, cdiff_cfg: PipelineConfig) -> None:
        """GC window must be wide enough to find candidates in a low-GC gene."""
        assert cdiff_cfg.gc_max - cdiff_cfg.gc_min >= 15.0, (
            "AT-rich targets need a GC window of at least 15 units; "
            "tcdB spans ~25-45% GC across the conserved catalytic domain"
        )

    def test_standard_lamp_tm_window(self, cdiff_cfg: PipelineConfig) -> None:
        """tcdB is a chromosomal DNA target; standard 60-65 degC LAMP window is correct."""
        assert cdiff_cfg.tm_min >= 60.0
        assert cdiff_cfg.tm_max <= 65.0

    def test_not_rt_lamp_tm_floor(self, cdiff_cfg: PipelineConfig) -> None:
        """C. difficile tcdB is a chromosomal DNA target; RT-LAMP tm_min floor is wrong here.

        The RT-LAMP one-step 63 degC floor (used for RNA-virus targets such as
        PRRSV, FMDV, AIV, NDV) is not needed for a bacterial chromosomal gene.
        Raising tm_min to the RT-LAMP floor would artificially reduce the
        candidate pool in this AT-rich gene without any benefit.
        """
        assert cdiff_cfg.tm_min < 63.5, (
            "C. diff tcdB is a chromosomal DNA target; tm_min must not be "
            "raised to the RT-LAMP one-step co-activity floor (63 degC)"
        )

    def test_tight_specificity_thresholds(self, cdiff_cfg: PipelineConfig) -> None:
        """Off-target thresholds must be tight to exclude C. sordellii lethal toxin.

        Clostridioides sordellii carries a lethal toxin (LT) with ~40% protein
        identity to tcdB.  Although nucleotide identity in the catalytic domain
        is lower, the 0.85/0.85 threshold is required to catch any residual
        cross-reactive primer candidates before wet-lab validation.
        """
        assert cdiff_cfg.min_identity_threshold >= 0.85, (
            "C. sordellii LT is the chief tcdB off-target; "
            ">= 0.85 identity threshold required to flag potential cross-reactants"
        )
        assert cdiff_cfg.min_coverage_threshold >= 0.85, (
            "Coverage threshold must be >= 0.85 to surface partial primer "
            "matches against C. sordellii lethal toxin homologs"
        )

    def test_max_sequences_spans_ribotype_diversity(self, cdiff_cfg: PipelineConfig) -> None:
        """Need enough sequences to span epidemic ribotypes for conserved-region analysis."""
        assert cdiff_cfg.max_sequences >= 20

    def test_amplicon_size_within_lamp_window(self, cdiff_cfg: PipelineConfig) -> None:
        assert 80 <= cdiff_cfg.f2_b2_min < cdiff_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, cdiff_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert cdiff_cfg.min_region_length >= cdiff_cfg.f2_b2_max + 40

    def test_output_dir_references_cdiff(self, cdiff_cfg: PipelineConfig) -> None:
        assert "cdiff" in str(cdiff_cfg.output_dir).lower()

    def test_identity_threshold_is_in_recommended_range(self, cdiff_cfg: PipelineConfig) -> None:
        assert 0.70 <= cdiff_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_in_recommended_range(self, cdiff_cfg: PipelineConfig) -> None:
        assert 0.70 <= cdiff_cfg.min_coverage_threshold <= 0.95
