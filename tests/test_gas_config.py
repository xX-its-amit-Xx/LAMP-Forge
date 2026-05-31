"""Validate the GAS speB point-of-care config (config/gas_speB.yaml).

Semantic checks specific to Group A Streptococcus (Streptococcus pyogenes)
detection via the speB gene (streptococcal cysteine protease exotoxin B):

- S. pyogenes taxon anchor (NCBI taxon 1314)
- speB entropy threshold appropriate for a well-conserved housekeeping-like gene
- GC window set for AT-rich S. pyogenes genome (~38.5% GC)
- Standard LAMP Tm window, not RT-LAMP (speB is a chromosomal DNA target)
- Tight off-target thresholds to exclude intra-genus streptococcal cross-reactivity

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_GAS_CONFIG = CONFIG_DIR / "gas_speB.yaml"


@pytest.fixture
def gas_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the GAS speB config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_GAS_CONFIG)


def test_gas_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_GAS_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestGasSpeBAConfig:
    """Semantic checks for the GAS speB point-of-care pharyngitis config.

    speB (streptococcal cysteine protease exotoxin B) is a GAS-specific gene
    present in virtually all Streptococcus pyogenes strains and absent from
    other streptococcal species.  It is the basis of published LAMP assays for
    rapid GAS pharyngitis detection at the point of care.  Detection drives the
    antibiotic-prescribing decision that prevents rheumatic fever and
    post-streptococcal glomerulonephritis.
    """

    def test_target_name(self, gas_cfg: PipelineConfig) -> None:
        assert gas_cfg.target_name == "gas_speB_poc"

    def test_gene_is_speb(self, gas_cfg: PipelineConfig) -> None:
        assert gas_cfg.gene == "speB"

    def test_taxon_is_streptococcus_pyogenes(self, gas_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Streptococcus pyogenes (NCBI taxon 1314).

        S. pyogenes is the causative agent of GAS pharyngitis, impetigo, and
        scarlet fever.  Taxon 1314 retrieves sequences across the globally
        dominant M-types (M1, M3, M12, M28, M89) needed for conserved-region
        analysis that spans clinical diversity.
        """
        assert gas_cfg.taxon_id == 1314

    def test_entropy_threshold_reflects_speb_conservation(self, gas_cfg: PipelineConfig) -> None:
        """speB is well-conserved within GAS; threshold must be in the housekeeping-gene range.

        speB is under strong pathogenicity-island selection pressure within
        S. pyogenes: the catalytic Cys/His dyad and prodomain cleavage site are
        near-invariant across M-types.  The entropy threshold must be tight
        enough to exploit these conserved blocks (< 0.30 bits) but not so strict
        as to yield no windows when inter-M-type synonymous drift is present
        (>= 0.10 bits).
        """
        assert gas_cfg.entropy_threshold >= 0.10, (
            "speB has some synonymous inter-M-type drift; entropy_threshold < 0.10 is too strict"
        )
        assert gas_cfg.entropy_threshold <= 0.30, (
            "speB catalytic residues are near-invariant; "
            "entropy_threshold > 0.30 admits poorly-conserved positions"
        )

    def test_gc_range_covers_at_rich_gas(self, gas_cfg: PipelineConfig) -> None:
        """GC window must accommodate S. pyogenes' AT-rich genome (~38.5% GC).

        S. pyogenes has one of the lowest genomic GC contents among beta-
        haemolytic streptococci (~38.5%).  The speB gene itself is ~37% GC.
        The gc_min must be low enough to allow primer candidates in the target
        (gc_min <= 35.0) while gc_max caps primers well below the GC-rich
        intra-genus off-targets (S. pneumoniae ~39.7%, but some streptococcal
        off-targets reach ~45%).
        """
        assert gas_cfg.gc_min <= 35.0, (
            "gc_min must be <= 35% to accommodate AT-rich speB codons "
            "in S. pyogenes (~38.5% genomic GC)"
        )
        assert gas_cfg.gc_max >= 50.0, (
            "gc_max must be >= 50% to allow primers in more GC-rich windows of the ~1.4 kb speB CDS"
        )

    def test_gc_window_is_wide_enough_for_at_rich_target(self, gas_cfg: PipelineConfig) -> None:
        """GC window must be wide enough to find primer candidates in an AT-rich gene."""
        assert gas_cfg.gc_max - gas_cfg.gc_min >= 15.0, (
            "AT-rich targets need a GC window of at least 15 units; "
            "speB spans ~30-50% GC across the CDS"
        )

    def test_max_sequences_spans_m_type_diversity(self, gas_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover key M-types for conserved-region analysis."""
        assert gas_cfg.max_sequences >= 20

    def test_standard_lamp_tm_window(self, gas_cfg: PipelineConfig) -> None:
        """speB is a chromosomal DNA target; standard 60-65 degC LAMP window is correct."""
        assert gas_cfg.tm_min >= 60.0
        assert gas_cfg.tm_max <= 65.0

    def test_not_rt_lamp(self, gas_cfg: PipelineConfig) -> None:
        """speB is a chromosomal DNA target; tm_min must not be raised to RT-LAMP floor.

        The RT-LAMP one-step 63 degC floor (used for PRRSV, FMDV, AIV, NDV RNA
        targets) is not needed here.  Raising it unnecessarily would bias primer
        design toward higher-Tm sequences, reducing candidates in this AT-rich gene.
        """
        assert gas_cfg.tm_min < 63.5, (
            "GAS / speB is a DNA target: tm_min should be at the standard "
            "LAMP 60 degC floor, not the RT-LAMP one-step co-activity floor"
        )

    def test_tight_specificity_thresholds(self, gas_cfg: PipelineConfig) -> None:
        """Off-target thresholds must be tight to exclude intra-genus streptococci.

        S. agalactiae, S. pneumoniae, and other streptococci share GC content
        and some conserved housekeeping gene blocks with S. pyogenes.  The
        0.85 x 0.85 threshold is required to catch false-positive primer
        matches before wet-lab validation.
        """
        assert gas_cfg.min_identity_threshold >= 0.85, (
            "Intra-genus streptococci require >= 0.85 identity threshold "
            "to detect cross-reactive primer candidates in silico"
        )
        assert gas_cfg.min_coverage_threshold >= 0.85, (
            "Coverage threshold must be >= 0.85 for intra-genus specificity checking"
        )

    def test_amplicon_size_within_lamp_window(self, gas_cfg: PipelineConfig) -> None:
        assert 80 <= gas_cfg.f2_b2_min < gas_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, gas_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert gas_cfg.min_region_length >= gas_cfg.f2_b2_max + 40

    def test_output_dir_references_gas(self, gas_cfg: PipelineConfig) -> None:
        assert "gas" in str(gas_cfg.output_dir).lower()

    def test_identity_threshold_is_in_recommended_range(self, gas_cfg: PipelineConfig) -> None:
        assert 0.70 <= gas_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_in_recommended_range(self, gas_cfg: PipelineConfig) -> None:
        assert 0.70 <= gas_cfg.min_coverage_threshold <= 0.95
