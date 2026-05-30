"""Validate the APB fthfs oilfield MIC monitoring config (config/apb_fthfs.yaml).

Semantic checks specific to acid-producing bacteria (APB) homoacetogens via fthfs:
- Moorella genus anchor (NCBI taxon 44417)
- Entropy threshold appropriate for fthfs (mid-range, like dsrB/mcrA)
- GC window covers Moorella and Acetobacterium GC range (~44-55%)
- Standard LAMP Tm window, not RT-LAMP (fthfs is a chromosomal DNA target)

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_APB_CONFIG = CONFIG_DIR / "apb_fthfs.yaml"


@pytest.fixture
def apb_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the APB fthfs config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_APB_CONFIG)


def test_apb_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_APB_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestApbFthfsConfig:
    """Semantic checks for the APB fthfs oilfield MIC monitoring config.

    fthfs (formyltetrahydrofolate synthetase) is the canonical molecular marker
    for homoacetogens -- acid-producing bacteria that reduce CO2 to acetic acid
    via the Wood-Ljungdahl pathway.  In oilfields, APB lower produced water pH
    and donate acetate/H2 to SRB, amplifying sulfidogenesis and MIC.
    """

    def test_target_name(self, apb_cfg: PipelineConfig) -> None:
        assert apb_cfg.target_name == "apb_fthfs_oilfield"

    def test_gene_is_fthfs(self, apb_cfg: PipelineConfig) -> None:
        assert apb_cfg.gene == "fthfs"

    def test_taxon_is_moorella_genus(self, apb_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Moorella genus (NCBI taxon 44417).

        Moorella is the dominant thermophilic homoacetogen in high-temperature
        oilfield reservoirs.  Supplementary Acetobacterium accessions in the
        config extend coverage to mesophilic surface-facility environments.
        """
        assert apb_cfg.taxon_id == 44417

    def test_entropy_threshold_reflects_functional_constraint(
        self, apb_cfg: PipelineConfig
    ) -> None:
        """fthfs ATP-grasp motifs are near-invariant; threshold in mid-range.

        fthfs is under catalytic constraint comparable to dsrB and mcrA: the
        ATP-grasp and GHMP-kinase functional domains are highly conserved across
        homoacetogen genera, but synonymous divergence is substantial.  The
        entropy threshold must be low enough to exploit the conserved motif
        blocks (< 0.40 bits) but not so strict as to yield no windows across
        the Moorella-Acetobacterium-Sporomusa phylogenetic range (>= 0.20 bits).
        """
        assert apb_cfg.entropy_threshold >= 0.20, (
            "fthfs diverges across homoacetogen genera; entropy_threshold < 0.20 "
            "is too strict and will yield no usable conserved windows"
        )
        assert apb_cfg.entropy_threshold <= 0.40, (
            "fthfs ATP-grasp and GHMP-kinase motifs are near-invariant; "
            "entropy_threshold > 0.40 admits poorly-conserved positions"
        )

    def test_gc_range_covers_homoacetogens(self, apb_cfg: PipelineConfig) -> None:
        """GC window must span Moorella (~55%) and Acetobacterium (~44%) GC content.

        Moorella thermoacetica has ~55% GC; Acetobacterium woodii has ~44% GC.
        The design window must cover both genera to produce usable primers when
        accessions from both are included in the sequence set.
        """
        assert apb_cfg.gc_min <= 44.0, (
            "gc_min must be <= 44% to accommodate Acetobacterium woodii (~44% GC)"
        )
        assert apb_cfg.gc_max >= 55.0, (
            "gc_max must be >= 55% to accommodate Moorella thermoacetica (~55% GC)"
        )

    def test_gc_window_is_sufficiently_wide(self, apb_cfg: PipelineConfig) -> None:
        """The gc window must be wide enough to span inter-genus GC diversity."""
        assert apb_cfg.gc_max - apb_cfg.gc_min >= 15.0, (
            "APB homoacetogens span ~44-55% GC; gc window must be at least 15 units wide"
        )

    def test_max_sequences_spans_genus_diversity(self, apb_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover key Moorella species and supplementary genera."""
        assert apb_cfg.max_sequences >= 20

    def test_standard_lamp_tm_window(self, apb_cfg: PipelineConfig) -> None:
        """fthfs is a chromosomal DNA target; standard 60-65 degC LAMP window is correct."""
        assert apb_cfg.tm_min >= 60.0
        assert apb_cfg.tm_max <= 65.0

    def test_not_rt_lamp(self, apb_cfg: PipelineConfig) -> None:
        """fthfs is a chromosomal DNA target; tm_min must not be raised to RT-LAMP floor.

        The RT-LAMP one-step 63 degC floor (used for PRRSV, FMDV, AIV RNA targets)
        is not needed here.  Using it unnecessarily would bias primer design toward
        higher-melting sequences, reducing the primer candidate pool.
        """
        assert apb_cfg.tm_min < 63.5, (
            "APB / fthfs is a DNA target: tm_min should be at the standard LAMP 60 degC, "
            "not the RT-LAMP one-step co-activity floor of 63 degC"
        )

    def test_amplicon_size_within_lamp_window(self, apb_cfg: PipelineConfig) -> None:
        assert 80 <= apb_cfg.f2_b2_min < apb_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, apb_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert apb_cfg.min_region_length >= apb_cfg.f2_b2_max + 40

    def test_output_dir_references_apb(self, apb_cfg: PipelineConfig) -> None:
        assert "apb" in str(apb_cfg.output_dir).lower()

    def test_identity_threshold_is_standard(self, apb_cfg: PipelineConfig) -> None:
        """Off-target identity threshold must be in the recommended range."""
        assert 0.70 <= apb_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_standard(self, apb_cfg: PipelineConfig) -> None:
        """Off-target coverage threshold must be in the recommended range."""
        assert 0.70 <= apb_cfg.min_coverage_threshold <= 0.95
