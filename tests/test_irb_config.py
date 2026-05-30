"""Validate the IRB omcA oilfield MIC monitoring config (config/irb_omcA.yaml).

Semantic checks specific to iron-reducing bacteria (IRB) via omcA:
- Shewanella genus anchor (NCBI taxon 22)
- Entropy threshold appropriate for omcA (between housekeeping and dsrB/mcrA)
- GC window covers Shewanella GC range (45-52%)
- Standard LAMP Tm window, not RT-LAMP (omcA is a chromosomal DNA target)

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_IRB_CONFIG = CONFIG_DIR / "irb_omcA.yaml"


@pytest.fixture
def irb_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the IRB omcA config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_IRB_CONFIG)


def test_irb_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_IRB_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestIrbOmcAConfig:
    """Semantic checks for the IRB omcA oilfield MIC monitoring config.

    omcA encodes the outer membrane decaheme cytochrome A in Shewanella spp.,
    the primary functional marker for iron-reducing bacteria (IRB) in oilfield
    produced water.  It sits between tight housekeeping-gene conservation
    (rpoB, invA) and loose polyphyletic functional-gene conservation (dsrB,
    mcrA), so its entropy threshold and GC window differ from both extremes.
    """

    def test_target_name(self, irb_cfg: PipelineConfig) -> None:
        assert irb_cfg.target_name == "irb_omcA_oilfield"

    def test_gene_is_omca(self, irb_cfg: PipelineConfig) -> None:
        assert irb_cfg.gene == "omcA"

    def test_taxon_is_shewanella_genus(self, irb_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Shewanella genus (NCBI taxon 22)."""
        assert irb_cfg.taxon_id == 22

    def test_entropy_threshold_reflects_functional_constraint(
        self, irb_cfg: PipelineConfig
    ) -> None:
        """omcA heme motifs are invariant; entropy threshold should be < dsrB/mcrA."""
        assert irb_cfg.entropy_threshold <= 0.35, (
            "omcA is under stronger catalytic constraint than dsrB/mcrA; "
            "entropy_threshold should be <= 0.35 to exploit conserved heme blocks"
        )
        assert irb_cfg.entropy_threshold >= 0.20, (
            "omcA diverges across Shewanella species; entropy_threshold < 0.20 "
            "is too strict and will yield no usable conserved windows"
        )

    def test_gc_range_covers_shewanella(self, irb_cfg: PipelineConfig) -> None:
        """Shewanella GC content is 45-52%; primer window must contain this range."""
        assert irb_cfg.gc_min <= 45.0, (
            "gc_min must be <= 45% to reach the low end of Shewanella GC range"
        )
        assert irb_cfg.gc_max >= 52.0, (
            "gc_max must be >= 52% to reach the high end of Shewanella GC range"
        )

    def test_max_sequences_spans_shewanella_diversity(self, irb_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover key oilfield Shewanella species."""
        assert irb_cfg.max_sequences >= 20

    def test_standard_lamp_tm_window(self, irb_cfg: PipelineConfig) -> None:
        """IRB targets are DNA; standard 60-65 degC LAMP window is correct (not RT-LAMP)."""
        assert irb_cfg.tm_min >= 60.0
        assert irb_cfg.tm_max <= 65.0

    def test_not_rt_lamp(self, irb_cfg: PipelineConfig) -> None:
        """omcA is a chromosomal DNA target; tm_min must not be raised to RT-LAMP floor."""
        assert irb_cfg.tm_min < 63.5, (
            "IRB / omcA is a DNA target: tm_min should be at the standard LAMP 60 degC, "
            "not the RT-LAMP one-step co-activity floor of 63 degC"
        )

    def test_amplicon_size_within_lamp_window(self, irb_cfg: PipelineConfig) -> None:
        assert 80 <= irb_cfg.f2_b2_min < irb_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, irb_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outside the F2-B2 span."""
        assert irb_cfg.min_region_length >= irb_cfg.f2_b2_max + 40

    def test_output_dir_references_irb(self, irb_cfg: PipelineConfig) -> None:
        assert "irb" in str(irb_cfg.output_dir).lower()
