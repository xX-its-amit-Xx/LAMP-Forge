"""Validate the universal 16S rRNA sample-adequacy control config (config/general_16s.yaml).

Semantic checks specific to the universal bacterial 16S rRNA control:
- Broadest possible taxonomic anchor (Bacteria, NCBI taxon 2)
- Entropy threshold is unusually strict (<=0.15) for universal primer placement
- max_sequences is high enough for cross-phylum conservation confirmation
- Standard LAMP Tm window, not RT-LAMP (16S rDNA is a chromosomal DNA target)
- Off-target identity threshold is tighter than standard (>=0.83) to catch
  potential archaeal rRNA cross-reactivity early

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_16S_CONFIG = CONFIG_DIR / "general_16s.yaml"


@pytest.fixture
def s16_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the universal 16S rRNA control config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_16S_CONFIG)


def test_16s_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_16S_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestUniversal16SConfig:
    """Semantic checks for the universal 16S rRNA sample-adequacy control config.

    The 16S rRNA gene is used here not as a specific microbial target but as an
    internal positive control: any bacterium in the sample will amplify, confirming
    that DNA extraction succeeded and sufficient bacterial DNA reached the reaction.

    This is the sixth channel in the oilfield MIC/souring panel (alongside SRB
    dsrB, methanogens mcrA, IRB omcA, APB fthfs, and NRB narG).  It is also the
    universal sample-adequacy control for farm biosecurity and point-of-care panels
    that detect bacterial pathogens.
    """

    def test_target_name(self, s16_cfg: PipelineConfig) -> None:
        assert s16_cfg.target_name == "universal_16S_control"

    def test_gene_is_16s(self, s16_cfg: PipelineConfig) -> None:
        """Gene must be set to the 16S ribosomal RNA annotation."""
        assert s16_cfg.gene is not None
        assert "16S" in s16_cfg.gene or "16s" in s16_cfg.gene.lower()

    def test_taxon_anchors_all_bacteria(self, s16_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Bacteria (NCBI taxon 2) for maximum inclusivity.

        All other LAMP-Forge assays anchor to a narrow clade (a genus or species).
        The 16S universal control anchors to the entire Bacteria domain to pull
        sequences from diverse phyla and confirm that the conserved primer sites
        are genuinely universal, not just conserved within a single clade.
        """
        assert s16_cfg.taxon_id == 2, (
            "universal 16S control must anchor to Bacteria (NCBI taxon 2); "
            "a narrower anchor will not confirm cross-phylum universality"
        )

    def test_entropy_threshold_is_strict_for_universality(
        self, s16_cfg: PipelineConfig
    ) -> None:
        """Entropy threshold must be unusually strict for universal primer placement.

        Unlike functional-gene markers (dsrB: 0.35, mcrA: 0.35, narG: 0.35),
        which tolerate synonymous divergence between genera, the 16S universal
        control requires primers at positions that are invariant across ALL bacteria.
        Entropy > 0.15 bits at a primer site means the position varies between
        phyla and will cause false-negatives for some organisms in the sample.

        The invariant blocks of 16S rRNA (Lane et al. 1985; Woese & Fox 1977)
        have entropy < 0.05 bits in large curated databases; 0.15 bits is a
        generous ceiling that still selects universally conserved positions.
        """
        assert s16_cfg.entropy_threshold <= 0.15, (
            "16S universal control needs invariant primer positions: "
            "entropy_threshold must be <= 0.15 bits for cross-phylum universality "
            f"(got {s16_cfg.entropy_threshold}; compare dsrB/narG at 0.35)"
        )

    def test_max_sequences_ensures_cross_phylum_diversity(
        self, s16_cfg: PipelineConfig
    ) -> None:
        """Enough sequences to sample major bacterial phyla.

        A universal 16S assay must be validated against Firmicutes,
        Proteobacteria, Actinobacteria, and Bacteroidetes at minimum (the four
        dominant phyla by 16S-based surveys).  Thirty sequences or more gives
        statistically meaningful diversity representation.
        """
        assert s16_cfg.max_sequences >= 30, (
            "universal 16S control requires >= 30 sequences for cross-phylum "
            "conservation validation across major bacterial phyla "
            f"(got {s16_cfg.max_sequences})"
        )

    def test_standard_lamp_tm_window(self, s16_cfg: PipelineConfig) -> None:
        """16S rDNA is a chromosomal DNA target; standard 60-65 degC window is correct."""
        assert s16_cfg.tm_min >= 60.0
        assert s16_cfg.tm_max <= 65.0

    def test_not_rt_lamp(self, s16_cfg: PipelineConfig) -> None:
        """16S rDNA is a DNA target; tm_min must not be raised to the RT-LAMP floor.

        The RT-LAMP one-step 63 degC floor (used for RNA-virus targets like
        PRRSV, FMDV, AIV, and NDV) is not appropriate for a DNA target.
        Raising tm_min unnecessarily narrows the primer candidate pool without
        any thermodynamic benefit for a DNA-based assay.
        """
        assert s16_cfg.tm_min < 63.5, (
            "16S rDNA is a DNA target: tm_min should be at the standard LAMP "
            "60 degC floor, not the RT-LAMP co-activity floor of 63 degC "
            f"(got {s16_cfg.tm_min})"
        )

    def test_gc_min_accommodates_at_rich_firmicutes(self, s16_cfg: PipelineConfig) -> None:
        """gc_min must accommodate AT-rich Firmicutes in the 16S conserved regions.

        Firmicutes 16S conserved flanking blocks are ~45% GC.  gc_min above 45%
        risks excluding primer candidates in these conserved blocks when Firmicutes
        sequences are included in the design set.
        """
        assert s16_cfg.gc_min <= 42.0, (
            "gc_min must be <= 42% to accommodate AT-rich Firmicutes in the "
            "16S conserved flanking regions (~45% GC)"
        )

    def test_gc_max_accommodates_gc_rich_actinobacteria(
        self, s16_cfg: PipelineConfig
    ) -> None:
        """gc_max must accommodate GC-rich Actinobacteria in the 16S conserved regions.

        Actinobacteria 16S conserved flanking blocks can reach ~57% GC
        (Streptomyces, Mycobacterium).  gc_max below 60% risks excluding primer
        candidates in these GC-richer conserved blocks.
        """
        assert s16_cfg.gc_max >= 60.0, (
            "gc_max must be >= 60% to accommodate GC-rich Actinobacteria in the "
            "16S conserved flanking regions (~57% GC)"
        )

    def test_gc_window_is_moderately_wide(self, s16_cfg: PipelineConfig) -> None:
        """GC window must span cross-phylum GC variation in the 16S conserved regions."""
        assert s16_cfg.gc_max - s16_cfg.gc_min >= 15.0, (
            "16S rRNA conserved regions span ~15 GC percentage points across "
            "major bacterial phyla; gc window must be at least 15 units wide"
        )

    def test_amplicon_size_within_lamp_window(self, s16_cfg: PipelineConfig) -> None:
        assert 80 <= s16_cfg.f2_b2_min < s16_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, s16_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 flanking primers outside F2-B2."""
        assert s16_cfg.min_region_length >= s16_cfg.f2_b2_max + 40

    def test_output_dir_references_16s(self, s16_cfg: PipelineConfig) -> None:
        assert "16s" in str(s16_cfg.output_dir).lower()

    def test_identity_threshold_is_tighter_than_standard(
        self, s16_cfg: PipelineConfig
    ) -> None:
        """Off-target identity threshold must be tighter than the standard 0.80.

        The 16S rRNA gene has conserved secondary-structure motifs shared between
        Bacteria and Archaea (stem-loop regions fold the same way in both domains).
        A standard 0.80 threshold may fail to flag archaeal cross-reactivity, which
        would be a false-positive in any assay intended to count only bacteria.
        A threshold of >= 0.83 provides additional sensitivity while remaining below
        the 0.90 threshold that would produce excessive noise from distant orthologs.
        """
        assert s16_cfg.min_identity_threshold >= 0.83, (
            "16S universal control should use tighter off-target identity threshold "
            "(>= 0.83) because bacteria and archaea share conserved rRNA structural "
            f"motifs (got {s16_cfg.min_identity_threshold})"
        )

    def test_coverage_threshold_is_standard(self, s16_cfg: PipelineConfig) -> None:
        """Off-target coverage threshold must be in the recommended range."""
        assert 0.70 <= s16_cfg.min_coverage_threshold <= 0.95
