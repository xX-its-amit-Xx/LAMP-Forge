"""Validate the CSFV NS5B on-farm biosecurity config (config/csfv_NS5B.yaml).

Semantic checks specific to Classical Swine Fever Virus (CSFV) via the
NS5B (RNA-dependent RNA polymerase) gene:

- CSFV taxon anchor (NCBI taxon 11096, Pestivirus A)
- RT-LAMP Tm floor (tm_min >= 63 degC) required for one-step RNA detection
- Entropy threshold appropriate for a pestivirus RdRp target
- GC window covers CSFV NS5B GC content (~51-55%) with inter-subtype margin
- Tight off-target specificity threshold (>= 0.85) required because BVDV
  (Bovine viral diarrhea virus) is a closely related pestivirus sharing
  ~50-60% NS5B nucleotide identity with CSFV
- Enough input sequences to span all three CSFV genotype groups (1/2/3)
- Amplicon geometry compatible with LAMP primer placement

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_CSFV_CONFIG = CONFIG_DIR / "csfv_NS5B.yaml"


@pytest.fixture
def csfv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the CSFV NS5B config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_CSFV_CONFIG)


def test_csfv_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_CSFV_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestCsfvNS5BConfig:
    """Semantic checks for the CSFV NS5B on-farm biosecurity config.

    Classical Swine Fever Virus (CSFV; hog cholera) is a WOAH-listed
    notifiable pestivirus (Flaviviridae, Pestivirus A) causing near-100%
    mortality in naive herds.  The NS5B gene (~1.8 kb) encodes the
    RNA-directed RNA polymerase and is the most conserved ORF across
    genotype groups 1/2/3 and subtypes 1.1-3.4.  CSFV is a positive-sense
    ssRNA virus; detection requires one-step RT-LAMP at 63-65 degC.
    """

    def test_target_name(self, csfv_cfg: PipelineConfig) -> None:
        assert csfv_cfg.target_name == "csfv_NS5B"

    def test_gene_is_NS5B(self, csfv_cfg: PipelineConfig) -> None:
        """NS5B (RNA-directed RNA polymerase) must be the target gene.

        NS5B is the most conserved ORF across all CSFV genotype groups (1/2/3)
        and the target of choice in published CSFV RT-PCR and RT-LAMP assays
        (Wang et al. 2018; Rodriguez-Prieto et al. 2017).  The alternative E2
        glycoprotein gene is too variable across subtypes for a pan-CSFV screen.
        """
        assert csfv_cfg.gene == "NS5B"

    def test_taxon_is_csfv(self, csfv_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Classical swine fever virus (NCBI taxon 11096).

        NCBI taxon 11096 covers Pestivirus A (Classical swine fever virus),
        including genotype 1 (Brescia, LOM), genotype 2 (Paderborn, Alfort),
        and genotype 3 (Alfort/187) strains spanning the full subtype diversity
        needed for a pan-CSFV conserved-window analysis.
        """
        assert csfv_cfg.taxon_id == 11096

    def test_rt_lamp_tm_floor_is_set(self, csfv_cfg: PipelineConfig) -> None:
        """CSFV is an RNA virus; tm_min must meet the RT-LAMP one-step co-activity floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx or WarmStart RTx at
        63-65 degC) requires that core primers have Tm >= 63 degC to maintain
        reverse-transcriptase co-activity.  This is the same floor used for
        PRRSV (Recipe 11), AIV (Recipe 12), FMDV (Recipe 14), and NDV (Recipe 18).
        """
        assert csfv_cfg.tm_min >= 63.0, (
            "CSFV is a positive-sense ssRNA pestivirus: tm_min must be >= 63 degC "
            "for one-step RT-LAMP co-activity with NEB RTx / Bst 2.0 WarmStart"
        )

    def test_tm_max_within_rt_lamp_window(self, csfv_cfg: PipelineConfig) -> None:
        """Tm max must stay within the RT-LAMP optimal window (63-65 degC)."""
        assert csfv_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_entropy_threshold_appropriate_for_rna_virus(self, csfv_cfg: PipelineConfig) -> None:
        """NS5B entropy threshold must reflect RNA-virus inter-subtype variability.

        The CSFV NS5B gene is under strong RdRp catalytic constraint (GDD motif,
        finger / palm / thumb subdomains), which makes it more conserved than
        most other pestivirus genes.  However, synonymous drift across the three
        genotype groups is still higher than for bacterial housekeeping genes.
        Thresholds below 0.20 bits may yield no windows spanning genotypes 1-3;
        above 0.40 bits admits positions too variable for cross-subtype primers.
        """
        assert csfv_cfg.entropy_threshold >= 0.20, (
            "CSFV NS5B diverges across genotype groups 1/2/3; "
            "entropy_threshold < 0.20 is too strict and will yield no usable windows"
        )
        assert csfv_cfg.entropy_threshold <= 0.40, (
            "CSFV NS5B is under RdRp catalytic constraint; "
            "entropy_threshold > 0.40 admits overly variable positions"
        )

    def test_gc_range_covers_csfv_ns5b(self, csfv_cfg: PipelineConfig) -> None:
        """GC window must span CSFV NS5B GC content (~51-55%) with subtype margin.

        CSFV NS5B is approximately 51-55% GC across genotype groups 1/2/3.
        The gc_min must be <= 46% and gc_max >= 60% to ensure the primer
        design window includes conserved NS5B blocks across all subtypes.
        An extra margin allows the primer3 engine to find qualifying primers
        in every conserved window across diverse strains.
        """
        assert csfv_cfg.gc_min <= 46.0, (
            "gc_min must be <= 46% to reach the GC-poor end of the CSFV NS5B primer space"
        )
        assert csfv_cfg.gc_max >= 60.0, (
            "gc_max must be >= 60% to reach the GC-rich end of the CSFV NS5B primer space"
        )

    def test_tight_specificity_threshold_for_pestivirus_family(
        self, csfv_cfg: PipelineConfig
    ) -> None:
        """Off-target identity threshold must be tight for the pestivirus family.

        Bovine viral diarrhea virus (BVDV, Pestivirus B) is the closest
        relative to CSFV in the off-target panel, sharing ~50-60% NS5B nt
        identity.  A threshold >= 0.85 is required to flag potential BVDV
        cross-reactivity while remaining strict enough to demand genuine CSFV
        sequence for a positive call.  Loose thresholds (< 0.80) will generate
        false flags for any BVDV-positive cattle in a mixed-species farm sample.
        """
        assert csfv_cfg.min_identity_threshold >= 0.85, (
            "BVDV (Pestivirus B) shares ~50-60% NS5B nt identity with CSFV; "
            "identity threshold must be >= 0.85 to flag potential cross-reactivity"
        )

    def test_max_sequences_spans_genotype_diversity(self, csfv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover genotype groups 1, 2, and 3 and major subtypes.

        CSFV has three genotype groups (1, 2, 3) with subtypes 1.1, 1.2, 1.3,
        2.1, 2.2, 2.3, and 3.1-3.4.  At least 30 sequences are required to
        represent this diversity for cross-genotype conservation analysis.
        Subtype 2.1 accounts for most contemporary Asian outbreaks; genotype 1
        includes the Brescia and LOM reference strains widely used in
        validation studies.
        """
        assert csfv_cfg.max_sequences >= 30, (
            "CSFV has three genotype groups with multiple subtypes; "
            "at least 30 sequences needed for cross-genotype conservation analysis"
        )

    def test_amplicon_size_within_lamp_window(self, csfv_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= csfv_cfg.f2_b2_min < csfv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, csfv_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert csfv_cfg.min_region_length >= csfv_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_csfv(self, csfv_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as a CSFV run."""
        assert "csfv" in str(csfv_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, csfv_cfg: PipelineConfig) -> None:
        assert 0.70 <= csfv_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, csfv_cfg: PipelineConfig) -> None:
        assert 0.70 <= csfv_cfg.min_coverage_threshold <= 0.95

    def test_tm_range_spans_rt_lamp_window(self, csfv_cfg: PipelineConfig) -> None:
        """Tm window must span >= 1 degC within the 63-65 degC RT-LAMP range."""
        assert csfv_cfg.tm_max - csfv_cfg.tm_min >= 1.0, (
            "tm_max - tm_min must be >= 1.0 degC to allow primer3 design flexibility"
        )

    def test_window_size_suitable_for_rna_virus(self, csfv_cfg: PipelineConfig) -> None:
        """Conservation window must be narrow enough to catch constrained RdRp blocks.

        For RNA viruses, a shorter sliding window (20-30 nt) better captures the
        structurally constrained catalytic motifs (e.g. GDD polymerase motif in
        NS5B) while ignoring the more variable flanking loops.  A window > 35 nt
        dilutes the signal from these tight conservation islands.
        """
        assert csfv_cfg.window_size <= 35, (
            "window_size > 35 may dilute narrow conservation peaks in CSFV NS5B; "
            "use <= 30 nt for RNA-virus RdRp targets"
        )
