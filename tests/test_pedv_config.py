"""Validate the PEDV N-gene on-farm biosecurity config (config/pedv_N_gene.yaml).

Semantic checks specific to Porcine Epidemic Diarrhea Virus (PEDV) via the
N (nucleocapsid) gene:

- PEDV taxon anchor (NCBI taxon 27317)
- RT-LAMP Tm floor (tm_min >= 63 degC) required for one-step RNA detection
- Entropy threshold appropriate for a porcine coronavirus RNA target
- GC window covers PEDV N-gene GC content (~38-42%) with inter-variant margin
- Tight off-target specificity threshold (>= 0.85) for the porcine coronavirus
  family (TGEV shares the alphacoronavirus N-gene fold; PDCoV is a deltacoronavirus
  co-circulating in the same farrowing-barn niche)
- Amplicon geometry compatible with LAMP primer placement
- Output directory is PEDV-identifiable

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_PEDV_CONFIG = CONFIG_DIR / "pedv_N_gene.yaml"


@pytest.fixture
def pedv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the PEDV N-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_PEDV_CONFIG)


def test_pedv_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_PEDV_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestPedvNGeneConfig:
    """Semantic checks for the PEDV N-gene on-farm biosecurity config.

    Porcine Epidemic Diarrhea Virus (PEDV) is a positive-sense ssRNA
    alphacoronavirus responsible for near-100% mortality in neonatal piglets.
    The N gene (~1326 nt) is the most conserved coding region across all PEDV
    genotypes (G1a, G1b classical strains; G2a, G2b virulent variants including
    USA-Ohio851-2013) and the basis of published RT-PCR and RT-LAMP assays.
    One-step RT-LAMP at 63-65 degC; co-circulates with TGEV and PDCoV.
    """

    def test_target_name(self, pedv_cfg: PipelineConfig) -> None:
        assert pedv_cfg.target_name == "pedv_N_gene"

    def test_gene_is_nucleocapsid(self, pedv_cfg: PipelineConfig) -> None:
        """N (nucleocapsid) must be the target gene.

        The N gene is the most conserved coding region across all PEDV
        genotypes and the basis of published RT-PCR and RT-LAMP assays.
        The spike (S) gene diverges >5% between classical and virulent strains
        and would require separate genotype-specific designs.
        """
        assert pedv_cfg.gene == "N"

    def test_taxon_is_pedv(self, pedv_cfg: PipelineConfig) -> None:
        """Anchor taxon must be PEDV (NCBI taxon 27317).

        NCBI taxon 27317 covers Porcine epidemic diarrhea virus, including
        the classical strain CV777 (AJ38368.1), the USA-Ohio851-2013 virulent
        variant, and contemporary G2b strains responsible for ongoing outbreaks.
        """
        assert pedv_cfg.taxon_id == 27317

    def test_rt_lamp_tm_floor_is_set(self, pedv_cfg: PipelineConfig) -> None:
        """PEDV is an RNA virus; tm_min must meet the RT-LAMP one-step co-activity floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx or WarmStart RTx at
        63-65 degC) requires that core primers have Tm >= 63 degC to maintain
        reverse-transcriptase co-activity.  This is the same floor used for
        PRRSV (Recipe 11), AIV (Recipe 12), FMDV (Recipe 14), and NDV (Recipe 18).
        """
        assert pedv_cfg.tm_min >= 63.0, (
            "PEDV is a positive-sense ssRNA coronavirus: tm_min must be >= 63 degC "
            "for one-step RT-LAMP co-activity with NEB RTx / Bst 2.0 WarmStart"
        )

    def test_tm_max_within_rt_lamp_window(self, pedv_cfg: PipelineConfig) -> None:
        """Tm max must stay within the RT-LAMP optimal window (63-65 degC)."""
        assert pedv_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_entropy_threshold_appropriate_for_rna_virus(self, pedv_cfg: PipelineConfig) -> None:
        """N-gene entropy threshold must reflect RNA-virus inter-variant variability.

        PEDV N gene is under nucleocapsid-assembly constraint, so it is more
        conserved than the spike gene across G1/G2 genotypes, but synonymous
        drift is still higher than for bacterial housekeeping genes.  Thresholds
        below 0.20 bits are too strict and will yield no windows spanning classical
        and virulent strains; above 0.40 bits admits overly variable positions.
        """
        assert pedv_cfg.entropy_threshold >= 0.20, (
            "PEDV N gene diverges across G1/G2 genotypes; "
            "entropy_threshold < 0.20 is too strict and will yield no usable windows"
        )
        assert pedv_cfg.entropy_threshold <= 0.40, (
            "PEDV N gene is under nucleocapsid assembly constraint; "
            "entropy_threshold > 0.40 admits overly variable positions"
        )

    def test_gc_range_covers_pedv_n_gene(self, pedv_cfg: PipelineConfig) -> None:
        """GC window must span PEDV N-gene GC content (~38-42%).

        The PEDV N gene is approximately 38-42% GC across all genotypes.
        The gc_min must be <= 38% and gc_max >= 42% to ensure the primer
        design window includes the conserved N-gene blocks.  Extra margin
        allows the primer3 engine to find qualifying primers in all windows.
        """
        assert pedv_cfg.gc_min <= 38.0, (
            "gc_min must be <= 38% to reach the GC-poor end of the PEDV N-gene primer space"
        )
        assert pedv_cfg.gc_max >= 42.0, (
            "gc_max must be >= 42% to reach the GC-rich end of the PEDV N-gene primer space"
        )

    def test_tight_specificity_threshold_for_coronavirus_family(
        self, pedv_cfg: PipelineConfig
    ) -> None:
        """Off-target identity threshold must be tight for the porcine coronavirus family.

        TGEV (Transmissible Gastroenteritis Virus, Alphacoronavirus 1) shares
        the alphacoronavirus N-gene fold and nucleocapsid-RNA interaction domains
        with PEDV.  PDCoV (Porcine Deltacoronavirus) is a deltacoronavirus that
        co-circulates in the same farrowing-barn niche.  An identity threshold
        >= 0.85 is required to flag potential cross-reactivity from these
        relatives while remaining strict enough to require genuine PEDV sequence.
        """
        assert pedv_cfg.min_identity_threshold >= 0.85, (
            "TGEV and PDCoV share the porcine enteric coronavirus niche; "
            "identity threshold must be >= 0.85 to flag potential cross-reactivity"
        )

    def test_max_sequences_spans_genotype_diversity(self, pedv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover G1a/G1b classical and G2a/G2b virulent variants.

        G1 (classical EU/Asian strains, e.g. CV777) and G2 (virulent North American
        and contemporary Asian strains, e.g. USA-Ohio851-2013, AH2012) share the
        N gene at ~96-99% nt identity.  At least 25 sequences are required to
        represent the diversity needed for a pan-PEDV conserved window.
        """
        assert pedv_cfg.max_sequences >= 25, (
            "PEDV has G1 and G2 genotype clusters; at least 25 sequences are needed "
            "to represent the genotype diversity for cross-variant conservation analysis"
        )

    def test_amplicon_size_within_lamp_window(self, pedv_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= pedv_cfg.f2_b2_min < pedv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, pedv_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert pedv_cfg.min_region_length >= pedv_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_pedv(self, pedv_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as a PEDV run."""
        assert "pedv" in str(pedv_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, pedv_cfg: PipelineConfig) -> None:
        assert 0.70 <= pedv_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, pedv_cfg: PipelineConfig) -> None:
        assert 0.70 <= pedv_cfg.min_coverage_threshold <= 0.95

    def test_tm_range_spans_rt_lamp_window(self, pedv_cfg: PipelineConfig) -> None:
        """Tm window must span >= 1 degC within the 63-65 degC RT-LAMP range."""
        assert pedv_cfg.tm_max - pedv_cfg.tm_min >= 1.0, (
            "tm_max - tm_min must be >= 1.0 degC to allow primer3 design flexibility"
        )
