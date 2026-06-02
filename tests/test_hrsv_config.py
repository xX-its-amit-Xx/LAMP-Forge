"""Validate the human RSV N-gene POC LAMP config (config/hrsv_N_gene.yaml).

Human Respiratory Syncytial Virus (hRSV) is the leading cause of acute lower
respiratory tract illness in infants and a major pathogen in elderly and
immunocompromised adults -- key populations at BioVind's target POC sites.

Semantic checks specific to hRSV via the N (nucleocapsid) gene:

- hRSV species taxon anchor (NCBI taxon 11250, Human orthopneumovirus / RSV-A)
- RT-LAMP Tm floor (tm_min >= 63 degC) required for one-step RNA detection
- Entropy threshold appropriate for a pneumovirus RNA target (~77% N-gene
  identity between RSV-A and RSV-B subtypes means a pan-RSV window requires
  a moderately relaxed threshold)
- GC window covers hRSV N-gene GC content (~37-42%) with inter-subtype margin
- Tight off-target specificity threshold (>= 0.85) due to human metapneumovirus
  (hMPV) sharing the Pneumoviridae family
- Amplicon geometry compatible with LAMP primer placement
- Output directory is hRSV-identifiable

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_HRSV_CONFIG = CONFIG_DIR / "hrsv_N_gene.yaml"

_RT_LAMP_FLOOR = 63.0  # degC; one-step RT-LAMP co-activity floor for NEB RTx


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def hrsv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the hRSV N-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_HRSV_CONFIG)


# ---------------------------------------------------------------------------
# Basic parse
# ---------------------------------------------------------------------------


def test_hrsv_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_HRSV_CONFIG)
    assert isinstance(cfg, PipelineConfig)


# ---------------------------------------------------------------------------
# Semantic checks
# ---------------------------------------------------------------------------


class TestHrsvNGeneConfig:
    """Semantic checks for the human RSV N-gene POC LAMP config.

    Human Respiratory Syncytial Virus is a negative-sense ssRNA pneumovirus
    (Pneumoviridae).  Two subtypes co-circulate: RSV-A (species anchor txid
    11250) and RSV-B (~77% N-gene nt identity to RSV-A).  A pan-RSV assay
    must detect both, making the N gene the correct target.  One-step RT-LAMP
    at 63-65 degC; closest off-target is human metapneumovirus (hMPV).
    """

    def test_target_name(self, hrsv_cfg: PipelineConfig) -> None:
        assert hrsv_cfg.target_name == "hrsv_N_gene"

    def test_gene_is_nucleocapsid(self, hrsv_cfg: PipelineConfig) -> None:
        """N (nucleocapsid) must be the target gene.

        The N gene (~1203 nt CDS) is the most conserved coding region across
        RSV-A and RSV-B subtypes (~77% nt identity) and is the basis of
        published one-step RT-LAMP assays for pan-RSV detection (Nie et al.
        2017; Zhang et al. 2021).
        """
        assert hrsv_cfg.gene == "N"

    def test_taxon_is_human_rsv(self, hrsv_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Human orthopneumovirus (RSV-A; NCBI txid 11250).

        NCBI taxon 11250 covers Human orthopneumovirus (RSV-A), the subtype
        used as the sequence anchor.  RSV-B (txid 11251) accessions should be
        added to the accessions list to capture cross-subtype conservation.
        """
        assert hrsv_cfg.taxon_id == 11250

    def test_rt_lamp_tm_floor_is_set(self, hrsv_cfg: PipelineConfig) -> None:
        """hRSV is RNA; tm_min must meet the RT-LAMP one-step co-activity floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx at 63-65 degC) requires
        core primers (F3/B3/FIP/BIP) to have Tm >= 63 degC to maintain
        reverse-transcriptase co-activity.  This is the same floor enforced for
        PRRSV (Recipe 11), AIV (Recipe 12), FMDV (Recipe 14), NDV (Recipe 18),
        PEDV (Recipe 27), and CSFV (Recipe 29).
        """
        assert hrsv_cfg.tm_min >= _RT_LAMP_FLOOR, (
            f"hRSV is negative-sense ssRNA: tm_min must be >= {_RT_LAMP_FLOOR} degC "
            "for one-step RT-LAMP co-activity with NEB RTx / Bst 2.0 WarmStart"
        )

    def test_tm_max_within_rt_lamp_window(self, hrsv_cfg: PipelineConfig) -> None:
        """Tm max must stay within the RT-LAMP optimal incubation window."""
        assert hrsv_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_tm_range_spans_at_least_one_degree(self, hrsv_cfg: PipelineConfig) -> None:
        """Tm window must span >= 1 degC to give primer3 design flexibility."""
        assert hrsv_cfg.tm_max - hrsv_cfg.tm_min >= 1.0, (
            "tm_max - tm_min must be >= 1.0 degC to allow primer3 design flexibility "
            "within the RT-LAMP window"
        )

    def test_entropy_threshold_appropriate_for_rna_virus(self, hrsv_cfg: PipelineConfig) -> None:
        """N-gene entropy threshold must accommodate RSV inter-subtype variability.

        RSV-A and RSV-B share ~77% N-gene nt identity.  A pan-RSV conservation
        analysis requires a threshold loose enough to find windows conserved across
        both subtypes but tight enough to exclude the highly variable attachment (G)
        gene positions.  Values below 0.20 bits are too strict for cross-subtype
        design; above 0.40 bits allows overly variable positions into primer regions.
        """
        assert hrsv_cfg.entropy_threshold >= 0.20, (
            "RSV-A and RSV-B share only ~77% N-gene nt identity; "
            "entropy_threshold < 0.20 is too strict for pan-RSV coverage"
        )
        assert hrsv_cfg.entropy_threshold <= 0.40, (
            "entropy_threshold > 0.40 admits overly variable N-gene positions; "
            "use 0.30 for a balance between cross-subtype coverage and primer stability"
        )

    def test_gc_range_covers_hrsv_n_gene(self, hrsv_cfg: PipelineConfig) -> None:
        """GC window must span hRSV N-gene GC content (~37-42% across subtypes).

        The hRSV N gene is approximately 37-42% GC (RSV is one of the more
        AT-rich human respiratory RNA viruses).  gc_min must be <= 37% and
        gc_max >= 42% to ensure the primer design engine can find qualifying
        primers in the conserved N-gene blocks for both RSV-A and RSV-B.
        """
        assert hrsv_cfg.gc_min <= 37.0, (
            "gc_min must be <= 37% to reach the GC-poor end of the hRSV N-gene primer space"
        )
        assert hrsv_cfg.gc_max >= 42.0, (
            "gc_max must be >= 42% to reach the GC-rich end of the hRSV N-gene primer space"
        )

    def test_tight_specificity_threshold_for_pneumovirus_family(
        self, hrsv_cfg: PipelineConfig
    ) -> None:
        """Off-target identity threshold must be tight for the Pneumoviridae family.

        Human metapneumovirus (hMPV; Pneumoviridae, Metapneumovirus genus) is
        the closest non-target relative and shares structural N-protein motifs
        with RSV.  An identity threshold >= 0.85 is required to flag potential
        cross-reactivity while remaining strict enough to require genuine hRSV
        N-gene sequence before flagging an off-target hit.
        """
        assert hrsv_cfg.min_identity_threshold >= 0.85, (
            "hMPV shares the Pneumoviridae family with hRSV; "
            "identity threshold must be >= 0.85 to flag potential N-gene cross-reactivity"
        )

    def test_max_sequences_spans_subtype_diversity(self, hrsv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover RSV-A clades and RSV-B BA lineages.

        RSV-A has three major clades (A1, A2, A2b); RSV-B circulates as BA
        lineages (BA7-BA10).  At least 30 sequences are required to represent
        the diversity needed for a reliable pan-RSV conserved window in the N gene.
        """
        assert hrsv_cfg.max_sequences >= 30, (
            "hRSV has RSV-A and RSV-B subtypes with internal clade structure; "
            "at least 30 sequences are needed to find cross-subtype conserved windows"
        )

    def test_amplicon_size_within_lamp_window(self, hrsv_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= hrsv_cfg.f2_b2_min < hrsv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, hrsv_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert hrsv_cfg.min_region_length >= hrsv_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_hrsv(self, hrsv_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as an hRSV run."""
        assert "hrsv" in str(hrsv_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, hrsv_cfg: PipelineConfig) -> None:
        assert 0.70 <= hrsv_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, hrsv_cfg: PipelineConfig) -> None:
        assert 0.70 <= hrsv_cfg.min_coverage_threshold <= 0.95

    def test_window_size_appropriate_for_rna_virus(self, hrsv_cfg: PipelineConfig) -> None:
        """Sliding window must be narrow enough to resolve short conserved blocks.

        RNA viruses have scattered conserved blocks rather than long uniform
        stretches.  A window size of 20-30 nt resolves the constrained structural
        and functional domains of the N protein RNA-binding region.
        """
        assert 15 <= hrsv_cfg.window_size <= 35, (
            "window_size should be 15-35 nt for RNA-virus N-gene conservation; "
            "wider windows smear short conserved blocks into false positives"
        )
