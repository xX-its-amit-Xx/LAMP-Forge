"""Validate the IBDV VP2 and IBV N-gene configs for the poultry biosecurity panel.

The four-target poultry LAMP panel (``lamp-forge poultry-risk``, Recipe 30)
covers AIV (M gene), NDV (M gene), IBDV (VP2), and IBV (N gene).  Configs for
AIV (``influenza_a_M.yaml``) and NDV (``ndv_M_gene.yaml``) are tested
elsewhere.  This module provides semantic config validation for the two configs
that complete the four-channel poultry panel:

    Config            Target                      NA type   Gene / marker
    ibdv_VP2.yaml     Inf. bursal disease virus   dsRNA     VP2 (outer capsid)
    ibv_N_gene.yaml   Inf. bronchitis virus        RNA+      N (nucleocapsid)

Both targets require RT-LAMP:
- IBDV is a dsRNA birnavirus; after chaotropic lysis the ssRNA template is
  accessible to RTx reverse transcriptase at 63-65 degC.
- IBV is a positive-sense ssRNA gammacoronavirus; standard one-step RT-LAMP
  (Bst 2.0 WarmStart + NEB RTx) at 63-65 degC applies.

Together with AIV and NDV, these two configs complete the four-target dedicated
poultry biosecurity panel for BioVind-style portable field deployment.

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"

_IBDV_CONFIG = CONFIG_DIR / "ibdv_VP2.yaml"
_IBV_CONFIG = CONFIG_DIR / "ibv_N_gene.yaml"

_RT_LAMP_FLOOR = 63.0  # degC; one-step RT-LAMP co-activity floor for NEB RTx


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ibdv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the IBDV VP2 config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_IBDV_CONFIG)


@pytest.fixture
def ibv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the IBV N-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_IBV_CONFIG)


# ---------------------------------------------------------------------------
# Smoke tests -- both configs must parse without error
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "config_path",
    [_IBDV_CONFIG, _IBV_CONFIG],
    ids=["ibdv_VP2", "ibv_N_gene"],
)
def test_poultry_config_parses(
    config_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each poultry config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(config_path)
    assert isinstance(cfg, PipelineConfig)
    assert cfg.target_name, f"{config_path.name}: target_name is empty"
    assert cfg.tm_min < cfg.tm_max, (
        f"{config_path.name}: tm_min={cfg.tm_min} >= tm_max={cfg.tm_max}"
    )
    assert cfg.f2_b2_min < cfg.f2_b2_max, (
        f"{config_path.name}: f2_b2_min={cfg.f2_b2_min} >= f2_b2_max={cfg.f2_b2_max}"
    )
    assert cfg.min_region_length >= cfg.f2_b2_max + 40, (
        f"{config_path.name}: min_region_length too short for f2_b2_max"
    )


# ---------------------------------------------------------------------------
# IBDV VP2 -- dsRNA birnavirus, RT-LAMP
# ---------------------------------------------------------------------------


class TestIbdvVp2Config:
    """Semantic checks for the Infectious Bursal Disease Virus VP2 RT-LAMP config.

    IBDV (Gumboro disease) is a dsRNA birnavirus causing immunosuppression and
    direct mortality in broiler and layer chicks.  Very virulent (vvIBDV)
    strains cause 60-70% direct mortality; classic strains cause
    immunosuppression that amplifies secondary viral and bacterial infections.
    VP2 (outer capsid protein, ~1350 nt) is the group-specific antigen and the
    basis of commercial IBDV ELISA kits and pan-IBDV RT-LAMP assays.  The
    conserved N-terminal structural domain is the primer landing zone; the
    hypervariable region (HVR, aa 206-350) is avoided.
    """

    def test_target_name(self, ibdv_cfg: PipelineConfig) -> None:
        assert ibdv_cfg.target_name == "ibdv_VP2"

    def test_gene_is_vp2(self, ibdv_cfg: PipelineConfig) -> None:
        """VP2 must be the target gene.

        VP2 is the major outer capsid protein and the sole gene carrying the
        group-specific antigen used by all IBDV detection assays.  VP3 (inner
        capsid) is structurally more conserved but not used in pan-IBDV
        diagnostics; VP2 is the correct target to match poultry-risk scoring.
        """
        assert ibdv_cfg.gene == "VP2"

    def test_taxon_is_ibdv(self, ibdv_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Infectious bursal disease virus (NCBI taxon 12227).

        NCBI taxon 12227 covers the Infectious bursal disease virus species,
        including classic strains (e.g. Faragher 52/70, D78), very virulent
        strains (UK661, HK46, Belgium strains), and variant strains
        (Delaware-E, GLS, variant A/B).  Anchoring here retrieves VP2 CDS
        records across the full pathotype spectrum for conservation analysis.
        """
        assert ibdv_cfg.taxon_id == 12227

    def test_rt_lamp_tm_floor_is_set(self, ibdv_cfg: PipelineConfig) -> None:
        """IBDV is dsRNA; RT is required and tm_min must meet the RT-LAMP floor.

        After chaotropic lysis, dsRNA strands are denatured and single-stranded
        RNA template is accessible to RTx reverse transcriptase.  One-step
        RT-LAMP (Bst 2.0 WarmStart + NEB RTx) at 63-65 degC requires core
        primers with Tm >= 63 degC to maintain RTx co-activity -- the same
        floor used for positive-sense ssRNA (PRRSV, FMDV, PEDV) and
        negative-sense ssRNA (AIV, NDV) targets.
        """
        assert ibdv_cfg.tm_min >= _RT_LAMP_FLOOR, (
            "IBDV is a dsRNA birnavirus requiring RT-LAMP: tm_min must be >= "
            f"{_RT_LAMP_FLOOR} degC for RTx co-activity with Bst 2.0 WarmStart "
            f"(got {ibdv_cfg.tm_min})"
        )

    def test_tm_max_within_rt_lamp_window(self, ibdv_cfg: PipelineConfig) -> None:
        """Tm max must stay within the one-step RT-LAMP incubation window."""
        assert ibdv_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_entropy_threshold_reflects_capsid_constraint(self, ibdv_cfg: PipelineConfig) -> None:
        """VP2 entropy threshold must reflect the strong capsid-assembly structural constraint.

        The N-terminal and C-terminal VP2 structural domains are under strong
        functional constraint (T=13 icosahedral capsid geometry).  The
        hypervariable region (HVR) is excluded from primer landing zones.
        A threshold of 0.20-0.30 bits captures the conserved structural blocks
        while excluding the HVR; below 0.15 is too strict (no windows); above
        0.35 admits positions that diverge between classic and vvIBDV strains.
        """
        assert ibdv_cfg.entropy_threshold >= 0.15, (
            "IBDV VP2 entropy_threshold < 0.15 bits is too strict; the capsid "
            "structural blocks still accumulate synonymous drift across strains"
        )
        assert ibdv_cfg.entropy_threshold <= 0.35, (
            "IBDV VP2 entropy_threshold > 0.35 bits admits hypervariable region "
            "positions that diverge between classic and vvIBDV strains"
        )

    def test_gc_range_covers_ibdv_vp2(self, ibdv_cfg: PipelineConfig) -> None:
        """GC window must span IBDV VP2 GC content (~44-48%).

        IBDV VP2 is ~44-48% GC (higher than most ssRNA animal viruses).
        A gc_min of <= 42% provides margin for AT-rich primer candidates
        at conserved block boundaries; gc_max >= 50% covers the GC-rich end.
        """
        assert ibdv_cfg.gc_min <= 42.0, (
            "gc_min must be <= 42% to include primer candidates at the AT-rich "
            "boundaries of IBDV VP2 conserved blocks"
        )
        assert ibdv_cfg.gc_max >= 50.0, (
            "gc_max must be >= 50% to cover the GC-rich end of the IBDV VP2 primer space"
        )

    def test_max_sequences_spans_pathotype_diversity(self, ibdv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover classic, vvIBDV, and variant strains.

        Pan-IBDV detection requires representation of classic strains
        (e.g. D78, Faragher), very virulent strains (UK661, HK46, Belgium),
        and variant strains (Delaware-E, GLS, variant A/B) from multiple
        continents.  At least 25 sequences span the necessary pathotype and
        geographic diversity.
        """
        assert ibdv_cfg.max_sequences >= 25, (
            "IBDV VP2 design needs at least 25 sequences to cover classic, "
            "very virulent, and variant pathotypes from different regions"
        )

    def test_amplicon_size_within_lamp_window(self, ibdv_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= ibdv_cfg.f2_b2_min < ibdv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, ibdv_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert ibdv_cfg.min_region_length >= ibdv_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_ibdv(self, ibdv_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as an IBDV run."""
        assert "ibdv" in str(ibdv_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, ibdv_cfg: PipelineConfig) -> None:
        assert 0.70 <= ibdv_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, ibdv_cfg: PipelineConfig) -> None:
        assert 0.70 <= ibdv_cfg.min_coverage_threshold <= 0.95

    def test_tm_range_allows_design_flexibility(self, ibdv_cfg: PipelineConfig) -> None:
        """Tm window must be wide enough for primer3 to find candidates."""
        assert ibdv_cfg.tm_max - ibdv_cfg.tm_min >= 1.0, (
            "tm_max - tm_min must be >= 1.0 degC to allow primer3 design flexibility"
        )


# ---------------------------------------------------------------------------
# IBV N gene -- positive-sense ssRNA gammacoronavirus, RT-LAMP
# ---------------------------------------------------------------------------


class TestIbvNGeneConfig:
    """Semantic checks for the Infectious Bronchitis Virus N-gene RT-LAMP config.

    IBV (Avian gammacoronavirus) is the most prevalent avian respiratory
    pathogen globally, with dozens of circulating serotypes.  The N gene
    (nucleocapsid protein, ~1230 nt) is the most conserved coding region
    across all IBV serotypes and the target of pan-IBV RT-PCR assays
    (Cavanagh et al. 1999).  IBV is a positive-sense ssRNA virus; one-step
    RT-LAMP (Bst 2.0 WarmStart + NEB RTx) at 63-65 degC applies.
    """

    def test_target_name(self, ibv_cfg: PipelineConfig) -> None:
        assert ibv_cfg.target_name == "ibv_N_gene"

    def test_gene_is_nucleocapsid(self, ibv_cfg: PipelineConfig) -> None:
        """N (nucleocapsid protein) must be the target gene.

        IBV S1 gene diverges at 40-85% nucleotide identity across serotypes,
        making it unsuitable for a pan-IBV detection assay.  The N gene is
        ~80-85% conserved between distant serotypes (Massachusetts vs. 793B
        vs. QX) and is the correct target for a first-pass pan-IBV screen.
        Serotype identification requires S1-gene sequencing downstream.
        """
        assert ibv_cfg.gene == "N"

    def test_taxon_is_ibv(self, ibv_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Avian infectious bronchitis virus (NCBI taxon 11120).

        NCBI taxon 11120 covers the Avian infectious bronchitis virus (IBV)
        species including Massachusetts (M41), 793B/4/91, QX, D388, Georgia
        08/1-1, Italy 02, and other globally distributed serotypes.  Sequences
        from multiple serotypes are needed to identify the conserved N-gene
        block common to all strains.
        """
        assert ibv_cfg.taxon_id == 11120

    def test_rt_lamp_tm_floor_is_set(self, ibv_cfg: PipelineConfig) -> None:
        """IBV is positive-sense ssRNA; tm_min must meet the RT-LAMP co-activity floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx at 63-65 degC) requires
        that core primers (F3/B3/FIP/BIP) have Tm >= 63 degC to maintain
        reverse-transcriptase co-activity.  This is the same floor used for
        PRRSV (Recipe 11), AIV (Recipe 12), FMDV (Recipe 14), NDV (Recipe 18),
        and PEDV (Recipe 27).
        """
        assert ibv_cfg.tm_min >= _RT_LAMP_FLOOR, (
            "IBV is a positive-sense ssRNA gammacoronavirus requiring RT-LAMP: "
            f"tm_min must be >= {_RT_LAMP_FLOOR} degC for RTx co-activity "
            f"with Bst 2.0 WarmStart (got {ibv_cfg.tm_min})"
        )

    def test_tm_max_within_rt_lamp_window(self, ibv_cfg: PipelineConfig) -> None:
        """Tm max must stay within the one-step RT-LAMP incubation window."""
        assert ibv_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_entropy_threshold_appropriate_for_ibv_serotype_diversity(
        self, ibv_cfg: PipelineConfig
    ) -> None:
        """N-gene entropy threshold must reflect IBV serotype diversity.

        IBV N gene is ~80-85% conserved between distant serotypes (Mass vs.
        793B vs. QX).  Synonymous drift between serotypes is higher than for
        PRRSV or NDV, which have fewer circulating lineages.  A threshold of
        0.20-0.40 bits captures the conserved N-gene core without excluding
        all windows across serotypes.
        """
        assert ibv_cfg.entropy_threshold >= 0.20, (
            "IBV N gene diverges across serotypes (Mass, 793B, QX, D388, etc.); "
            "entropy_threshold < 0.20 bits is too strict and will yield no usable windows"
        )
        assert ibv_cfg.entropy_threshold <= 0.40, (
            "IBV N gene is under nucleocapsid-assembly constraint; "
            "entropy_threshold > 0.40 admits overly variable positions"
        )

    def test_gc_range_covers_ibv_n_gene(self, ibv_cfg: PipelineConfig) -> None:
        """GC window must span IBV N-gene GC content (~37-42%).

        IBV N gene is approximately 37-42% GC (gammacoronavirus genomes are
        AT-leaning compared with mammalian betacoronaviruses).  A gc_min of
        <= 37% provides margin for AT-rich primer candidates; gc_max >= 45%
        covers the GC-rich end of N-gene primer space across serotypes.
        """
        assert ibv_cfg.gc_min <= 37.0, (
            "gc_min must be <= 37% to accommodate AT-rich IBV N-gene primer "
            "candidates at the boundaries of conserved blocks"
        )
        assert ibv_cfg.gc_max >= 45.0, (
            "gc_max must be >= 45% to cover the GC-rich end of IBV N-gene primer space"
        )

    def test_gc_window_wide_enough_for_serotype_diversity(self, ibv_cfg: PipelineConfig) -> None:
        """GC window must be wide enough for IBV N-gene inter-serotype GC variation."""
        assert ibv_cfg.gc_max - ibv_cfg.gc_min >= 15.0, (
            "IBV N-gene GC content varies across serotypes; gc window must be "
            "at least 15 units wide to ensure primer3 finds candidates in all blocks"
        )

    def test_max_sequences_spans_serotype_diversity(self, ibv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover major globally circulating serotypes.

        IBV has at least 6 major globally circulating serotype groups:
        Massachusetts, 793B/4/91, QX, D388, Georgia 08, and Italy 02.
        Additional regional serotypes circulate in Asia, the Middle East, and
        Africa.  At least 30 sequences are needed to represent the serotype
        diversity required for pan-IBV N-gene conservation analysis.
        """
        assert ibv_cfg.max_sequences >= 30, (
            "IBV has many circulating serotypes (Mass, 793B, QX, D388, Georgia 08, "
            "etc.); at least 30 sequences are needed for pan-IBV conservation analysis"
        )

    def test_tight_specificity_for_coronavirus_family(self, ibv_cfg: PipelineConfig) -> None:
        """Off-target identity threshold must be tight for the gammacoronavirus family.

        Avian coronavirus sequences from avian metapneumovirus (AMPV) and
        duck coronaviruses may share short N-gene motifs.  A threshold >= 0.83
        is required to flag potential cross-reactivity from these relatives.
        """
        assert ibv_cfg.min_identity_threshold >= 0.83, (
            "Gammacoronavirus N-gene structural motifs are shared across avian "
            "coronaviruses; identity threshold must be >= 0.83 to flag cross-reactivity"
        )

    def test_amplicon_size_within_lamp_window(self, ibv_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= ibv_cfg.f2_b2_min < ibv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, ibv_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert ibv_cfg.min_region_length >= ibv_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_ibv(self, ibv_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as an IBV run."""
        assert "ibv" in str(ibv_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, ibv_cfg: PipelineConfig) -> None:
        assert 0.70 <= ibv_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, ibv_cfg: PipelineConfig) -> None:
        assert 0.70 <= ibv_cfg.min_coverage_threshold <= 0.95
