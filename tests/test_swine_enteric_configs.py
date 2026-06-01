"""Validate the two remaining swine neonatal enteric LAMP configs.

The swine neonatal enteric panel (``lamp-forge swine-enteric-risk``) monitors
three pathogens -- PEDV (config/pedv_N_gene.yaml, tested separately in
test_pedv_config.py), PDCoV, and Porcine Rotavirus A.  This module provides
semantic config validation for the two configs that complete the three-channel
panel:

    Config                  Target                    NA type  Gene / marker
    pdcov_N_gene.yaml       Porcine deltacoronavirus  RNA      N (nucleocapsid)
    porcine_rota_a.yaml     Porcine Rotavirus A       dsRNA    VP6 (inner capsid)

Both targets require RT-LAMP (PDCoV is positive-sense ssRNA; Rotavirus A is
dsRNA and, after heat-denaturation or chaotropic lysis, primed by RTx at
63-65 degC).  Together with PEDV they form the complete swine neonatal enteric
diagnostic panel for BioVind-style portable deployment.

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"

_PDCOV_CONFIG = CONFIG_DIR / "pdcov_N_gene.yaml"
_ROTA_CONFIG = CONFIG_DIR / "porcine_rota_a.yaml"

_RT_LAMP_FLOOR = 63.0  # degC; one-step RT-LAMP co-activity floor for NEB RTx


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pdcov_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the PDCoV N-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_PDCOV_CONFIG)


@pytest.fixture
def rota_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the Porcine Rotavirus A VP6 config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_ROTA_CONFIG)


# ---------------------------------------------------------------------------
# Smoke tests -- both configs must parse without error
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "config_path",
    [_PDCOV_CONFIG, _ROTA_CONFIG],
    ids=["pdcov", "porcine_rota_a"],
)
def test_swine_enteric_config_parses(
    config_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each swine enteric config must load to a PipelineConfig without raising ConfigError."""
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
# PDCoV N gene -- RNA virus (deltacoronavirus), RT-LAMP
# ---------------------------------------------------------------------------


class TestPdcovNGeneConfig:
    """Semantic checks for the Porcine Deltacoronavirus (PDCoV) N-gene RT-LAMP config.

    PDCoV (Porcine deltacoronavirus; Swine enteric deltacoronavirus) is a
    positive-sense ssRNA deltacoronavirus causing diarrhea and vomiting in
    neonatal piglets (40-80% CFR in piglets under 1 week).  The N gene is the
    most conserved coding region across North American and Asian strains and
    is the basis of published PDCoV RT-PCR and RT-LAMP assays.  One-step
    RT-LAMP at 63-65 degC.  PDCoV is phylogenetically distinct from PEDV
    (Alphacoronavirus 1); the two viruses co-circulate in farrowing barns.
    """

    def test_target_name(self, pdcov_cfg: PipelineConfig) -> None:
        assert pdcov_cfg.target_name == "pdcov_N_gene"

    def test_gene_is_nucleocapsid(self, pdcov_cfg: PipelineConfig) -> None:
        """N (nucleocapsid) must be the target gene."""
        assert pdcov_cfg.gene == "N"

    def test_taxon_is_pdcov(self, pdcov_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Porcine deltacoronavirus (NCBI taxon 1542401).

        NCBI taxon 1542401 covers all PDCoV (Swine enteric deltacoronavirus)
        strains including the USA-IL2014 prototype, Korean KNU14-04, Chinese
        CHN-GD-2015, and Vietnamese-2016 strains that together represent the
        phylogenetic breadth of circulating PDCoV lineages.
        """
        assert pdcov_cfg.taxon_id == 1542401

    def test_rt_lamp_tm_floor_is_set(self, pdcov_cfg: PipelineConfig) -> None:
        """PDCoV is a positive-sense ssRNA virus; tm_min must meet the RT-LAMP floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx at 63-65 degC) requires
        core primers with Tm >= 63 degC to maintain reverse-transcriptase
        co-activity.  This is the same floor used for PEDV (Recipe 27), PRRSV
        (Recipe 11), FMDV (Recipe 14), AIV (Recipe 12), and NDV (Recipe 18).
        """
        assert pdcov_cfg.tm_min >= _RT_LAMP_FLOOR, (
            "PDCoV is a positive-sense ssRNA virus: tm_min must be >= "
            f"{_RT_LAMP_FLOOR} degC for one-step RT-LAMP co-activity "
            f"(got {pdcov_cfg.tm_min})"
        )

    def test_tm_max_within_rt_lamp_window(self, pdcov_cfg: PipelineConfig) -> None:
        """Tm max must stay within the RT-LAMP optimal window."""
        assert pdcov_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_entropy_threshold_appropriate_for_rna_virus(self, pdcov_cfg: PipelineConfig) -> None:
        """N-gene entropy threshold must reflect RNA-virus inter-variant variability.

        PDCoV N gene is constrained by nucleocapsid-assembly function but
        accumulates synonymous drift across North American and Asian clades.
        Thresholds below 0.20 bits are too strict and exclude all usable windows;
        above 0.40 bits admits positions that vary across strains.
        """
        assert pdcov_cfg.entropy_threshold >= 0.20, (
            "PDCoV N gene diverges across North American and Asian clades; "
            "entropy_threshold < 0.20 bits is too strict for a cross-clade design"
        )
        assert pdcov_cfg.entropy_threshold <= 0.40, (
            "PDCoV N gene is under nucleocapsid-assembly constraint; "
            "entropy_threshold > 0.40 admits overly variable positions"
        )

    def test_gc_range_covers_pdcov_n_gene(self, pdcov_cfg: PipelineConfig) -> None:
        """GC window must span PDCoV N-gene GC content (~50-55%).

        PDCoV N gene is approximately 50-55% GC across circulating strains.
        The design window must encompass this range with margin for the primer3
        engine to find qualifying primers at conserved block boundaries.
        """
        assert pdcov_cfg.gc_min <= 45.0, (
            "gc_min must be <= 45% to ensure primer candidates are found in all "
            "conserved PDCoV N-gene blocks"
        )
        assert pdcov_cfg.gc_max >= 55.0, (
            "gc_max must be >= 55% to cover the GC-rich end of the PDCoV N-gene primer space"
        )

    def test_tight_specificity_threshold_for_coronavirus_family(
        self, pdcov_cfg: PipelineConfig
    ) -> None:
        """Off-target identity threshold must be tight for the porcine coronavirus family.

        TGEV (Alphacoronavirus 1) and PEDV co-circulate in the same farrowing-barn
        niche.  A threshold >= 0.85 is required to flag potential cross-reactivity
        from these relatives while distinguishing true PDCoV sequence.
        """
        assert pdcov_cfg.min_identity_threshold >= 0.85, (
            "TGEV and PEDV share the porcine enteric coronavirus niche; "
            "identity threshold must be >= 0.85 to flag potential cross-reactivity"
        )

    def test_max_sequences_spans_geographic_diversity(self, pdcov_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover USA, Asian (China/Korea/Thailand) strains.

        PDCoV strains from the USA (prototype USA-IL2014), China (CHN-GD-2015),
        South Korea (KNU14-04), Thailand, and other regions show moderate
        phylogenetic divergence.  At least 25 sequences are needed to represent
        the inter-clade diversity for conservation analysis.
        """
        assert pdcov_cfg.max_sequences >= 25, (
            "PDCoV strains span North American and Asian clades; at least 25 "
            "sequences are required to represent the cross-clade diversity"
        )

    def test_amplicon_size_within_lamp_window(self, pdcov_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= pdcov_cfg.f2_b2_min < pdcov_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, pdcov_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert pdcov_cfg.min_region_length >= pdcov_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_pdcov(self, pdcov_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as a PDCoV run."""
        assert "pdcov" in str(pdcov_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, pdcov_cfg: PipelineConfig) -> None:
        assert 0.70 <= pdcov_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, pdcov_cfg: PipelineConfig) -> None:
        assert 0.70 <= pdcov_cfg.min_coverage_threshold <= 0.95


# ---------------------------------------------------------------------------
# Porcine Rotavirus A VP6 -- dsRNA virus, RT-LAMP after denaturation
# ---------------------------------------------------------------------------


class TestPorcineRotaAConfig:
    """Semantic checks for the Porcine Rotavirus A VP6 RT-LAMP config.

    Porcine Rotavirus A is the most prevalent cause of neonatal piglet diarrhea
    worldwide (10-40% CFR, all G/P genotype combinations).  VP6 (segment 6,
    ~1356 nt) encodes the inner capsid protein and the group A-specific antigen;
    it is the most conserved Rotavirus A segment and the target of all commercial
    ELISA kits and the basis of VP6-RT-LAMP assays.  RT-LAMP at 63-65 degC
    works on dsRNA after endogenous denaturation during chaotropic lysis.
    """

    def test_target_name(self, rota_cfg: PipelineConfig) -> None:
        assert rota_cfg.target_name == "porcine_rota_a"

    def test_gene_is_vp6(self, rota_cfg: PipelineConfig) -> None:
        """VP6 must be the target gene.

        VP6 is the group A-specific antigen and the most conserved segment
        across all Rotavirus A genotypes, making it the basis of all standard
        detection assays.  VP4 and VP7 are under immune selection and diverge
        across the numerous G and P genotype combinations circulating in pigs.
        """
        assert rota_cfg.gene == "VP6"

    def test_taxon_is_rotavirus_a(self, rota_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Rotavirus A (NCBI taxon 10359).

        NCBI taxon 10359 covers the Rotavirus A species, which includes all
        porcine G and P genotype combinations (G1-G5, G9, G11; P[6], P[7],
        P[13]).  Porcine rotaviruses are not assigned a separate species-level
        taxon node; anchoring at Rotavirus A with gene=VP6 retrieves porcine,
        bovine, and human VP6 sequences from NCBI, providing a broad
        conservation anchor that is then refined by the porcine-focused
        accession set.
        """
        assert rota_cfg.taxon_id == 10359

    def test_rt_lamp_tm_floor_is_set(self, rota_cfg: PipelineConfig) -> None:
        """Rotavirus A is dsRNA; RT is required and tm_min must meet the RT-LAMP floor.

        Heat denaturation of dsRNA (standard in chaotropic RNA extraction kits)
        releases single-stranded RNA template.  RTx reverse transcriptase then
        primes at 63-65 degC, co-active with Bst 2.0 WarmStart.  tm_min must
        be >= 63 degC to maintain RTx co-activity, the same requirement as
        positive-sense ssRNA (PRRSV, FMDV) and negative-sense ssRNA (AIV, NDV)
        targets.
        """
        assert rota_cfg.tm_min >= _RT_LAMP_FLOOR, (
            "Porcine Rotavirus A is a dsRNA virus requiring RT-LAMP: tm_min must "
            f"be >= {_RT_LAMP_FLOOR} degC for RTx co-activity with Bst 2.0 WarmStart "
            f"(got {rota_cfg.tm_min})"
        )

    def test_tm_max_within_rt_lamp_window(self, rota_cfg: PipelineConfig) -> None:
        """Tm max must stay within the RT-LAMP optimal window."""
        assert rota_cfg.tm_max <= 65.0, (
            "tm_max > 65 degC exceeds the optimal one-step RT-LAMP incubation temperature"
        )

    def test_entropy_threshold_appropriate_for_vp6(self, rota_cfg: PipelineConfig) -> None:
        """VP6 entropy threshold must reflect its high conservation within Rotavirus A.

        VP6 is the most conserved Rotavirus A segment; its inner-capsid structural
        constraint limits synonymous drift more than the outer-capsid VP4 and VP7
        genes.  A threshold of 0.15-0.35 bits is appropriate: below 0.15 is too
        strict and admits no windows; above 0.35 admits positions that vary across
        G/P genotype combinations.
        """
        assert rota_cfg.entropy_threshold >= 0.15, (
            "VP6 entropy_threshold < 0.15 bits is too strict; no usable conserved "
            "windows will be identified across porcine G/P genotype combinations"
        )
        assert rota_cfg.entropy_threshold <= 0.35, (
            "VP6 is the most conserved Rotavirus A gene; entropy_threshold > 0.35 "
            "admits overly variable positions outside the group A-specific antigen core"
        )

    def test_gc_range_covers_rotavirus_vp6(self, rota_cfg: PipelineConfig) -> None:
        """GC window must span porcine Rotavirus A VP6 GC content (~40-45%).

        Porcine Rotavirus A VP6 is approximately 40-45% GC.  A gc_min of <= 38%
        provides margin for AT-rich primer candidates at block boundaries, and
        gc_max >= 45% accommodates the full GC range.
        """
        assert rota_cfg.gc_min <= 38.0, (
            "gc_min must be <= 38% to include AT-rich primer candidates at the "
            "boundaries of the VP6 conserved blocks"
        )
        assert rota_cfg.gc_max >= 45.0, (
            "gc_max must be >= 45% to cover the GC-rich end of the VP6 primer space"
        )

    def test_gc_window_is_wide_enough(self, rota_cfg: PipelineConfig) -> None:
        """GC window must be wide enough for VP6 inter-genotype diversity."""
        assert rota_cfg.gc_max - rota_cfg.gc_min >= 15.0, (
            "VP6 GC content spans ~10 percentage points across porcine genotypes; "
            "gc window must be at least 15 units wide to ensure primer3 finds candidates"
        )

    def test_max_sequences_spans_genotype_diversity(self, rota_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover major porcine G and P genotype combinations.

        Porcine Rotavirus A circulates as at least 7 G genotypes (G1-G5, G9,
        G11) and 3 P genotypes (P[6], P[7], P[13]) in North America, Europe,
        and Asia.  At least 25 sequences are required to represent the diversity
        needed for a pan-porcine-RVA conserved VP6 window.
        """
        assert rota_cfg.max_sequences >= 25, (
            "Porcine Rotavirus A has diverse G and P genotypes; at least 25 "
            "sequences are required to identify a conserved pan-porcine VP6 window"
        )

    def test_amplicon_size_within_lamp_window(self, rota_cfg: PipelineConfig) -> None:
        """Amplicon size must stay within the LAMP reaction geometry."""
        assert 80 <= rota_cfg.f2_b2_min < rota_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, rota_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 outer primers flanking F2-B2."""
        assert rota_cfg.min_region_length >= rota_cfg.f2_b2_max + 40, (
            "min_region_length must exceed f2_b2_max by at least 40 nt "
            "to accommodate F3 and B3 outer primer binding sites"
        )

    def test_output_dir_references_rota(self, rota_cfg: PipelineConfig) -> None:
        """Output directory must be identifiable as a Porcine Rotavirus A run."""
        assert "rota" in str(rota_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, rota_cfg: PipelineConfig) -> None:
        assert 0.70 <= rota_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, rota_cfg: PipelineConfig) -> None:
        assert 0.70 <= rota_cfg.min_coverage_threshold <= 0.95

    def test_tm_range_spans_rt_lamp_window(self, rota_cfg: PipelineConfig) -> None:
        """Tm window must span >= 1 degC to allow primer3 design flexibility."""
        assert rota_cfg.tm_max - rota_cfg.tm_min >= 1.0, (
            "tm_max - tm_min must be >= 1.0 degC to allow primer3 design flexibility"
        )
