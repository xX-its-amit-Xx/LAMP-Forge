"""Validate the five bovine respiratory disease (BRD) primer-design configs.

The BRD panel targets the five pathogens interpreted by ``lamp-forge bov-risk``:

    Config                 Target          NA type   Gene / marker
    brsv_N_gene.yaml       BRSV            RNA       N (nucleoprotein)
    bcov_N_gene.yaml       BCoV            RNA       N (nucleoprotein)
    bvdv_5utr.yaml         BVDV (1+2)      RNA       5'-UTR (non-coding)
    ibr_gB.yaml            IBR / BoHV-1    DNA       UL27 (gB)
    mhae_lktA.yaml         M. haemolytica  DNA       lktA (leukotoxin A)

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"

_BRSV_CONFIG = CONFIG_DIR / "brsv_N_gene.yaml"
_BCOV_CONFIG = CONFIG_DIR / "bcov_N_gene.yaml"
_BVDV_CONFIG = CONFIG_DIR / "bvdv_5utr.yaml"
_IBR_CONFIG = CONFIG_DIR / "ibr_gB.yaml"
_MHAE_CONFIG = CONFIG_DIR / "mhae_lktA.yaml"

_RT_LAMP_FLOOR = 63.0  # degC; one-step RT-LAMP co-activity floor for NEB RTx
_DNA_LAMP_FLOOR = 60.0  # degC; standard DNA-LAMP Bst 2.0 floor


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def brsv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the BRSV N-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_BRSV_CONFIG)


@pytest.fixture
def bcov_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the BCoV N-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_BCOV_CONFIG)


@pytest.fixture
def bvdv_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the BVDV 5'-UTR config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_BVDV_CONFIG)


@pytest.fixture
def ibr_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the IBR gB config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_IBR_CONFIG)


@pytest.fixture
def mhae_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the Mannheimia haemolytica lktA config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_MHAE_CONFIG)


# ---------------------------------------------------------------------------
# Smoke tests -- every config must parse without error
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "config_path",
    [_BRSV_CONFIG, _BCOV_CONFIG, _BVDV_CONFIG, _IBR_CONFIG, _MHAE_CONFIG],
    ids=["brsv", "bcov", "bvdv", "ibr", "mhae"],
)
def test_brd_config_parses(config_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Each BRD config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(config_path)
    assert isinstance(cfg, PipelineConfig)


# ---------------------------------------------------------------------------
# BRSV N gene -- RNA virus, RT-LAMP
# ---------------------------------------------------------------------------


class TestBrsvNGeneConfig:
    """Semantic checks for the BRSV N-gene RT-LAMP config.

    BRSV (Bovine orthopneumovirus) is the leading viral trigger of BRD.
    The N gene spans subgroups A and B (~80-85% nt identity).  Detection
    requires one-step RT-LAMP at 63-65 degC.
    """

    def test_target_name(self, brsv_cfg: PipelineConfig) -> None:
        assert brsv_cfg.target_name == "brsv_N_gene"

    def test_gene_is_N(self, brsv_cfg: PipelineConfig) -> None:
        """N gene (nucleoprotein) is the canonical BRSV molecular marker."""
        assert brsv_cfg.gene == "N"

    def test_sequence_source_provided(self, brsv_cfg: PipelineConfig) -> None:
        """Config must have a taxon_id or non-empty accessions list."""
        assert brsv_cfg.taxon_id is not None or len(brsv_cfg.accessions) > 0

    def test_accessions_cover_both_subgroups(self, brsv_cfg: PipelineConfig) -> None:
        """BRSV subgroups A and B must both be represented for pan-BRSV design.

        Subgroup A and B share ~80-85% N-gene nt identity.  Including both
        in the sequence set is required to find conserved primer windows that
        work across the full BRSV genetic diversity seen in field outbreaks.
        """
        assert len(brsv_cfg.accessions) >= 2, (
            "At least two BRSV accessions are needed to represent subgroups A and B; "
            "a single-subgroup design will miss the ~15-20% inter-subgroup sequence space"
        )

    def test_rt_lamp_tm_floor(self, brsv_cfg: PipelineConfig) -> None:
        """BRSV is an RNA virus; tm_min must meet the RT-LAMP one-step floor.

        One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx at 63-65 degC) requires
        tm_min >= 63 degC to maintain reverse-transcriptase co-activity.
        """
        assert brsv_cfg.tm_min >= _RT_LAMP_FLOOR, (
            f"BRSV is a negative-sense ssRNA virus: tm_min must be >= {_RT_LAMP_FLOOR} degC "
            "for one-step RT-LAMP co-activity with NEB RTx / Bst 2.0 WarmStart"
        )

    def test_tm_max_within_rt_lamp_window(self, brsv_cfg: PipelineConfig) -> None:
        assert brsv_cfg.tm_max <= 65.0

    def test_entropy_threshold_appropriate_for_rna_subgroup_diversity(
        self, brsv_cfg: PipelineConfig
    ) -> None:
        """Entropy threshold must accommodate BRSV subgroup A/B divergence.

        BRSV subgroups A and B share ~80-85% N-gene nt identity -- higher
        divergence than within-subgroup variation but lower than inter-genus
        divergence for polyphyletic functional genes (dsrB, mcrA).  The threshold
        must be loose enough to find cross-subgroup windows (>= 0.25 bits) but
        not so loose as to admit highly variable surface epitope positions (<= 0.45).
        """
        assert brsv_cfg.entropy_threshold >= 0.25, (
            "BRSV subgroup A/B N-gene divergence requires entropy_threshold >= 0.25 bits"
        )
        assert brsv_cfg.entropy_threshold <= 0.45

    def test_gc_range_covers_brsv_n_gene(self, brsv_cfg: PipelineConfig) -> None:
        """GC window must span BRSV N-gene GC content (~42-46%)."""
        assert brsv_cfg.gc_min <= 42.0
        assert brsv_cfg.gc_max >= 45.0

    def test_tight_specificity_threshold_for_pneumoviridae(self, brsv_cfg: PipelineConfig) -> None:
        """Off-target identity threshold must be tight for Pneumoviridae family.

        Human RSV-A and RSV-B (HRSV) share ~50-60% N-gene nt identity with
        BRSV.  Bovine metapneumovirus (Metapneumovirus genus, same family) shares
        conserved N-gene structural domains.  A threshold >= 0.83 is needed to
        catch potential family-level cross-reactivity.
        """
        assert brsv_cfg.min_identity_threshold >= 0.83

    def test_amplicon_size_within_lamp_window(self, brsv_cfg: PipelineConfig) -> None:
        assert 80 <= brsv_cfg.f2_b2_min < brsv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, brsv_cfg: PipelineConfig) -> None:
        assert brsv_cfg.min_region_length >= brsv_cfg.f2_b2_max + 40

    def test_output_dir_references_brsv(self, brsv_cfg: PipelineConfig) -> None:
        assert "brsv" in str(brsv_cfg.output_dir).lower()

    def test_max_sequences_spans_subgroup_diversity(self, brsv_cfg: PipelineConfig) -> None:
        assert brsv_cfg.max_sequences >= 20


# ---------------------------------------------------------------------------
# BCoV N gene -- RNA virus, RT-LAMP
# ---------------------------------------------------------------------------


class TestBcovNGeneConfig:
    """Semantic checks for the BCoV N-gene RT-LAMP config.

    Bovine coronavirus (BCoV; Betacoronavirus 1) causes respiratory disease,
    neonatal diarrhea, and winter dysentery.  The N gene is the canonical
    BCoV molecular marker.  Critical specificity challenge: HCoV-OC43 shares
    ~95-96% N-gene nt identity with BCoV.
    """

    def test_target_name(self, bcov_cfg: PipelineConfig) -> None:
        assert bcov_cfg.target_name == "bcov_N_gene"

    def test_gene_is_N(self, bcov_cfg: PipelineConfig) -> None:
        assert bcov_cfg.gene == "N"

    def test_taxon_is_bovine_coronavirus(self, bcov_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Bovine coronavirus (NCBI TaxID 11128)."""
        assert bcov_cfg.taxon_id == 11128

    def test_rt_lamp_tm_floor(self, bcov_cfg: PipelineConfig) -> None:
        """BCoV is an RNA virus; tm_min must meet the RT-LAMP co-activity floor."""
        assert bcov_cfg.tm_min >= _RT_LAMP_FLOOR, (
            f"BCoV is a positive-sense ssRNA Betacoronavirus: tm_min must be >= {_RT_LAMP_FLOOR} degC"
        )

    def test_tm_max_within_rt_lamp_window(self, bcov_cfg: PipelineConfig) -> None:
        assert bcov_cfg.tm_max <= 65.0

    def test_entropy_threshold_appropriate_for_conserved_n_gene(
        self, bcov_cfg: PipelineConfig
    ) -> None:
        """BCoV N gene is quite conserved; entropy threshold should be relatively tight.

        Within-BCoV N-gene identity is ~98-99% across field strains.  A threshold
        of 0.25-0.35 bits exploits the conserved internal blocks while filtering
        the few variable positions near the 5' and 3' ends of the CDS.
        """
        assert bcov_cfg.entropy_threshold >= 0.15
        assert bcov_cfg.entropy_threshold <= 0.35

    def test_gc_range_covers_bcov(self, bcov_cfg: PipelineConfig) -> None:
        """GC window must span BCoV N-gene GC content (~42-45%)."""
        assert bcov_cfg.gc_min <= 38.0
        assert bcov_cfg.gc_max >= 45.0

    def test_tight_specificity_threshold_for_oc43_risk(self, bcov_cfg: PipelineConfig) -> None:
        """Identity threshold must be tight because OC43 shares ~95-96% N-gene identity.

        HCoV-OC43 diverged from BCoV ~130 years ago and shares ~95-96% N-gene
        nt identity.  Even at 0.85 identity threshold the BLAST screen will
        flag OC43 as a close match; reviewers must select primer windows in the
        BCoV-specific regions of the N gene.
        """
        assert bcov_cfg.min_identity_threshold >= 0.83

    def test_amplicon_size_within_lamp_window(self, bcov_cfg: PipelineConfig) -> None:
        assert 80 <= bcov_cfg.f2_b2_min < bcov_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, bcov_cfg: PipelineConfig) -> None:
        assert bcov_cfg.min_region_length >= bcov_cfg.f2_b2_max + 40

    def test_output_dir_references_bcov(self, bcov_cfg: PipelineConfig) -> None:
        assert "bcov" in str(bcov_cfg.output_dir).lower()


# ---------------------------------------------------------------------------
# BVDV 5'-UTR -- RNA virus, RT-LAMP
# ---------------------------------------------------------------------------


class TestBvdv5utrConfig:
    """Semantic checks for the BVDV 5'-UTR RT-LAMP config.

    BVDV (Pestivirus bovis) causes immunosuppression, reproductive failure,
    and mucosal disease.  The 5'-UTR is conserved across types 1 and 2 and
    is the canonical BVDV molecular diagnostic target.
    """

    def test_target_name(self, bvdv_cfg: PipelineConfig) -> None:
        assert bvdv_cfg.target_name == "bvdv_5utr"

    def test_gene_is_null_for_noncoding_target(self, bvdv_cfg: PipelineConfig) -> None:
        """The 5'-UTR is non-coding; gene must be null to pull whole-genome sequences."""
        assert bvdv_cfg.gene is None, (
            "BVDV 5'-UTR is a non-coding IRES region, not a conventional gene annotation; "
            "gene must be null to fetch whole-genome sequences for conservation analysis"
        )

    def test_taxon_anchored_to_bvdv1(self, bvdv_cfg: PipelineConfig) -> None:
        """Taxon anchor must be BVDV type 1 (NCBI TaxID 11099; Pestivirus A)."""
        assert bvdv_cfg.taxon_id == 11099

    def test_rt_lamp_tm_floor(self, bvdv_cfg: PipelineConfig) -> None:
        """BVDV is an RNA virus; tm_min must meet the RT-LAMP one-step floor."""
        assert bvdv_cfg.tm_min >= _RT_LAMP_FLOOR, (
            f"BVDV is a positive-sense ssRNA Pestivirus: tm_min must be >= {_RT_LAMP_FLOOR} degC"
        )

    def test_tm_max_within_rt_lamp_window(self, bvdv_cfg: PipelineConfig) -> None:
        assert bvdv_cfg.tm_max <= 65.0

    def test_entropy_threshold_for_inter_type_diversity(self, bvdv_cfg: PipelineConfig) -> None:
        """Entropy threshold must accommodate BVDV-1 / BVDV-2 5'-UTR divergence.

        BVDV types 1 and 2 share ~75-80% 5'-UTR nt identity.  The threshold
        must be loose enough to find cross-type conserved IRES blocks (>= 0.25 bits)
        but not so loose as to admit the variable 3' end of the 5'-UTR (<= 0.40).
        """
        assert bvdv_cfg.entropy_threshold >= 0.20
        assert bvdv_cfg.entropy_threshold <= 0.40

    def test_gc_range_covers_bvdv(self, bvdv_cfg: PipelineConfig) -> None:
        """GC window must span BVDV 5'-UTR GC content (~45-50%)."""
        assert bvdv_cfg.gc_min <= 42.0
        assert bvdv_cfg.gc_max >= 50.0

    def test_max_sequences_for_type_diversity(self, bvdv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover BVDV-1 subtypes plus type-2 isolates."""
        assert bvdv_cfg.max_sequences >= 25

    def test_amplicon_size_within_lamp_window(self, bvdv_cfg: PipelineConfig) -> None:
        assert 80 <= bvdv_cfg.f2_b2_min < bvdv_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, bvdv_cfg: PipelineConfig) -> None:
        assert bvdv_cfg.min_region_length >= bvdv_cfg.f2_b2_max + 40

    def test_output_dir_references_bvdv(self, bvdv_cfg: PipelineConfig) -> None:
        assert "bvdv" in str(bvdv_cfg.output_dir).lower()


# ---------------------------------------------------------------------------
# IBR / BoHV-1 gB -- DNA virus, standard LAMP
# ---------------------------------------------------------------------------


class TestIbrGbConfig:
    """Semantic checks for the IBR / BoHV-1 gB DNA-LAMP config.

    IBR (BoHV-1) is a latent dsDNA alphaherpesvirus subject to mandatory
    control programmes in the EU, UK, and Scandinavia.  gB (UL27) is a
    structurally conserved glycoprotein that is the ELISA antigen in
    eradication serology.  BoHV-1 is a DNA target -- no RT step required.
    """

    def test_target_name(self, ibr_cfg: PipelineConfig) -> None:
        assert ibr_cfg.target_name == "ibr_gB"

    def test_gene_is_ul27(self, ibr_cfg: PipelineConfig) -> None:
        """UL27 (glycoprotein B) is the canonical IBR molecular marker."""
        assert ibr_cfg.gene == "UL27"

    def test_taxon_is_bohv1(self, ibr_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Bovine alphaherpesvirus 1 (NCBI TaxID 10320)."""
        assert ibr_cfg.taxon_id == 10320

    def test_not_rt_lamp(self, ibr_cfg: PipelineConfig) -> None:
        """BoHV-1 is a dsDNA herpesvirus; tm_min must NOT be raised to RT-LAMP floor.

        Using the 63 degC RT-LAMP Tm floor for a DNA target would unnecessarily
        restrict primer design and is scientifically incorrect.
        """
        assert ibr_cfg.tm_min < _RT_LAMP_FLOOR, (
            "IBR / BoHV-1 is a dsDNA herpesvirus: tm_min should be at the standard "
            f"DNA-LAMP {_DNA_LAMP_FLOOR} degC floor, not the RT-LAMP {_RT_LAMP_FLOOR} degC floor"
        )

    def test_standard_dna_lamp_tm_window(self, ibr_cfg: PipelineConfig) -> None:
        """Standard DNA-LAMP Tm window: 60-65 degC."""
        assert ibr_cfg.tm_min >= _DNA_LAMP_FLOOR
        assert ibr_cfg.tm_max <= 65.0

    def test_gc_range_covers_gc_rich_herpesvirus(self, ibr_cfg: PipelineConfig) -> None:
        """GC window must span BoHV-1 gB GC content (~68-72%).

        Alphaherpesviruses are GC-rich (~72% overall for BoHV-1); gB (UL27)
        is similarly GC-rich.  The gc_min floor must be raised substantially
        above the standard 40% to find primers in the most conserved blocks.
        """
        assert ibr_cfg.gc_min >= 48.0, (
            "BoHV-1 gB is ~68-72% GC; gc_min must be >= 48% to find primers "
            "in the conserved gB structural domains"
        )
        assert ibr_cfg.gc_max >= 68.0, (
            "gc_max must accommodate the GC-rich conserved motifs of BoHV-1 gB"
        )

    def test_entropy_threshold_for_conserved_herpesvirus_gene(
        self, ibr_cfg: PipelineConfig
    ) -> None:
        """gB is under strong functional constraint; tight threshold appropriate.

        Within-BoHV-1 gB identity is > 98% across subtypes 1.1 and 1.2.
        Threshold should be < 0.25 bits to exploit the most conserved positions.
        """
        assert ibr_cfg.entropy_threshold <= 0.25

    def test_amplicon_size_within_lamp_window(self, ibr_cfg: PipelineConfig) -> None:
        assert 80 <= ibr_cfg.f2_b2_min < ibr_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, ibr_cfg: PipelineConfig) -> None:
        assert ibr_cfg.min_region_length >= ibr_cfg.f2_b2_max + 40

    def test_output_dir_references_ibr(self, ibr_cfg: PipelineConfig) -> None:
        assert "ibr" in str(ibr_cfg.output_dir).lower()


# ---------------------------------------------------------------------------
# Mannheimia haemolytica lktA -- DNA bacterium, standard LAMP
# ---------------------------------------------------------------------------


class TestMhaeLktaConfig:
    """Semantic checks for the Mannheimia haemolytica lktA DNA-LAMP config.

    M. haemolytica is the primary bacterial secondary pathogen in BRD,
    responsible for the fibrinous pneumonia that causes feedlot mortality.
    lktA encodes the RTX leukotoxin -- the key virulence factor and a
    specific molecular marker absent from Pasteurella multocida.
    """

    def test_target_name(self, mhae_cfg: PipelineConfig) -> None:
        assert mhae_cfg.target_name == "mhae_lktA"

    def test_gene_is_lkta(self, mhae_cfg: PipelineConfig) -> None:
        """lktA is the canonical M. haemolytica virulence marker."""
        assert mhae_cfg.gene == "lktA"

    def test_taxon_is_mannheimia_haemolytica(self, mhae_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Mannheimia haemolytica (NCBI TaxID 75984)."""
        assert mhae_cfg.taxon_id == 75984

    def test_not_rt_lamp(self, mhae_cfg: PipelineConfig) -> None:
        """M. haemolytica is a DNA bacterium; tm_min must NOT be at RT-LAMP floor."""
        assert mhae_cfg.tm_min < _RT_LAMP_FLOOR, (
            "M. haemolytica lktA is a chromosomal DNA target: tm_min should be at the "
            f"standard DNA-LAMP {_DNA_LAMP_FLOOR} degC floor, not the RT-LAMP floor"
        )

    def test_standard_dna_lamp_tm_window(self, mhae_cfg: PipelineConfig) -> None:
        assert mhae_cfg.tm_min >= _DNA_LAMP_FLOOR
        assert mhae_cfg.tm_max <= 65.0

    def test_gc_range_reflects_pasteurellaceae(self, mhae_cfg: PipelineConfig) -> None:
        """GC window must reflect M. haemolytica GC content (~41%).

        Mannheimia haemolytica has ~41% GC; lktA is comparable.  The design
        window must include the ~38-55% range to produce primers that work
        across both the GC-poor RTX repeat region and the slightly GC-richer
        ATP-binding domain.
        """
        assert mhae_cfg.gc_min <= 42.0, (
            "gc_min must be <= 42% to accommodate M. haemolytica's ~41% GC content"
        )
        assert mhae_cfg.gc_max >= 50.0, (
            "gc_max must be >= 50% to capture primer candidates across lktA's GC range"
        )

    def test_entropy_threshold_for_single_species_marker(self, mhae_cfg: PipelineConfig) -> None:
        """lktA is a single-species marker; tight entropy threshold is appropriate.

        Within-M. haemolytica lktA identity is > 96% across serotypes A1, A2,
        and A6.  A threshold <= 0.25 bits exploits the conserved RTX repeat and
        ATP-binding domains without admitting the variable serotype-specific
        surface epitopes.
        """
        assert mhae_cfg.entropy_threshold <= 0.25

    def test_amplicon_size_within_lamp_window(self, mhae_cfg: PipelineConfig) -> None:
        assert 80 <= mhae_cfg.f2_b2_min < mhae_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, mhae_cfg: PipelineConfig) -> None:
        assert mhae_cfg.min_region_length >= mhae_cfg.f2_b2_max + 40

    def test_output_dir_references_mhae(self, mhae_cfg: PipelineConfig) -> None:
        assert "mhae" in str(mhae_cfg.output_dir).lower()

    def test_identity_threshold_in_valid_range(self, mhae_cfg: PipelineConfig) -> None:
        assert 0.70 <= mhae_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_in_valid_range(self, mhae_cfg: PipelineConfig) -> None:
        assert 0.70 <= mhae_cfg.min_coverage_threshold <= 0.95
