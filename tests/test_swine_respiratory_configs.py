"""Validate the two new swine respiratory PRDC LAMP configs.

The swine respiratory / PRDC panel (``lamp-forge swine-respiratory-risk``)
monitors three pathogens.  PRRSV already has its own config tested elsewhere
(config/prrsv_orf7.yaml); this module validates the two new configs that
complete the three-channel PRDC panel:

    Config                       Target                    NA type  Gene / marker
    pcv2_cap.yaml                Porcine circovirus 2      DNA      cap (ORF2)
    mycoplasma_hyop_p97.yaml     Mycoplasma hyopneumoniae  DNA      p97 adhesin

Both are standard LAMP (DNA targets; no reverse transcription required).
Together with PRRSV (RT-LAMP; prrsv_orf7.yaml) they form the complete
three-channel PRDC panel for BioVind portable deployment.

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"

_PCV2_CONFIG = CONFIG_DIR / "pcv2_cap.yaml"
_MHP_CONFIG = CONFIG_DIR / "mycoplasma_hyop_p97.yaml"

_STANDARD_LAMP_FLOOR = 60.0  # degC; no RT enzyme required for DNA targets


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pcv2_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the PCV2 cap-gene config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_PCV2_CONFIG)


@pytest.fixture
def mhp_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the Mycoplasma hyopneumoniae p97 config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_MHP_CONFIG)


# ---------------------------------------------------------------------------
# Smoke tests -- both configs must parse without error
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "config_path",
    [_PCV2_CONFIG, _MHP_CONFIG],
    ids=["pcv2_cap", "mycoplasma_hyop_p97"],
)
def test_swine_respiratory_config_parses(
    config_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each PRDC config must load to a PipelineConfig without raising ConfigError."""
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
# Target taxonomy IDs
# ---------------------------------------------------------------------------


class TestTaxonIds:
    def test_pcv2_taxon_id(self, pcv2_cfg: PipelineConfig) -> None:
        assert pcv2_cfg.taxon_id == 85413  # Porcine circovirus 2

    def test_mhp_taxon_id(self, mhp_cfg: PipelineConfig) -> None:
        assert mhp_cfg.taxon_id == 2104  # Mycoplasma hyopneumoniae


# ---------------------------------------------------------------------------
# Target gene names
# ---------------------------------------------------------------------------


class TestGeneNames:
    def test_pcv2_gene(self, pcv2_cfg: PipelineConfig) -> None:
        assert pcv2_cfg.gene == "cap"

    def test_mhp_gene(self, mhp_cfg: PipelineConfig) -> None:
        assert mhp_cfg.gene == "p97"


# ---------------------------------------------------------------------------
# Standard LAMP temperature windows (DNA targets -- no RT-LAMP needed)
# ---------------------------------------------------------------------------


class TestPrimerTmRanges:
    def test_pcv2_tm_min_is_standard_lamp(self, pcv2_cfg: PipelineConfig) -> None:
        assert pcv2_cfg.tm_min >= _STANDARD_LAMP_FLOOR

    def test_mhp_tm_min_is_standard_lamp(self, mhp_cfg: PipelineConfig) -> None:
        assert mhp_cfg.tm_min >= _STANDARD_LAMP_FLOOR

    def test_pcv2_tm_max_above_tm_min(self, pcv2_cfg: PipelineConfig) -> None:
        assert pcv2_cfg.tm_max > pcv2_cfg.tm_min

    def test_mhp_tm_max_above_tm_min(self, mhp_cfg: PipelineConfig) -> None:
        assert mhp_cfg.tm_max > mhp_cfg.tm_min


# ---------------------------------------------------------------------------
# GC content windows
# ---------------------------------------------------------------------------


class TestGcWindows:
    def test_pcv2_gc_window_valid(self, pcv2_cfg: PipelineConfig) -> None:
        # PCV2 cap gene is GC-rich; gc_min should be >= 40
        assert pcv2_cfg.gc_min >= 40.0
        assert pcv2_cfg.gc_max <= 70.0
        assert pcv2_cfg.gc_max > pcv2_cfg.gc_min

    def test_mhp_gc_window_valid(self, mhp_cfg: PipelineConfig) -> None:
        # Mhp genome is AT-rich (~28% GC); gc_min should accommodate this
        assert mhp_cfg.gc_min < 40.0
        assert mhp_cfg.gc_max <= 65.0
        assert mhp_cfg.gc_max > mhp_cfg.gc_min


# ---------------------------------------------------------------------------
# Amplicon size bounds
# ---------------------------------------------------------------------------


class TestAmpliconSize:
    def test_pcv2_amplicon_size(self, pcv2_cfg: PipelineConfig) -> None:
        assert pcv2_cfg.f2_b2_min >= 100
        assert pcv2_cfg.f2_b2_max <= 200
        assert pcv2_cfg.f2_b2_max > pcv2_cfg.f2_b2_min

    def test_mhp_amplicon_size(self, mhp_cfg: PipelineConfig) -> None:
        assert mhp_cfg.f2_b2_min >= 100
        assert mhp_cfg.f2_b2_max <= 200
        assert mhp_cfg.f2_b2_max > mhp_cfg.f2_b2_min


# ---------------------------------------------------------------------------
# Conservation parameters
# ---------------------------------------------------------------------------


class TestConservationParams:
    def test_pcv2_entropy_threshold_allows_genotype_diversity(
        self, pcv2_cfg: PipelineConfig
    ) -> None:
        # PCV2 genotypes vary ~15-25%; threshold should be >= 0.25 to be
        # practical (tighter thresholds would find zero conserved windows)
        assert pcv2_cfg.entropy_threshold >= 0.25

    def test_mhp_min_region_length(self, mhp_cfg: PipelineConfig) -> None:
        assert mhp_cfg.min_region_length >= 200

    def test_pcv2_min_region_length(self, pcv2_cfg: PipelineConfig) -> None:
        assert pcv2_cfg.min_region_length >= 200
