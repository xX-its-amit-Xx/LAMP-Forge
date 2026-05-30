"""Validate that every bundled example config in config/ parses correctly.

These tests catch schema drift: if load_config's validation rules change, a
bundled YAML that no longer satisfies the new rules will fail here before a
user discovers the mismatch at runtime.  All tests are pure Python — no
external binaries (MAFFT, BLAST) are required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_YAML_CONFIGS = sorted(CONFIG_DIR.glob("*.yaml"))


@pytest.mark.parametrize(
    "config_path",
    _YAML_CONFIGS,
    ids=[p.stem for p in _YAML_CONFIGS],
)
def test_example_config_parses(
    config_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each bundled YAML must load to a PipelineConfig without raising ConfigError."""
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


class TestFmdvConfig:
    """Specific semantic checks for the FMDV RT-LAMP config."""

    @pytest.fixture
    def fmdv_cfg(self, monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
        monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
        return load_config(CONFIG_DIR / "fmdv_3dpol.yaml")

    def test_target_name(self, fmdv_cfg: PipelineConfig) -> None:
        assert fmdv_cfg.target_name == "fmdv_3Dpol"

    def test_taxon_id_is_fmdv(self, fmdv_cfg: PipelineConfig) -> None:
        assert fmdv_cfg.taxon_id == 12110

    def test_gene_is_3d(self, fmdv_cfg: PipelineConfig) -> None:
        assert fmdv_cfg.gene == "3D"

    def test_rt_lamp_tm_floor(self, fmdv_cfg: PipelineConfig) -> None:
        """RT-LAMP one-step protocol requires tm_min >= 63.0 degC."""
        assert fmdv_cfg.tm_min >= 63.0, (
            "FMDV is an RNA target: tm_min must be >= 63.0 for one-step RT-LAMP co-activity"
        )

    def test_max_sequences_spans_all_serotypes(self, fmdv_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover all 7 FMDV serotypes."""
        assert fmdv_cfg.max_sequences >= 35

    def test_specificity_thresholds_are_tight(self, fmdv_cfg: PipelineConfig) -> None:
        """SVD (swine vesicular disease) is a picornavirus; use tight thresholds."""
        assert fmdv_cfg.min_identity_threshold >= 0.85
        assert fmdv_cfg.min_coverage_threshold >= 0.85

    def test_conservation_entropy_allows_rna_virus_drift(
        self, fmdv_cfg: PipelineConfig
    ) -> None:
        """RNA viruses vary more than bacteria; entropy threshold should be >= 0.25."""
        assert fmdv_cfg.entropy_threshold >= 0.25

    def test_amplicon_size_within_lamp_window(self, fmdv_cfg: PipelineConfig) -> None:
        assert 80 <= fmdv_cfg.f2_b2_min < fmdv_cfg.f2_b2_max <= 200

    def test_output_dir_is_fmdv(self, fmdv_cfg: PipelineConfig) -> None:
        assert "fmdv" in str(fmdv_cfg.output_dir).lower()
