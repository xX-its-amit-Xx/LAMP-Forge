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

    def test_conservation_entropy_allows_rna_virus_drift(self, fmdv_cfg: PipelineConfig) -> None:
        """RNA viruses vary more than bacteria; entropy threshold should be >= 0.25."""
        assert fmdv_cfg.entropy_threshold >= 0.25

    def test_amplicon_size_within_lamp_window(self, fmdv_cfg: PipelineConfig) -> None:
        assert 80 <= fmdv_cfg.f2_b2_min < fmdv_cfg.f2_b2_max <= 200

    def test_output_dir_is_fmdv(self, fmdv_cfg: PipelineConfig) -> None:
        assert "fmdv" in str(fmdv_cfg.output_dir).lower()


class TestSrbDsrBConfig:
    """Semantic checks for the SRB dsrB oilfield souring config.

    dsrB is a polyphyletic functional gene (present across five phyla of
    sulfate-reducers) so its config must tolerate higher sequence entropy
    than single-genus housekeeping-gene targets.
    """

    @pytest.fixture
    def srb_cfg(self, monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
        monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
        return load_config(CONFIG_DIR / "srb_dsrB.yaml")

    def test_target_name(self, srb_cfg: PipelineConfig) -> None:
        assert srb_cfg.target_name == "srb_dsrB_oilfield"

    def test_gene_is_dsrb(self, srb_cfg: PipelineConfig) -> None:
        assert srb_cfg.gene == "dsrB"

    def test_taxon_anchored_to_desulfovibrio(self, srb_cfg: PipelineConfig) -> None:
        """Anchor genus for broad SRB design is Desulfovibrio (taxon 872)."""
        assert srb_cfg.taxon_id == 872

    def test_entropy_threshold_allows_polyphyletic_diversity(self, srb_cfg: PipelineConfig) -> None:
        """dsrB spans five phyla; entropy threshold must be >= 0.30 bits."""
        assert srb_cfg.entropy_threshold >= 0.30, (
            "dsrB is functionally but not tightly sequence-conserved across SRB phyla; "
            "entropy_threshold must be >= 0.30 to yield conserved windows"
        )

    def test_max_sequences_samples_srb_diversity(self, srb_cfg: PipelineConfig) -> None:
        """Broad SRB diversity requires a large sequence set."""
        assert srb_cfg.max_sequences >= 30

    def test_gc_range_covers_proteobacteria_and_firmicutes(self, srb_cfg: PipelineConfig) -> None:
        """GC range must be broad enough for SRB genera spanning ~45-60% GC."""
        assert srb_cfg.gc_min <= 40.0
        assert srb_cfg.gc_max >= 60.0

    def test_standard_lamp_tm_window(self, srb_cfg: PipelineConfig) -> None:
        """SRB targets are DNA; standard 60-65 degC LAMP window is correct."""
        assert srb_cfg.tm_min >= 60.0
        assert srb_cfg.tm_max <= 65.0

    def test_output_dir_references_srb(self, srb_cfg: PipelineConfig) -> None:
        assert "srb" in str(srb_cfg.output_dir).lower()


class TestMethanogenMcrAConfig:
    """Semantic checks for the methanogen mcrA oilfield souring config.

    mcrA is the universal methanogen marker spanning all five archaeal
    methanogen orders.  The config must tolerate archaeal GC-range diversity
    and relaxed entropy thresholds compared to bacterial housekeeping genes.
    """

    @pytest.fixture
    def mcra_cfg(self, monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
        monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
        return load_config(CONFIG_DIR / "methanogen_mcrA.yaml")

    def test_target_name(self, mcra_cfg: PipelineConfig) -> None:
        assert mcra_cfg.target_name == "methanogen_mcrA"

    def test_gene_is_mcra(self, mcra_cfg: PipelineConfig) -> None:
        assert mcra_cfg.gene == "mcrA"

    def test_taxon_is_archaea_domain(self, mcra_cfg: PipelineConfig) -> None:
        """Anchor taxon is Archaea domain (NCBI taxon 2157) for broad coverage."""
        assert mcra_cfg.taxon_id == 2157

    def test_entropy_threshold_allows_order_level_diversity(self, mcra_cfg: PipelineConfig) -> None:
        """mcrA spans five methanogen orders; entropy threshold must be >= 0.30."""
        assert mcra_cfg.entropy_threshold >= 0.30, (
            "mcrA is under catalytic constraint but diverges across archaeal orders; "
            "entropy_threshold must be >= 0.30 to capture conserved catalytic windows"
        )

    def test_gc_range_spans_archaeal_diversity(self, mcra_cfg: PipelineConfig) -> None:
        """Archaea span ~30-65% GC; the design window must be at least 25 units wide."""
        assert mcra_cfg.gc_max - mcra_cfg.gc_min >= 25.0, (
            "Archaeal GC content varies widely; gc window must be wide enough to find primers"
        )

    def test_gc_min_allows_at_rich_archaea(self, mcra_cfg: PipelineConfig) -> None:
        """AT-rich methanogens (e.g. Methanococcus jannaschii, 31% GC) need gc_min <= 35."""
        assert mcra_cfg.gc_min <= 35.0

    def test_max_sequences_spans_five_orders(self, mcra_cfg: PipelineConfig) -> None:
        """Five methanogen orders require a broad sequence set."""
        assert mcra_cfg.max_sequences >= 30

    def test_standard_lamp_tm_window(self, mcra_cfg: PipelineConfig) -> None:
        """Methanogen targets are DNA; standard 60-65 degC window is correct."""
        assert mcra_cfg.tm_min >= 60.0
        assert mcra_cfg.tm_max <= 65.0

    def test_output_dir_references_methanogen(self, mcra_cfg: PipelineConfig) -> None:
        assert "methanogen" in str(mcra_cfg.output_dir).lower()
