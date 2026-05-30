"""Validate the NRB narG oilfield souring management config (config/nrb_narG.yaml).

Semantic checks specific to nitrate-reducing bacteria (NRB) via narG:
- Paracoccus genus anchor (NCBI taxon 265)
- Entropy threshold appropriate for narG (polyphyletic functional gene)
- GC window covers GC-rich NRB genera (Paracoccus ~67%, Pseudomonas ~63%)
- Standard LAMP Tm window, not RT-LAMP (narG is a chromosomal DNA target)

All tests are pure Python -- no external binaries (MAFFT, BLAST) required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.config import load_config
from lamp_forge.types import PipelineConfig

CONFIG_DIR = Path(__file__).parent.parent / "config"
_NRB_CONFIG = CONFIG_DIR / "nrb_narG.yaml"


@pytest.fixture
def nrb_cfg(monkeypatch: pytest.MonkeyPatch) -> PipelineConfig:
    """Load the NRB narG config with a CI-safe NCBI email."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    return load_config(_NRB_CONFIG)


def test_nrb_config_parses(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config must load to a PipelineConfig without raising ConfigError."""
    monkeypatch.setenv("NCBI_EMAIL", "ci@lamp-forge.example")
    cfg = load_config(_NRB_CONFIG)
    assert isinstance(cfg, PipelineConfig)


class TestNrbNarGConfig:
    """Semantic checks for the NRB narG oilfield souring management config.

    narG (nitrate reductase alpha subunit) is the canonical functional-gene
    marker for nitrate-reducing bacteria in environmental monitoring.  In
    oilfields, NRB are the target organism for nitrate injection biocide
    programs: injected nitrate stimulates NRB to outcompete SRB for electron
    donors and produce nitrite, a direct inhibitor of sulfate reductase.
    Monitoring NRB:SRB balance via narG:dsrB provides a real-time efficacy
    readout for souring-control programs.
    """

    def test_target_name(self, nrb_cfg: PipelineConfig) -> None:
        assert nrb_cfg.target_name == "nrb_narG_oilfield"

    def test_gene_is_narG(self, nrb_cfg: PipelineConfig) -> None:
        assert nrb_cfg.gene == "narG"

    def test_taxon_is_paracoccus_genus(self, nrb_cfg: PipelineConfig) -> None:
        """Anchor taxon must be Paracoccus genus (NCBI taxon 265).

        Paracoccus is the model respiratory denitrifier genus and primary
        narG reference clade.  The config supplements this anchor with
        accessions from key oilfield NRB genera (Pseudomonas stutzeri,
        Thiobacillus denitrificans, Thauera aromatica) for field coverage.
        """
        assert nrb_cfg.taxon_id == 265

    def test_entropy_threshold_reflects_polyphyletic_marker(self, nrb_cfg: PipelineConfig) -> None:
        """narG is a polyphyletic functional gene; threshold must be in the
        range used for comparable functional markers (dsrB: 0.35, mcrA: 0.35).

        The molybdopterin-guanine dinucleotide (MGD) binding Cys/Met residues
        of NarG are near-invariant across all respiratory nitrate reducers, but
        synonymous divergence between genera (Paracoccus vs Pseudomonas vs
        Thiobacillus) is high.  Thresholds below 0.25 bits are too strict for
        inter-genus primer placement; above 0.45 bits admits poorly-conserved
        positions that cannot support specific primers.
        """
        assert nrb_cfg.entropy_threshold >= 0.25, (
            "narG diverges across NRB genera; entropy_threshold < 0.25 will "
            "yield no usable conserved windows for a multi-genus input set"
        )
        assert nrb_cfg.entropy_threshold <= 0.45, (
            "narG MGD-binding motifs are near-invariant; entropy_threshold "
            "> 0.45 admits positions that cannot support specific LAMP primers"
        )

    def test_gc_max_covers_gc_rich_nrb_genera(self, nrb_cfg: PipelineConfig) -> None:
        """gc_max must accommodate GC-rich NRB genera.

        Dominant oilfield NRB genera are GC-rich:
          Paracoccus denitrificans ~67% GC
          Pseudomonas stutzeri     ~63% GC
          Thiobacillus denitrificans ~66% GC
        gc_max must be >= 65% to allow primer design in conserved narG blocks.
        """
        assert nrb_cfg.gc_max >= 65.0, (
            "gc_max must be >= 65% to accommodate GC-rich NRB genera "
            "(Paracoccus ~67%, Pseudomonas ~63%, Thiobacillus ~66%)"
        )

    def test_gc_min_provides_primer_design_headroom(self, nrb_cfg: PipelineConfig) -> None:
        """gc_min must be low enough for primer design across GC diversity."""
        assert nrb_cfg.gc_min <= 45.0, (
            "gc_min must be <= 45% to provide headroom for primer design "
            "in conserved narG windows across GC-rich NRB genera"
        )

    def test_gc_window_is_sufficiently_wide(self, nrb_cfg: PipelineConfig) -> None:
        """GC window must span inter-genus GC diversity for narG."""
        assert nrb_cfg.gc_max - nrb_cfg.gc_min >= 20.0, (
            "NRB genera span significant GC diversity; gc window must be "
            "at least 20 GC units wide to support multi-genus primer design"
        )

    def test_max_sequences_spans_genus_diversity(self, nrb_cfg: PipelineConfig) -> None:
        """Need enough sequences to cover Paracoccus species and supplementary genera."""
        assert nrb_cfg.max_sequences >= 20

    def test_standard_lamp_tm_window(self, nrb_cfg: PipelineConfig) -> None:
        """narG is a chromosomal DNA target; standard 60-65 degC LAMP window is correct."""
        assert nrb_cfg.tm_min >= 60.0
        assert nrb_cfg.tm_max <= 65.0

    def test_not_rt_lamp(self, nrb_cfg: PipelineConfig) -> None:
        """narG is a chromosomal DNA target; tm_min must not be raised to RT-LAMP floor.

        The RT-LAMP one-step 63 degC floor (used for PRRSV, FMDV, AIV RNA
        targets) is not appropriate for a DNA target like narG.  Using it
        unnecessarily narrows the primer candidate pool without benefit.
        """
        assert nrb_cfg.tm_min < 63.5, (
            "NRB / narG is a DNA target: tm_min should be at the standard "
            "LAMP 60 degC, not the RT-LAMP co-activity floor of 63 degC"
        )

    def test_amplicon_size_within_lamp_window(self, nrb_cfg: PipelineConfig) -> None:
        assert 80 <= nrb_cfg.f2_b2_min < nrb_cfg.f2_b2_max <= 200

    def test_min_region_length_fits_amplicon(self, nrb_cfg: PipelineConfig) -> None:
        """min_region_length must leave room for F3/B3 flanking primers."""
        assert nrb_cfg.min_region_length >= nrb_cfg.f2_b2_max + 40

    def test_output_dir_references_nrb(self, nrb_cfg: PipelineConfig) -> None:
        assert "nrb" in str(nrb_cfg.output_dir).lower()

    def test_identity_threshold_is_standard(self, nrb_cfg: PipelineConfig) -> None:
        """Off-target identity threshold must be in the recommended range."""
        assert 0.70 <= nrb_cfg.min_identity_threshold <= 0.95

    def test_coverage_threshold_is_standard(self, nrb_cfg: PipelineConfig) -> None:
        """Off-target coverage threshold must be in the recommended range."""
        assert 0.70 <= nrb_cfg.min_coverage_threshold <= 0.95
