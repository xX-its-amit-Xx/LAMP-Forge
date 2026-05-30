"""Shared pytest fixtures.

A note on test data: we keep fixtures small and synthetic where possible so
the test suite runs offline and fast. The Mycobacterium fixture is a tiny
hand-built MSA of synthetic sequences with a deliberately conserved 300 bp
middle region — enough to exercise the full pipeline without hitting NCBI.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from lamp_forge.types import PipelineConfig

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    """Path to bundled test FASTAs."""
    return FIXTURES_DIR


@pytest.fixture
def conserved_msa() -> MultipleSeqAlignment:
    """A small synthetic MSA with a clearly conserved middle region.

    Construction: 5 sequences, 600 bp each. The middle 300 bp (positions
    150-450) is identical across all sequences. Outside that range, each
    position is randomised across {A,C,G,T} but with a 60% bias to the
    consensus base — entropy will be ~0.97 bits there vs. 0.0 in the
    conserved middle.
    """
    import random

    random.seed(42)
    consensus = (
        "ATGCGTAACCGGTTATCGATCGATCGATGCATGCATCGATCGATCGGTACGTACGTAGCATGCATGCAT"  # 0-67
        "CGATCGATCGATCGATCGATCGATGCATGCATCGATCGATCGATGCATCGATGCATGCATCGAT"  # 68-131
        "GCATCGTAGCATCGTAGCATCGCAT"  # 132-155
        + "A" * 50
        + "CGATCGTGACGTACGATCGATCGTAGCATCGATGCATCGTAGCATGCATCGATCGATGCATCGATGCAT"  # 200-269
        + "CGATCGTAGCATCGATCGTAGCATGCATGCATCGATGCATCGTAGCATCGATGCATCGTAGCATGCAT"  # 270-339
        + "GCATCGATCGATGCATGCATCGATGCATCGATCGTAGCATCGATGCATCGATCGTAGCATCGATGCAT"  # 340-409
        + "C" * 40
        + "CGTAGCATCGTAGCATCGATGCATCGTAGCATCGATGCATCGTAGCATGCATCGATGCATCGATGCAT"  # 450-517
        + "CGATCGATGCATCGATCGTAGCATCGATGCATCGATCGTAGCATGCATCGATGCATCGATCGTAGCAT"  # 518-585
        + "AGCATGCATGC"  # 586-596
    )[:600].ljust(600, "A")

    def _mutate_position(base: str, conserved: bool) -> str:
        if conserved:
            return base
        if random.random() < 0.6:
            return base
        return random.choice([b for b in "ACGT" if b != base])

    records: list[SeqRecord] = []
    for idx in range(5):
        seq_chars = []
        for i, b in enumerate(consensus):
            conserved = 150 <= i < 450
            seq_chars.append(_mutate_position(b, conserved))
        records.append(SeqRecord(Seq("".join(seq_chars)), id=f"seq_{idx:02d}"))

    return MultipleSeqAlignment(records)


@pytest.fixture
def conserved_msa_fasta(conserved_msa: MultipleSeqAlignment, tmp_path: Path) -> Path:
    """Same synthetic MSA, persisted as an aligned FASTA file."""
    from Bio import SeqIO

    path = tmp_path / "msa.fasta"
    SeqIO.write(conserved_msa, path, "fasta")
    return path


@pytest.fixture
def example_config(tmp_path: Path) -> PipelineConfig:
    """Minimal PipelineConfig pointing at writeable tmp paths."""
    cache = tmp_path / "cache"
    output = tmp_path / "results"
    off_targets = tmp_path / "off_targets"
    cache.mkdir(exist_ok=True)
    output.mkdir(exist_ok=True)
    off_targets.mkdir(exist_ok=True)
    return PipelineConfig(
        target_name="test_target",
        taxon_id=1773,
        accessions=[],
        gene="rpoB",
        max_sequences=5,
        email="test@example.com",
        off_target_dir=off_targets,
        min_identity_threshold=0.80,
        min_coverage_threshold=0.80,
        window_size=30,
        entropy_threshold=0.25,
        min_region_length=200,
        tm_min=58.0,  # slight relaxation for synthetic fixtures
        tm_max=68.0,
        tm_match_tolerance=3.0,
        gc_min=30.0,
        gc_max=70.0,
        hairpin_dg_threshold=-3.0,
        dimer_dg_threshold=-6.0,
        f2_b2_min=120,
        f2_b2_max=160,
        output_dir=output,
        top_n=5,
        generate_html=False,
        generate_csv=True,
        cache_dir=cache,
    )


@pytest.fixture
def off_target_fastas(tmp_path: Path) -> Path:
    """A directory of synthetic off-target FASTAs for specificity tests."""
    d = tmp_path / "off_targets"
    d.mkdir(exist_ok=True)
    # off-target #1: shares a 25 nt stretch with our consensus to provoke a hit
    (d / "off_target_1.fasta").write_text(
        ">off1\n"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        "ATGCGTAACCGGTTATCGATCGATCGAT"  # ~28nt that may match
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
    )
    # off-target #2: random, should not match
    (d / "off_target_2.fasta").write_text(
        ">off2\n" + "CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA" * 5 + "\n"
    )
    return d


@pytest.fixture
def has_mafft() -> bool:
    """Skip marker: True if the ``mafft`` binary is on PATH."""
    return shutil.which("mafft") is not None


@pytest.fixture
def has_blast() -> bool:
    """Skip marker: True if both ``makeblastdb`` and ``blastn`` are on PATH."""
    return shutil.which("makeblastdb") is not None and shutil.which("blastn") is not None
