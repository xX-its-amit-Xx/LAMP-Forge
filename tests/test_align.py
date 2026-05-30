"""Tests for the MAFFT wrapper.

The interesting tests need ``mafft`` on PATH; the rest exercise file-loading
and validation paths that work offline.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from lamp_forge.align import (
    load_alignment,
    run_mafft,
)


class TestLoadAlignment:
    def test_loads_valid_alignment(self, conserved_msa_fasta: Path) -> None:
        msa = load_alignment(conserved_msa_fasta)
        assert len(msa) >= 2

    def test_single_sequence_rejected(self, tmp_path: Path) -> None:
        path = tmp_path / "single.fasta"
        SeqIO.write([SeqRecord(Seq("ACGT"), id="s1")], path, "fasta")
        with pytest.raises(ValueError, match="only 1"):
            load_alignment(path)

    def test_uneven_lengths_rejected(self, tmp_path: Path) -> None:
        path = tmp_path / "uneven.fasta"
        SeqIO.write(
            [
                SeqRecord(Seq("ACGT"), id="s1"),
                SeqRecord(Seq("ACGTACGT"), id="s2"),
            ],
            path,
            "fasta",
        )
        with pytest.raises(ValueError, match="unequal length"):
            load_alignment(path)


class TestRunMafft:
    def test_missing_input_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            run_mafft(tmp_path / "missing.fa", tmp_path / "out.fa")

    def test_runs_if_mafft_available(
        self, has_mafft: bool, tmp_path: Path, fixtures_dir: Path
    ) -> None:
        if not has_mafft:
            pytest.skip("MAFFT not on PATH")
        out = tmp_path / "aligned.fasta"
        run_mafft(fixtures_dir / "mycobacterium_small.fasta", out)
        assert out.exists()
        records = list(SeqIO.parse(out, "fasta"))
        assert len(records) >= 3
        # All sequences should have the same length after alignment.
        lengths = {len(r.seq) for r in records}
        assert len(lengths) == 1
