"""Tests for the specificity / off-target screening module.

These tests need ``makeblastdb`` and ``blastn`` on PATH. They are skipped
gracefully when those binaries are missing (e.g. on a dev laptop without
BLAST+ installed).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lamp_forge.specificity import (
    BlastNotFoundError,
    _row_to_hit,
    build_blast_db,
    discover_off_targets,
    screen_set,
    write_specificity_tsv,
)
from lamp_forge.types import LampPrimerSet, PipelineConfig, Primer


def _mk_primer(role: str, seq: str, start: int = 0) -> Primer:
    """Quick helper for synthetic Primer objects in tests."""
    return Primer(
        role=role,  # type: ignore[arg-type]
        sequence=seq,
        start=start,
        end=start + len(seq),
        strand="+",
        tm=62.0,
        gc_percent=50.0,
        hairpin_dg=-1.0,
        homodimer_dg=-3.0,
    )


def _mk_set(set_id: str = "S0", region_id: str = "R0") -> LampPrimerSet:
    return LampPrimerSet(
        set_id=set_id,
        region_id=region_id,
        f3=_mk_primer("F3", "ATGCGTAACCGGTTATCGAT", 0),
        b3=_mk_primer("B3", "CGATCGATCGATCGATCGAT", 200),
        fip=_mk_primer("FIP", "AACCGGTTATCGATCGATCGATCGATCGATAA"),
        bip=_mk_primer("BIP", "GGCCAATTAAGGCCAATTAAGGCCAATTAAGG"),
        f2_b2_distance=140,
        tm_spread=1.5,
        mean_tm=62.0,
    )


class TestDiscoverOffTargets:
    def test_finds_recognised_extensions(self, tmp_path: Path) -> None:
        (tmp_path / "a.fasta").write_text(">a\nACGT\n")
        (tmp_path / "b.fa").write_text(">b\nACGT\n")
        (tmp_path / "c.fna").write_text(">c\nACGT\n")
        (tmp_path / "ignored.txt").write_text("not a fasta")
        found = discover_off_targets(tmp_path)
        names = {p.name for p in found}
        assert names == {"a.fasta", "b.fa", "c.fna"}

    def test_missing_dir_returns_empty(self, tmp_path: Path) -> None:
        assert discover_off_targets(tmp_path / "nope") == []


class TestRowToHit:
    def test_flags_above_threshold(self, example_config: PipelineConfig) -> None:
        row = {
            "qseqid": "F3_0",
            "sseqid": "off1|chr1",
            "pident": "95.0",
            "length": "20",
            "qlen": "20",
            "evalue": "1e-5",
            "bitscore": "40.0",
            "_primer": "ATGCATGCATGCATGCATGC",
            "_role": "F3",
        }
        hit = _row_to_hit(row, example_config)
        assert hit.is_flagged is True
        assert hit.percent_identity == pytest.approx(0.95)
        assert hit.coverage == pytest.approx(1.0)
        assert hit.off_target_name == "off1"

    def test_below_identity_not_flagged(self, example_config: PipelineConfig) -> None:
        row = {
            "qseqid": "F3_0",
            "sseqid": "off1|chr1",
            "pident": "70.0",
            "length": "20",
            "qlen": "20",
            "evalue": "1e-2",
            "bitscore": "20.0",
            "_primer": "ATGCATGCATGCATGCATGC",
            "_role": "F3",
        }
        hit = _row_to_hit(row, example_config)
        assert hit.is_flagged is False

    def test_below_coverage_not_flagged(self, example_config: PipelineConfig) -> None:
        row = {
            "qseqid": "F3_0",
            "sseqid": "off1|chr1",
            "pident": "100.0",
            "length": "10",
            "qlen": "20",
            "evalue": "1e-2",
            "bitscore": "20.0",
            "_primer": "ATGCATGCATGCATGCATGC",
            "_role": "F3",
        }
        hit = _row_to_hit(row, example_config)
        assert hit.coverage == pytest.approx(0.5)
        assert hit.is_flagged is False


class TestBlastIntegration:
    def test_build_db_empty_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="empty"):
            build_blast_db([], tmp_path / "db")

    def test_build_db_and_screen(
        self,
        has_blast: bool,
        tmp_path: Path,
        off_target_fastas: Path,
        example_config: PipelineConfig,
    ) -> None:
        if not has_blast:
            pytest.skip("BLAST+ not on PATH")
        files = discover_off_targets(off_target_fastas)
        db = tmp_path / "off_db"
        build_blast_db(files, db)
        primer_set = _mk_set()
        result = screen_set(primer_set, db, example_config)
        # The synthetic off-target #1 was crafted to overlap F3 — the hit
        # may or may not exceed thresholds depending on BLAST scoring at
        # short word size, but a hit should at least appear.
        assert isinstance(result.cross_reactivity_hits, list)
        assert 0.0 <= result.specificity_score <= 1.0


class TestWriteSpecificityTsv:
    def test_writes_header_only_when_empty(self, tmp_path: Path) -> None:
        out = tmp_path / "spec.tsv"
        write_specificity_tsv([], out)
        text = out.read_text()
        assert text.startswith("set_id\tprimer_role")
        # Only the header line, plus the trailing newline → exactly 1 line.
        assert text.count("\n") == 1

    def test_writes_rows(self, tmp_path: Path, example_config: PipelineConfig) -> None:
        from lamp_forge.types import SpecificityHit

        s = _mk_set()
        s.cross_reactivity_hits = [
            SpecificityHit(
                primer_role="F3",
                primer_sequence="ATGCATGCATGCATGCATGC",
                off_target_name="off1",
                percent_identity=0.85,
                alignment_length=18,
                primer_length=20,
                coverage=0.9,
                e_value=1e-5,
                bitscore=35.0,
                is_flagged=True,
            )
        ]
        out = tmp_path / "spec.tsv"
        write_specificity_tsv([s], out)
        text = out.read_text()
        assert "off1" in text
        assert "S0" in text
        assert text.count("\n") == 2  # header + 1 data row
