"""Tests for serialisation and HTML/JSON/CSV report writers."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from lamp_forge.conserve import compute_track, find_conserved_regions
from lamp_forge.report import (
    plot_conservation,
    plot_specificity_heatmap,
    render_html,
    write_conservation_tsv,
    write_csv,
    write_json,
    write_manifest,
)
from lamp_forge.types import LampPrimerSet, PipelineConfig, Primer


def _mk_primer(role: str, seq: str = "ATGCATGCATGCATGCATGC") -> Primer:
    return Primer(
        role=role,  # type: ignore[arg-type]
        sequence=seq,
        start=0,
        end=len(seq),
        strand="+",
        tm=62.0,
        gc_percent=50.0,
        hairpin_dg=-1.0,
        homodimer_dg=-3.0,
    )


def _mk_set(set_id: str = "S0", region_id: str = "R0") -> LampPrimerSet:
    s = LampPrimerSet(
        set_id=set_id,
        region_id=region_id,
        f3=_mk_primer("F3"),
        b3=_mk_primer("B3"),
        fip=_mk_primer("FIP", "AACCGGTTATCGATCGATCGATCGATCGATAA"),
        bip=_mk_primer("BIP", "GGCCAATTAAGGCCAATTAAGGCCAATTAAGG"),
        f2_b2_distance=140,
        tm_spread=1.5,
        mean_tm=62.0,
        conservation_score=0.85,
        specificity_score=0.9,
        thermodynamic_score=0.8,
        composite_score=0.86,
    )
    return s


class TestWriteJson:
    def test_writes_parseable_json(self, tmp_path: Path) -> None:
        path = tmp_path / "p.json"
        write_json([_mk_set(), _mk_set("S1")], path)
        payload = json.loads(path.read_text())
        assert payload["version"] == "0.1.0"
        assert len(payload["primer_sets"]) == 2
        assert payload["primer_sets"][0]["set_id"] == "S0"

    def test_empty_list(self, tmp_path: Path) -> None:
        path = tmp_path / "p.json"
        write_json([], path)
        payload = json.loads(path.read_text())
        assert payload["primer_sets"] == []


class TestWriteCsv:
    def test_one_row_per_primer(self, tmp_path: Path) -> None:
        path = tmp_path / "p.csv"
        write_csv([_mk_set()], path)
        with path.open() as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        # 4 primers (F3, B3, FIP, BIP); LF/LB are None on the test set.
        assert len(rows) == 4
        roles = {r["primer_role"] for r in rows}
        assert roles == {"F3", "B3", "FIP", "BIP"}
        for r in rows:
            assert r["set_id"] == "S0"
            assert r["primer_name"] == f"S0_{r['primer_role']}"


class TestPlotConservation:
    def test_returns_data_uri(self, conserved_msa_fasta: Path) -> None:
        from lamp_forge.align import load_alignment

        msa = load_alignment(conserved_msa_fasta)
        track = compute_track(msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.20, min_region_length=100)
        uri = plot_conservation(track, regions, 0.20)
        assert uri.startswith("data:image/png;base64,")
        assert len(uri) > 200  # crude sanity that there's actually a PNG


class TestPlotSpecificityHeatmap:
    def test_empty_sets_returns_placeholder(self) -> None:
        uri = plot_specificity_heatmap([], top_n=10)
        assert uri.startswith("data:image/png;base64,")


class TestWriteConservationTsv:
    def test_one_row_per_position(self, conserved_msa_fasta: Path, tmp_path: Path) -> None:
        from lamp_forge.align import load_alignment

        msa = load_alignment(conserved_msa_fasta)
        track = compute_track(msa, window_size=30)
        out = tmp_path / "c.tsv"
        write_conservation_tsv(track, out)
        lines = out.read_text().strip().splitlines()
        # header + width rows
        assert len(lines) == len(track.raw_entropy) + 1


class TestRenderHtml:
    def test_renders_with_minimal_data(
        self,
        conserved_msa_fasta: Path,
        example_config: PipelineConfig,
        tmp_path: Path,
    ) -> None:
        from lamp_forge.align import load_alignment

        msa = load_alignment(conserved_msa_fasta)
        track = compute_track(msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.20, min_region_length=100)
        # Need at least one region for the report to render a detail block,
        # but the HTML should render even with zero sets.
        sets = []
        if regions:
            s = _mk_set(region_id=regions[0].region_id)
            # Place primers inside the region for the ladder plot.
            r = regions[0]
            s.f3 = Primer("F3", s.f3.sequence, r.start, r.start + 20, "+", 62, 50, -1, -3)
            s.b3 = Primer("B3", s.b3.sequence, r.end - 20, r.end, "-", 62, 50, -1, -3)
            s.fip = Primer(
                "FIP",
                s.fip.sequence,
                r.start + 30,
                r.start + 30 + len(s.fip.sequence),
                "+",
                62,
                50,
                -1,
                -3,
            )
            s.bip = Primer(
                "BIP",
                s.bip.sequence,
                r.end - 30 - len(s.bip.sequence),
                r.end - 30,
                "-",
                62,
                50,
                -1,
                -3,
            )
            sets = [s]
        out = tmp_path / "report.html"
        render_html(sets, regions, track, example_config, 5, 0, out)
        assert out.exists()
        html = out.read_text()
        assert "<title>LAMP-Forge report" in html
        assert "test_target" in html


class TestWriteManifest:
    def test_includes_version_and_config(
        self, example_config: PipelineConfig, tmp_path: Path
    ) -> None:
        path = tmp_path / "m.json"
        write_manifest(example_config, {"count": 5, "records": []}, 2, 7, path)
        payload = json.loads(path.read_text())
        assert payload["lamp_forge_version"] == "0.1.0"
        assert payload["n_conserved_regions"] == 2
        assert payload["n_primer_sets"] == 7
        assert payload["config"]["target_name"] == "test_target"
