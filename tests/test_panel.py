"""Tests for multiplex panel compatibility analysis.

The selection logic is exercised with an injected ΔG function so it is fast and
deterministic (no primer3 dependency). One end-to-end CLI test uses the real
primer3-backed engine to confirm the wiring.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from lamp_forge.panel import (
    PanelCandidateSet,
    PanelOligo,
    cross_dimer_matrix,
    load_candidate_sets,
    render_panel_html,
    run_panel,
    select_panel,
    write_panel_csv,
    write_panel_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _primer_dict(role: str, seq: str) -> dict:
    """Minimal serialised Primer matching report.write_json output."""
    return {
        "role": role,
        "sequence": seq,
        "start": 0,
        "end": len(seq),
        "strand": "+",
        "tm": 62.0,
        "gc_percent": 50.0,
        "hairpin_dg": -1.0,
        "homodimer_dg": -2.0,
        "inner_subseq": None,
        "outer_subseq": None,
    }


def _write_sets_json(path: Path, sets: list[tuple[str, float, dict[str, str]]]) -> None:
    """Write a primer_sets.json. ``sets`` = list of (set_id, score, {field: seq})."""
    fields = ("f3", "b3", "fip", "bip", "lf", "lb")
    entries = []
    for set_id, score, field_seqs in sets:
        entry: dict = {"set_id": set_id, "region_id": "R000", "composite_score": score}
        for f in fields:
            seq = field_seqs.get(f)
            entry[f] = _primer_dict(f.upper(), seq) if seq else None
        entries.append(entry)
    payload = {"version": "0.1.0", "generated_at": "2026-01-01T00:00:00Z", "primer_sets": entries}
    path.write_text(json.dumps(payload), encoding="utf-8")


def _make_dg(dg_map: dict[frozenset[str], float], default: float = -2.0):
    """Build a deterministic ΔG function from an explicit sequence-pair map."""

    def dg(a: str, b: str) -> float:
        return dg_map.get(frozenset((a, b)), default)

    return dg


def _candidate(target: str, set_id: str, score: float, seqs: list[str]) -> PanelCandidateSet:
    oligos = tuple(
        PanelOligo(target=target, set_id=set_id, role=f"R{i}", sequence=s)
        for i, s in enumerate(seqs)
    )
    return PanelCandidateSet(target=target, set_id=set_id, composite_score=score, oligos=oligos)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


class TestLoadCandidateSets:
    def test_sorts_by_composite_and_truncates(self, tmp_path: Path) -> None:
        p = tmp_path / "sets.json"
        _write_sets_json(
            p,
            [
                ("S0", 0.50, {"f3": "AAAA", "b3": "CCCC"}),
                ("S1", 0.90, {"f3": "GGGG", "b3": "TTTT"}),
                ("S2", 0.70, {"f3": "ACGT", "b3": "TGCA"}),
            ],
        )
        cands = load_candidate_sets(p, "TGT", top_n=2)
        assert [c.set_id for c in cands] == ["S1", "S2"]  # best-first, truncated
        assert cands[0].oligos[0].target == "TGT"

    def test_four_primer_set_has_four_oligos(self, tmp_path: Path) -> None:
        p = tmp_path / "sets.json"
        _write_sets_json(
            p, [("S0", 0.8, {"f3": "AAAA", "b3": "CCCC", "fip": "GGGG", "bip": "TTTT"})]
        )
        cands = load_candidate_sets(p, "TGT")
        assert len(cands[0].oligos) == 4

    def test_six_primer_set_includes_loops(self, tmp_path: Path) -> None:
        p = tmp_path / "sets.json"
        _write_sets_json(
            p,
            [
                (
                    "S0",
                    0.8,
                    {
                        "f3": "AAAA",
                        "b3": "CCCC",
                        "fip": "GGGG",
                        "bip": "TTTT",
                        "lf": "ACAC",
                        "lb": "GTGT",
                    },
                )
            ],
        )
        cands = load_candidate_sets(p, "TGT")
        roles = {o.role for o in cands[0].oligos}
        assert roles == {"F3", "B3", "FIP", "BIP", "LF", "LB"}

    def test_raises_when_no_usable_sets(self, tmp_path: Path) -> None:
        p = tmp_path / "empty.json"
        p.write_text(json.dumps({"primer_sets": []}), encoding="utf-8")
        with pytest.raises(ValueError, match="no usable primer sets"):
            load_candidate_sets(p, "TGT")


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


class TestSelectPanel:
    def test_requires_two_targets(self) -> None:
        cands = {"A": [_candidate("A", "A0", 0.9, ["AAAA"])]}
        with pytest.raises(ValueError, match="at least 2 targets"):
            select_panel(cands, dg_func=_make_dg({}))

    def test_avoids_the_bad_pairing(self) -> None:
        # The A1/B1 pairing is a terrible dimer; everything else is mild (-2 default).
        a1 = _candidate("A", "A1", 0.9, ["a1x"])
        a2 = _candidate("A", "A2", 0.5, ["a2x"])
        b1 = _candidate("B", "B1", 0.9, ["b1x"])
        b2 = _candidate("B", "B2", 0.5, ["b2x"])
        dg = _make_dg({frozenset(("a1x", "b1x")): -12.0})
        result = select_panel({"A": [a1, a2], "B": [b1, b2]}, dimer_dg_threshold=-5.0, dg_func=dg)
        assert result.search_mode == "exhaustive"
        assert result.worst_dg == -2.0  # the bad -12 pairing was avoided
        assert result.is_clean
        # Whatever was chosen, it must not be the A1+B1 combination.
        assert not (result.selection["A"].set_id == "A1" and result.selection["B"].set_id == "B1")

    def test_flags_when_every_combo_is_bad(self) -> None:
        # Both A candidates clash with the only B candidate.
        a1 = _candidate("A", "A1", 0.9, ["a1x"])
        a2 = _candidate("A", "A2", 0.5, ["a2x"])
        b1 = _candidate("B", "B1", 0.9, ["b1x"])
        dg = _make_dg({frozenset(("a1x", "b1x")): -9.0, frozenset(("a2x", "b1x")): -8.0})
        result = select_panel({"A": [a1, a2], "B": [b1]}, dimer_dg_threshold=-5.0, dg_func=dg)
        assert not result.is_clean
        # Picks the less-bad of the two (A2: -8 beats A1: -9).
        assert result.worst_dg == -8.0
        assert result.selection["A"].set_id == "A2"
        assert len(result.flagged) == 1

    def test_greedy_path_reaches_optimum(self) -> None:
        # Force greedy by capping combinations below the product (2*2=4).
        a1 = _candidate("A", "A1", 0.9, ["a1x"])  # top composite → greedy start
        a2 = _candidate("A", "A2", 0.5, ["a2x"])
        b1 = _candidate("B", "B1", 0.9, ["b1x"])  # top composite → greedy start
        b2 = _candidate("B", "B2", 0.5, ["b2x"])
        dg = _make_dg({frozenset(("a1x", "b1x")): -12.0})  # the greedy start pair is the bad one
        result = select_panel(
            {"A": [a1, a2], "B": [b1, b2]},
            dimer_dg_threshold=-5.0,
            dg_func=dg,
            max_combinations=1,
        )
        assert result.search_mode == "greedy"
        assert result.worst_dg == -2.0  # single swap escapes the bad start
        assert result.is_clean


# ---------------------------------------------------------------------------
# Matrix + outputs
# ---------------------------------------------------------------------------


class TestMatrixAndOutputs:
    def test_cross_dimer_matrix_masks_intra_and_diagonal(self) -> None:
        import math

        a = _candidate("A", "A1", 0.9, ["a1x", "a1y"])  # two oligos, same target
        b = _candidate("B", "B1", 0.9, ["b1x"])
        result = select_panel({"A": [a], "B": [b]}, dg_func=_make_dg({}), dimer_dg_threshold=-5.0)
        from lamp_forge.panel import selected_oligos

        oligos = selected_oligos(result)
        m = cross_dimer_matrix(oligos, dg_func=_make_dg({}))
        n = len(oligos)
        assert n == 3
        for i in range(n):
            assert math.isnan(m[i][i])  # diagonal masked
        # the two A oligos are intra-target → masked
        assert math.isnan(m[0][1])
        # A↔B entries are real
        assert m[0][2] == -2.0

    def test_write_json_and_csv(self, tmp_path: Path) -> None:
        a1 = _candidate("A", "A1", 0.9, ["a1x"])
        b1 = _candidate("B", "B1", 0.8, ["b1x"])
        dg = _make_dg({frozenset(("a1x", "b1x")): -9.0})
        result = select_panel({"A": [a1], "B": [b1]}, dimer_dg_threshold=-5.0, dg_func=dg)

        jpath = tmp_path / "panel.json"
        write_panel_json(result, jpath)
        data = json.loads(jpath.read_text())
        assert data["is_clean"] is False
        assert data["worst_inter_assay_dg"] == -9.0
        assert len(data["selection"]) == 2
        assert len(data["flagged_cross_dimers"]) == 1
        assert data["flagged_cross_dimers"][0]["dg_kcal_mol"] == -9.0

        cpath = tmp_path / "panel.csv"
        write_panel_csv(result, cpath)
        lines = cpath.read_text().strip().splitlines()
        assert lines[0].startswith("target,set_id,role")
        assert len(lines) == 1 + 2  # header + 2 oligos

    def test_render_html_writes_file(self, tmp_path: Path) -> None:
        a1 = _candidate("A", "A1", 0.9, ["a1x"])
        b1 = _candidate("B", "B1", 0.8, ["b1x"])
        result = select_panel({"A": [a1], "B": [b1]}, dg_func=_make_dg({}))
        out = tmp_path / "panel.html"
        render_panel_html(result, out, dg_func=_make_dg({}))
        html = out.read_text(encoding="utf-8")
        assert "multiplex panel report" in html
        assert "A1" in html and "B1" in html


class TestRunPanel:
    def test_end_to_end_writes_all_artifacts(self, tmp_path: Path) -> None:
        a = tmp_path / "a.json"
        b = tmp_path / "b.json"
        _write_sets_json(a, [("A0", 0.9, {"f3": "a_f3", "b3": "a_b3"})])
        _write_sets_json(b, [("B0", 0.9, {"f3": "b_f3", "b3": "b_b3"})])
        out = tmp_path / "panel_out"
        result = run_panel(
            {"A": a, "B": b},
            out,
            top_per_target=3,
            dimer_dg_threshold=-5.0,
            generate_html=True,
            dg_func=_make_dg({}),
        )
        assert result.is_clean
        assert (out / "panel.json").exists()
        assert (out / "panel_primers.csv").exists()
        assert (out / "panel_report.html").exists()


# ---------------------------------------------------------------------------
# CLI (real primer3 engine)
# ---------------------------------------------------------------------------


class TestPanelCLI:
    def _write_real_sets(self, path: Path, prefix: str) -> None:
        # Realistic ~20mers so the real primer3 heterodimer engine has something
        # to chew on; exact ΔG values don't matter to the assertions.
        _write_sets_json(
            path,
            [
                (
                    f"{prefix}_S0",
                    0.9,
                    {
                        "f3": "GACCTGACATTCGGATCACG",
                        "b3": "TTGCAGTCGAACTGGTACGA",
                        "fip": "CAGTTCGACTGCAAGGTACCGATCGAATGTCAGGTC",
                        "bip": "GTACCAGTTCGACTGCAATGCGTGATCCGAATGTCA",
                    },
                ),
            ],
        )

    def test_panel_command_happy_path(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        a = tmp_path / "a.json"
        b = tmp_path / "b.json"
        self._write_real_sets(a, "A")
        self._write_real_sets(b, "B")
        out = tmp_path / "panel"

        runner = CliRunner()
        res = runner.invoke(
            cli,
            ["panel", "--set", f"A={a}", "--set", f"B={b}", "--out", str(out), "--no-html"],
            obj={},
        )
        assert res.exit_code == 0, res.output
        assert (out / "panel.json").exists()
        assert (out / "panel_primers.csv").exists()

    def test_panel_command_rejects_bad_set_spec(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        a = tmp_path / "a.json"
        self._write_real_sets(a, "A")
        runner = CliRunner()
        # Missing '=' in the second spec.
        res = runner.invoke(cli, ["panel", "--set", f"A={a}", "--set", "B"], obj={})
        assert res.exit_code == 2
        assert "LABEL=PATH" in res.output

    def test_panel_command_requires_two_targets(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        a = tmp_path / "a.json"
        self._write_real_sets(a, "A")
        runner = CliRunner()
        res = runner.invoke(cli, ["panel", "--set", f"A={a}"], obj={})
        assert res.exit_code == 2
        assert "at least 2" in res.output
