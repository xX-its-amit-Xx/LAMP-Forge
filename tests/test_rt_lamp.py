"""Tests for the RT-LAMP compatibility checker.

All tests are pure Python — no external binaries required.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from lamp_forge.rt_lamp import (
    RtLampCheckResult,
    RtLampParams,
    TargetNucleicAcid,
    check_primer_set_for_rt_lamp,
    check_primer_sets_for_rt_lamp,
    write_rt_check_csv,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_set_dict(
    set_id: str = "region_01_set_001",
    f3_tm: float = 64.0,
    b3_tm: float = 64.0,
    fip_tm: float = 64.0,
    bip_tm: float = 64.0,
    lf_tm: float | None = 62.0,
    lb_tm: float | None = 62.0,
) -> dict[str, Any]:
    """Minimal primer set dict mirroring the primer_sets.json schema."""
    d: dict[str, Any] = {
        "set_id": set_id,
        "composite_score": 0.9,
        "f3": {"role": "F3", "sequence": "ATCGATCGATCG", "tm": f3_tm},
        "b3": {"role": "B3", "sequence": "CGATCGATCGAT", "tm": b3_tm},
        "fip": {
            "role": "FIP",
            "sequence": "ATCGATCGATCGATCGATCGATCGATCGATCG",
            "tm": fip_tm,
        },
        "bip": {
            "role": "BIP",
            "sequence": "CGATCGATCGATCGATCGATCGATCGATCGAT",
            "tm": bip_tm,
        },
    }
    if lf_tm is not None:
        d["lf"] = {"role": "LF", "sequence": "ATCGATCGATCG", "tm": lf_tm}
    if lb_tm is not None:
        d["lb"] = {"role": "LB", "sequence": "CGATCGATCGAT", "tm": lb_tm}
    return d


# ---------------------------------------------------------------------------
# RtLampParams
# ---------------------------------------------------------------------------


class TestRtLampParams:
    def test_defaults(self) -> None:
        p = RtLampParams()
        assert p.rt_min_tm == 63.0
        assert p.rt_max_tm == 65.0
        assert p.loop_min_tm == 60.0

    def test_custom_values(self) -> None:
        p = RtLampParams(rt_min_tm=62.0, rt_max_tm=66.0, loop_min_tm=59.0)
        assert p.rt_min_tm == 62.0

    def test_raises_when_min_ge_max(self) -> None:
        with pytest.raises(ValueError, match="rt_min_tm"):
            RtLampParams(rt_min_tm=65.0, rt_max_tm=65.0)

    def test_raises_when_loop_min_above_rt_min(self) -> None:
        with pytest.raises(ValueError, match="loop_min_tm"):
            RtLampParams(rt_min_tm=63.0, loop_min_tm=64.0)

    def test_loop_min_equal_to_rt_min_is_valid(self) -> None:
        p = RtLampParams(rt_min_tm=63.0, loop_min_tm=63.0)
        assert p.loop_min_tm == 63.0


# ---------------------------------------------------------------------------
# TargetNucleicAcid
# ---------------------------------------------------------------------------


class TestTargetNucleicAcid:
    def test_rna_value(self) -> None:
        assert TargetNucleicAcid.RNA == "rna"

    def test_dna_value(self) -> None:
        assert TargetNucleicAcid.DNA == "dna"

    def test_from_string_rna(self) -> None:
        assert TargetNucleicAcid("rna") is TargetNucleicAcid.RNA

    def test_from_string_dna(self) -> None:
        assert TargetNucleicAcid("dna") is TargetNucleicAcid.DNA


# ---------------------------------------------------------------------------
# check_primer_set_for_rt_lamp
# ---------------------------------------------------------------------------


class TestCheckPrimerSetForRtLamp:
    def test_all_in_range_is_optimized(self) -> None:
        d = _make_set_dict(f3_tm=64.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d)
        assert result.is_rt_optimized is True
        assert result.n_core_below_rt_min == 0
        assert result.warnings == ()
        assert result.recommended_action == ""

    def test_one_core_below_rt_min(self) -> None:
        d = _make_set_dict(f3_tm=61.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d)
        assert result.is_rt_optimized is False
        assert result.n_core_below_rt_min == 1
        assert len(result.warnings) == 1
        assert "F3" in result.warnings[0]
        assert "61.0" in result.warnings[0]

    def test_two_core_primers_below_rt_min(self) -> None:
        d = _make_set_dict(f3_tm=61.0, b3_tm=61.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_core_below_rt_min == 2
        assert result.is_rt_optimized is False
        assert "2 core primers" in result.recommended_action

    def test_loop_primer_below_loop_min_does_not_affect_core_count(self) -> None:
        d = _make_set_dict(
            f3_tm=64.0,
            b3_tm=64.0,
            fip_tm=64.0,
            bip_tm=64.0,
            lf_tm=59.0,
            lb_tm=62.0,
        )
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_core_below_rt_min == 0
        assert result.is_rt_optimized is True
        assert len(result.warnings) == 1
        assert "LF" in result.warnings[0]

    def test_primer_above_max_generates_warning(self) -> None:
        d = _make_set_dict(f3_tm=66.0)
        result = check_primer_set_for_rt_lamp(d)
        assert any("F3" in w and "66.0" in w for w in result.warnings)

    def test_primer_above_max_does_not_set_core_below(self) -> None:
        d = _make_set_dict(f3_tm=66.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_core_below_rt_min == 0

    def test_four_primer_set_without_loops(self) -> None:
        d = _make_set_dict(lf_tm=None, lb_tm=None)
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_primers == 4

    def test_six_primer_set_with_loops(self) -> None:
        d = _make_set_dict()
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_primers == 6

    def test_n_in_rt_range_counts_correctly(self) -> None:
        # All 4 core in range, 2 loops in range (62.0 >= 60.0) => 6 in range
        d = _make_set_dict(f3_tm=64.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0, lf_tm=62.0, lb_tm=62.0)
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_in_rt_range == 6

    def test_n_in_rt_range_excludes_out_of_range(self) -> None:
        # f3 is below core min, so 5 in range
        d = _make_set_dict(f3_tm=61.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0, lf_tm=62.0, lb_tm=62.0)
        result = check_primer_set_for_rt_lamp(d)
        assert result.n_in_rt_range == 5

    def test_dna_target_skips_tm_check(self) -> None:
        d = _make_set_dict(f3_tm=61.0, b3_tm=61.0, fip_tm=61.0, bip_tm=61.0)
        result = check_primer_set_for_rt_lamp(d, target_na=TargetNucleicAcid.DNA)
        assert result.is_rt_optimized is True
        assert result.n_core_below_rt_min == 0
        assert result.warnings == ()
        assert "DNA" in result.recommended_action

    def test_dna_target_counts_all_in_range(self) -> None:
        d = _make_set_dict()
        result = check_primer_set_for_rt_lamp(d, target_na=TargetNucleicAcid.DNA)
        assert result.n_in_rt_range == result.n_primers

    def test_set_id_from_dict(self) -> None:
        d = _make_set_dict(set_id="region_02_set_005")
        result = check_primer_set_for_rt_lamp(d)
        assert result.set_id == "region_02_set_005"

    def test_missing_set_id_uses_unknown(self) -> None:
        d: dict[str, Any] = {
            "f3": {"role": "F3", "sequence": "ATCG", "tm": 64.0},
        }
        result = check_primer_set_for_rt_lamp(d)
        assert result.set_id == "unknown"

    def test_target_na_stored(self) -> None:
        d = _make_set_dict()
        result = check_primer_set_for_rt_lamp(d, target_na=TargetNucleicAcid.RNA)
        assert result.target_na is TargetNucleicAcid.RNA

    def test_custom_params(self) -> None:
        params = RtLampParams(rt_min_tm=62.0, rt_max_tm=66.0, loop_min_tm=60.0)
        d = _make_set_dict(f3_tm=62.5, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d, params=params)
        assert result.is_rt_optimized is True

    def test_recommended_action_empty_when_optimized(self) -> None:
        d = _make_set_dict()
        result = check_primer_set_for_rt_lamp(d)
        assert result.is_rt_optimized is True
        assert result.recommended_action == ""

    def test_recommended_action_mentions_min_tm(self) -> None:
        params = RtLampParams(rt_min_tm=64.0, rt_max_tm=66.0)
        d = _make_set_dict(f3_tm=62.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d, params=params)
        assert "64.0" in result.recommended_action

    def test_returns_frozen_result(self) -> None:
        d = _make_set_dict()
        result = check_primer_set_for_rt_lamp(d)
        assert isinstance(result, RtLampCheckResult)
        with pytest.raises(AttributeError):
            result.set_id = "modified"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# check_primer_sets_for_rt_lamp
# ---------------------------------------------------------------------------


class TestCheckPrimerSetsForRtLamp:
    def _make_sets(self) -> list[dict[str, Any]]:
        return [
            _make_set_dict("set_001", f3_tm=64.0),
            _make_set_dict("set_002", f3_tm=61.0),
            _make_set_dict("set_003", f3_tm=64.0),
        ]

    def test_returns_one_result_per_set(self) -> None:
        results = check_primer_sets_for_rt_lamp(self._make_sets())
        assert len(results) == 3

    def test_top_n_limits_results(self) -> None:
        results = check_primer_sets_for_rt_lamp(self._make_sets(), top_n=2)
        assert len(results) == 2

    def test_top_n_zero_returns_empty(self) -> None:
        results = check_primer_sets_for_rt_lamp(self._make_sets(), top_n=0)
        assert results == []

    def test_top_n_none_returns_all(self) -> None:
        results = check_primer_sets_for_rt_lamp(self._make_sets(), top_n=None)
        assert len(results) == 3

    def test_target_na_propagated(self) -> None:
        results = check_primer_sets_for_rt_lamp(self._make_sets(), target_na=TargetNucleicAcid.DNA)
        assert all(r.target_na is TargetNucleicAcid.DNA for r in results)

    def test_correct_sets_are_optimized(self) -> None:
        results = check_primer_sets_for_rt_lamp(self._make_sets())
        assert results[0].is_rt_optimized is True
        assert results[1].is_rt_optimized is False
        assert results[2].is_rt_optimized is True

    def test_empty_input_returns_empty(self) -> None:
        results = check_primer_sets_for_rt_lamp([])
        assert results == []


# ---------------------------------------------------------------------------
# write_rt_check_csv
# ---------------------------------------------------------------------------


class TestWriteRtCheckCsv:
    def _make_results(self) -> list[RtLampCheckResult]:
        sets = [
            _make_set_dict("set_001", f3_tm=64.0),
            _make_set_dict("set_002", f3_tm=61.0),
        ]
        return check_primer_sets_for_rt_lamp(sets)

    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv(self._make_results(), out)
        assert out.exists()

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "rt_check.csv"
        write_rt_check_csv(self._make_results(), out)
        assert out.exists()

    def test_header_columns(self, tmp_path: Path) -> None:
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv(self._make_results(), out)
        rows = list(csv.reader(out.open()))
        assert rows[0] == [
            "set_id",
            "target_na",
            "n_primers",
            "n_in_rt_range",
            "n_core_below_rt_min",
            "is_rt_optimized",
            "warnings",
            "recommended_action",
        ]

    def test_row_count(self, tmp_path: Path) -> None:
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv(self._make_results(), out)
        rows = list(csv.reader(out.open()))
        assert len(rows) == 3  # header + 2 data rows

    def test_set_id_in_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv(self._make_results(), out)
        reader = csv.DictReader(out.open())
        ids = [row["set_id"] for row in reader]
        assert ids == ["set_001", "set_002"]

    def test_is_rt_optimized_column(self, tmp_path: Path) -> None:
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv(self._make_results(), out)
        reader = csv.DictReader(out.open())
        rows = list(reader)
        assert rows[0]["is_rt_optimized"] == "True"
        assert rows[1]["is_rt_optimized"] == "False"

    def test_warnings_joined_by_semicolon(self, tmp_path: Path) -> None:
        d = _make_set_dict("s", f3_tm=61.0, b3_tm=61.0)
        result = check_primer_set_for_rt_lamp(d)
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv([result], out)
        reader = csv.DictReader(out.open())
        row = next(reader)
        assert "; " in row["warnings"] or len(result.warnings) <= 1

    def test_empty_warnings_when_clean(self, tmp_path: Path) -> None:
        d = _make_set_dict("s", f3_tm=64.0, b3_tm=64.0, fip_tm=64.0, bip_tm=64.0)
        result = check_primer_set_for_rt_lamp(d)
        out = tmp_path / "rt_check.csv"
        write_rt_check_csv([result], out)
        reader = csv.DictReader(out.open())
        row = next(reader)
        assert row["warnings"] == ""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestRtCheckCli:
    def _make_primer_sets_json(self, tmp_path: Path) -> Path:
        sets = [
            _make_set_dict("set_001", f3_tm=64.0),
            _make_set_dict("set_002", f3_tm=61.0),
        ]
        payload = {"primer_sets": sets}
        p = tmp_path / "primer_sets.json"
        p.write_text(json.dumps(payload))
        return p

    def test_basic_invocation(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "rt-check",
                "--input",
                str(self._make_primer_sets_json(tmp_path)),
                "--na-type",
                "rna",
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output

    def test_output_contains_set_ids(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["rt-check", "--input", str(self._make_primer_sets_json(tmp_path))],
            obj={},
        )
        assert "set_001" in result.output
        assert "set_002" in result.output

    def test_dna_target_shows_not_required(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "rt-check",
                "--input",
                str(self._make_primer_sets_json(tmp_path)),
                "--na-type",
                "dna",
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert "RT" in result.output or "DNA" in result.output

    def test_writes_csv_when_requested(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        out = tmp_path / "rt_check.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "rt-check",
                "--input",
                str(self._make_primer_sets_json(tmp_path)),
                "--out-csv",
                str(out),
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_top_n_limits_output(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "rt-check",
                "--input",
                str(self._make_primer_sets_json(tmp_path)),
                "--top-n",
                "1",
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert "set_001" in result.output
        assert "set_002" not in result.output

    def test_summary_line_appears(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["rt-check", "--input", str(self._make_primer_sets_json(tmp_path))],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert "optimized" in result.output.lower() or "RT-LAMP" in result.output
