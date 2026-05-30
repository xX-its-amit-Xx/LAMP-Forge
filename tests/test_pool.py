"""Tests for the multiplex primer-pool calculator.

All tests are pure Python — no external binaries.  Pool volumes are
deterministic arithmetic so exact floating-point comparisons with tolerances
are safe.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from lamp_forge.pool import (
    _ROLE_CONC_UM,
    _ROLE_ORDER,
    PoolEntry,
    PoolingParams,
    PoolPlan,
    build_pool_plan,
    write_pool_csv,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_ROLES_FULL = list(_ROLE_ORDER)  # ["F3", "B3", "FIP", "BIP", "LF", "LB"]
_ROLE_CONC_FULL = sum(_ROLE_CONC_UM.values())  # 44.0 uM for a 6-primer set


def _make_oligos(prefix: str = "") -> list[dict[str, str]]:
    """Return minimal oligo dicts for all six LAMP primer roles."""
    seqs = {
        "F3": "ACGTACGTACGTACGT",
        "B3": "TGCATGCATGCATGCA",
        "FIP": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "BIP": "TGCATGCATGCATGCATGCATGCATGCATGCATGCA",
        "LF": "GGCCGGCCGGCC",
        "LB": "CCGGCCGGCCGG",
    }
    return [{"role": r, "sequence": (prefix + seqs[r])[-len(seqs[r]) :]} for r in _ROLES_FULL]


def _make_selection(
    targets: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Build a minimal selection list as it appears in panel.json."""
    if targets is None:
        targets = ["SRB"]
    return [
        {
            "target": t,
            "set_id": "region_01_set_001",
            "composite_score": 0.9,
            "oligos": _make_oligos(),
        }
        for t in targets
    ]


# ---------------------------------------------------------------------------
# PoolingParams
# ---------------------------------------------------------------------------


class TestPoolingParams:
    def test_defaults(self) -> None:
        p = PoolingParams()
        assert p.stock_conc_um == 100.0
        assert p.total_pool_volume_ul == 500.0

    def test_custom_values(self) -> None:
        p = PoolingParams(stock_conc_um=200.0, total_pool_volume_ul=1000.0)
        assert p.stock_conc_um == 200.0
        assert p.total_pool_volume_ul == 1000.0

    def test_rejects_zero_stock_conc(self) -> None:
        with pytest.raises(ValueError, match="stock_conc_um"):
            PoolingParams(stock_conc_um=0.0)

    def test_rejects_negative_stock_conc(self) -> None:
        with pytest.raises(ValueError, match="stock_conc_um"):
            PoolingParams(stock_conc_um=-10.0)

    def test_rejects_zero_total_volume(self) -> None:
        with pytest.raises(ValueError, match="total_pool_volume_ul"):
            PoolingParams(total_pool_volume_ul=0.0)

    def test_rejects_negative_total_volume(self) -> None:
        with pytest.raises(ValueError, match="total_pool_volume_ul"):
            PoolingParams(total_pool_volume_ul=-1.0)

    def test_very_high_stock_conc_valid(self) -> None:
        p = PoolingParams(stock_conc_um=500.0)
        assert p.stock_conc_um == 500.0


# ---------------------------------------------------------------------------
# build_pool_plan — basic construction
# ---------------------------------------------------------------------------


class TestBuildPoolPlanBasic:
    def test_returns_pool_plan(self) -> None:
        selection = _make_selection()
        plan = build_pool_plan(selection)
        assert isinstance(plan, PoolPlan)

    def test_n_targets_correct(self) -> None:
        selection = _make_selection(["SRB", "MCR"])
        plan = build_pool_plan(selection, PoolingParams(stock_conc_um=200.0))
        assert plan.n_targets == 2

    def test_n_primers_six_per_target(self) -> None:
        selection = _make_selection()
        plan = build_pool_plan(selection)
        assert plan.n_primers == 6

    def test_total_volume_stored(self) -> None:
        params = PoolingParams(total_pool_volume_ul=800.0)
        plan = build_pool_plan(_make_selection(), params)
        assert plan.total_pool_volume_ul == 800.0

    def test_all_entries_are_pool_entry(self) -> None:
        plan = build_pool_plan(_make_selection())
        assert all(isinstance(e, PoolEntry) for e in plan.entries)

    def test_raises_on_empty_selection(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            build_pool_plan([])

    def test_default_params_used_when_none(self) -> None:
        plan = build_pool_plan(_make_selection(), params=None)
        assert plan.total_pool_volume_ul == 500.0


# ---------------------------------------------------------------------------
# build_pool_plan — volume arithmetic
# ---------------------------------------------------------------------------


class TestBuildPoolPlanVolumes:
    def test_vol_stock_formula_fip(self) -> None:
        # vol = target_conc * pool_vol / stock_conc
        # FIP: 16 uM * 500 / 100 = 80 uL
        params = PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0)
        plan = build_pool_plan(_make_selection(["SRB"]), params)
        fip_entry = next(e for e in plan.entries if e.primer_role == "FIP")
        assert abs(fip_entry.vol_stock_ul - 80.0) < 1e-9

    def test_vol_stock_formula_f3(self) -> None:
        # F3: 2 uM * 500 / 100 = 10 uL
        params = PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0)
        plan = build_pool_plan(_make_selection(["SRB"]), params)
        f3_entry = next(e for e in plan.entries if e.primer_role == "F3")
        assert abs(f3_entry.vol_stock_ul - 10.0) < 1e-9

    def test_water_volume_balances_total(self) -> None:
        params = PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0)
        plan = build_pool_plan(_make_selection(["SRB"]), params)
        total = sum(e.vol_stock_ul for e in plan.entries) + plan.water_volume_ul
        assert abs(total - plan.total_pool_volume_ul) < 1e-9

    def test_water_volume_positive_single_target(self) -> None:
        # 1 target: 44 uM sum < 100 uM stock -> water > 0
        plan = build_pool_plan(_make_selection(["SRB"]))
        assert plan.water_volume_ul > 0.0

    def test_water_volume_single_target_correct(self) -> None:
        # stock = 100 uM, vol = 500 uL, sum_conc = 44 uM
        # water = (1 - 44/100) * 500 = 280 uL
        params = PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0)
        plan = build_pool_plan(_make_selection(["SRB"]), params)
        expected_water = 500.0 * (1.0 - _ROLE_CONC_FULL / 100.0)
        assert abs(plan.water_volume_ul - expected_water) < 1e-9

    def test_two_targets_100um_stock_positive_water(self) -> None:
        # 2 targets: 88 uM sum < 100 uM stock -> water > 0
        params = PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0)
        plan = build_pool_plan(_make_selection(["SRB", "MCR"]), params)
        assert plan.water_volume_ul > 0.0

    def test_higher_stock_conc_gives_more_water(self) -> None:
        selection = _make_selection(["SRB"])
        plan_100 = build_pool_plan(selection, PoolingParams(stock_conc_um=100.0))
        plan_200 = build_pool_plan(selection, PoolingParams(stock_conc_um=200.0))
        assert plan_200.water_volume_ul > plan_100.water_volume_ul

    def test_larger_total_volume_scales_stock_volumes(self) -> None:
        s = _make_selection(["SRB"])
        plan_500 = build_pool_plan(s, PoolingParams(total_pool_volume_ul=500.0))
        plan_1000 = build_pool_plan(s, PoolingParams(total_pool_volume_ul=1000.0))
        # All primer volumes should be exactly double.
        for e500, e1000 in zip(plan_500.entries, plan_1000.entries, strict=True):
            assert abs(e1000.vol_stock_ul - 2.0 * e500.vol_stock_ul) < 1e-9

    def test_three_targets_100um_stock_raises_infeasible(self) -> None:
        # 3 * 44 = 132 uM > 100 uM stock -> infeasible
        params = PoolingParams(stock_conc_um=100.0)
        with pytest.raises(ValueError, match="too low"):
            build_pool_plan(_make_selection(["SRB", "MCR", "16S"]), params)

    def test_error_message_names_min_stock(self) -> None:
        params = PoolingParams(stock_conc_um=100.0)
        with pytest.raises(ValueError, match="132"):
            build_pool_plan(_make_selection(["SRB", "MCR", "16S"]), params)

    def test_three_targets_200um_stock_feasible(self) -> None:
        params = PoolingParams(stock_conc_um=200.0)
        plan = build_pool_plan(_make_selection(["SRB", "MCR", "16S"]), params)
        assert plan.n_targets == 3
        assert plan.water_volume_ul >= 0.0


# ---------------------------------------------------------------------------
# build_pool_plan — entry ordering and content
# ---------------------------------------------------------------------------


class TestBuildPoolPlanEntries:
    def test_role_order_within_target(self) -> None:
        plan = build_pool_plan(_make_selection(["SRB"]))
        roles = [e.primer_role for e in plan.entries]
        assert roles == list(_ROLE_ORDER)

    def test_targets_grouped_together(self) -> None:
        params = PoolingParams(stock_conc_um=200.0)
        plan = build_pool_plan(_make_selection(["SRB", "MCR"]), params)
        srb_roles = [e.primer_role for e in plan.entries if e.target_label == "SRB"]
        mcr_roles = [e.primer_role for e in plan.entries if e.target_label == "MCR"]
        assert srb_roles == list(_ROLE_ORDER)
        assert mcr_roles == list(_ROLE_ORDER)

    def test_primer_name_contains_target_and_role(self) -> None:
        plan = build_pool_plan(_make_selection(["SRB"]))
        for e in plan.entries:
            assert e.target_label in e.primer_name
            assert e.primer_role in e.primer_name

    def test_sequence_uppercased(self) -> None:
        selection = [
            {
                "target": "SRB",
                "set_id": "s001",
                "composite_score": 0.9,
                "oligos": [
                    {"role": "F3", "sequence": "acgtacgt"},
                    {"role": "B3", "sequence": "tgcatgca"},
                    {"role": "FIP", "sequence": "acgtacgtacgtacgtacgt"},
                    {"role": "BIP", "sequence": "tgcatgcatgcatgcatgca"},
                ],
            }
        ]
        plan = build_pool_plan(selection)
        for e in plan.entries:
            assert e.sequence == e.sequence.upper()

    def test_missing_loop_primers_skipped(self) -> None:
        selection = [
            {
                "target": "SRB",
                "set_id": "s001",
                "composite_score": 0.9,
                "oligos": [
                    {"role": "F3", "sequence": "ACGTACGT"},
                    {"role": "B3", "sequence": "TGCATGCA"},
                    {"role": "FIP", "sequence": "ACGTACGTACGTACGTACGT"},
                    {"role": "BIP", "sequence": "TGCATGCATGCATGCATGCA"},
                ],
            }
        ]
        plan = build_pool_plan(selection)
        assert plan.n_primers == 4
        roles = {e.primer_role for e in plan.entries}
        assert "LF" not in roles
        assert "LB" not in roles

    def test_target_conc_matches_role_spec(self) -> None:
        plan = build_pool_plan(_make_selection(["SRB"]))
        for e in plan.entries:
            assert abs(e.target_conc_um - _ROLE_CONC_UM[e.primer_role]) < 1e-9

    def test_stock_conc_stored_in_entries(self) -> None:
        params = PoolingParams(stock_conc_um=200.0)
        plan = build_pool_plan(_make_selection(["SRB"]), params)
        for e in plan.entries:
            assert e.stock_conc_um == 200.0


# ---------------------------------------------------------------------------
# write_pool_csv
# ---------------------------------------------------------------------------


class TestWritePoolCsv:
    def _plan(self, targets: list[str] | None = None, stock: float = 100.0) -> PoolPlan:
        t = targets or ["SRB"]
        params = PoolingParams(stock_conc_um=stock)
        return build_pool_plan(_make_selection(t), params)

    def test_creates_file(self, tmp_path: Path) -> None:
        out = tmp_path / "pool.csv"
        write_pool_csv(self._plan(), out)
        assert out.exists()

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "sub" / "dir" / "pool.csv"
        write_pool_csv(self._plan(), out)
        assert out.exists()

    def test_header_columns(self, tmp_path: Path) -> None:
        out = tmp_path / "pool.csv"
        write_pool_csv(self._plan(), out)
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert header == [
            "target",
            "role",
            "primer_name",
            "sequence",
            "stock_conc_um",
            "target_conc_um",
            "vol_stock_ul",
        ]

    def test_data_rows_equal_n_primers_plus_water(self, tmp_path: Path) -> None:
        plan = self._plan()
        out = tmp_path / "pool.csv"
        write_pool_csv(plan, out)
        with out.open() as fh:
            rows = list(csv.reader(fh))
        # header + n_primers + 1 WATER row
        assert len(rows) == 1 + plan.n_primers + 1

    def test_water_row_last(self, tmp_path: Path) -> None:
        out = tmp_path / "pool.csv"
        write_pool_csv(self._plan(), out)
        with out.open() as fh:
            rows = list(csv.reader(fh))
        assert rows[-1][0] == "WATER"

    def test_vol_stock_numeric(self, tmp_path: Path) -> None:
        out = tmp_path / "pool.csv"
        write_pool_csv(self._plan(), out)
        with out.open() as fh:
            data = list(csv.DictReader(fh))
        non_water = [r for r in data if r["target"] != "WATER"]
        for row in non_water:
            val = float(row["vol_stock_ul"])
            assert val > 0.0

    def test_water_volume_numeric_in_csv(self, tmp_path: Path) -> None:
        plan = self._plan()
        out = tmp_path / "pool.csv"
        write_pool_csv(plan, out)
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        water_row = next(r for r in rows if r["target"] == "WATER")
        csv_water = float(water_row["vol_stock_ul"])
        assert abs(csv_water - plan.water_volume_ul) < 1e-3

    def test_sequences_uppercase_in_csv(self, tmp_path: Path) -> None:
        out = tmp_path / "pool.csv"
        write_pool_csv(self._plan(), out)
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        for row in rows:
            seq = row["sequence"]
            if seq:
                assert seq == seq.upper()

    def test_two_target_csv_has_correct_row_count(self, tmp_path: Path) -> None:
        plan = self._plan(targets=["SRB", "MCR"], stock=200.0)
        out = tmp_path / "pool.csv"
        write_pool_csv(plan, out)
        with out.open() as fh:
            rows = list(csv.reader(fh))
        assert len(rows) == 1 + 12 + 1  # header + 2*6 primers + water


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------


class TestPoolCli:
    def _write_panel_json(self, path: Path, targets: list[str]) -> None:
        """Write a minimal panel.json fixture."""
        selection = []
        seqs = {
            "F3": "ACGTACGTACGTACGT",
            "B3": "TGCATGCATGCATGCA",
            "FIP": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "BIP": "TGCATGCATGCATGCATGCATGCATGCATGCATGCA",
            "LF": "GGCCGGCCGGCC",
            "LB": "CCGGCCGGCCGG",
        }
        for t in targets:
            selection.append(
                {
                    "target": t,
                    "set_id": "region_01_set_001",
                    "composite_score": 0.9,
                    "oligos": [{"role": r, "sequence": seqs[r]} for r in _ROLE_ORDER],
                }
            )
        payload = {
            "version": "0.1.0",
            "is_clean": True,
            "worst_inter_assay_dg": -3.5,
            "selection": selection,
        }
        path.write_text(json.dumps(payload), encoding="utf-8")

    def test_basic_invocation_single_target(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB"])
        out = tmp_path / "pool.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["pool", "--panel", str(panel), "--out", str(out)],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_output_mentions_targets_and_primers(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB"])
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["pool", "--panel", str(panel)],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert "target" in result.output.lower() or "1" in result.output
        assert "primer" in result.output.lower() or "6" in result.output

    def test_default_output_path(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB"])
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["pool", "--panel", str(panel)],
            obj={},
        )
        assert result.exit_code == 0, result.output
        default_out = tmp_path / "pool_sheet.csv"
        assert default_out.exists()

    def test_two_targets_200um_stock(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB", "MCR"])
        out = tmp_path / "pool.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "pool",
                "--panel",
                str(panel),
                "--stock-conc",
                "200",
                "--out",
                str(out),
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        with out.open() as fh:
            rows = list(csv.reader(fh))
        # header + 12 primer rows + 1 water row
        assert len(rows) == 14

    def test_infeasible_stock_exits_nonzero(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB", "MCR", "16S"])
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "pool",
                "--panel",
                str(panel),
                "--stock-conc",
                "100",
            ],
            obj={},
        )
        assert result.exit_code != 0

    def test_stock_conc_override(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB"])
        out = tmp_path / "pool.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "pool",
                "--panel",
                str(panel),
                "--stock-conc",
                "200",
                "--out",
                str(out),
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        non_water = [r for r in rows if r["target"] != "WATER"]
        for row in non_water:
            assert row["stock_conc_um"] == "200.0"

    def test_total_volume_override(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        panel = tmp_path / "panel.json"
        self._write_panel_json(panel, ["SRB"])
        out = tmp_path / "pool.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "pool",
                "--panel",
                str(panel),
                "--total-volume",
                "1000",
                "--out",
                str(out),
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        # FIP should be 16 * 1000 / 100 = 160 uL
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        fip_row = next(r for r in rows if r["role"] == "FIP")
        assert abs(float(fip_row["vol_stock_ul"]) - 160.0) < 0.01
