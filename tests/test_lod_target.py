"""Tests for the reverse-LOD calculator (lod_target module).

All tests are pure-Python and require no external binaries.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from lamp_forge.lod_target import (
    extraction_sweep,
    required_effective_vol_ul,
    write_sweep_csv,
)

# ---------------------------------------------------------------------------
# required_effective_vol_ul
# ---------------------------------------------------------------------------


class TestRequiredEffectiveVol:
    def test_lod_95_formula(self) -> None:
        # lambda_LOD_95 = -ln(0.05) = 2.99573...
        # V_eff = lambda * 1000 / LOD_target
        # For LOD = 100 copies/mL: V_eff = 2996 / 100 = 29.96 uL
        v_eff = required_effective_vol_ul(100.0, detection_probability=0.95)
        expected = -math.log(0.05) * 1000.0 / 100.0
        assert abs(v_eff - expected) < 1e-9

    def test_lod_99(self) -> None:
        v_eff = required_effective_vol_ul(50.0, detection_probability=0.99)
        expected = -math.log(0.01) * 1000.0 / 50.0
        assert abs(v_eff - expected) < 1e-9

    def test_halving_lod_doubles_volume(self) -> None:
        v1 = required_effective_vol_ul(100.0)
        v2 = required_effective_vol_ul(50.0)
        assert abs(v2 / v1 - 2.0) < 1e-9

    def test_invalid_zero_lod(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            required_effective_vol_ul(0.0)

    def test_invalid_negative_lod(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            required_effective_vol_ul(-10.0)

    def test_invalid_probability_zero(self) -> None:
        with pytest.raises(ValueError):
            required_effective_vol_ul(100.0, detection_probability=0.0)

    def test_invalid_probability_one(self) -> None:
        with pytest.raises(ValueError):
            required_effective_vol_ul(100.0, detection_probability=1.0)


# ---------------------------------------------------------------------------
# extraction_sweep
# ---------------------------------------------------------------------------


class TestExtractionSweep:
    def test_returns_correct_number_of_points(self) -> None:
        pts = extraction_sweep(
            100.0,
            sample_volumes_ul=(200.0, 1000.0),
            efficiencies=(0.40, 0.60),
        )
        assert len(pts) == 4  # 2 volumes × 2 efficiencies

    def test_sorted_by_lod_ascending(self) -> None:
        pts = extraction_sweep(100.0)
        lods = [p.lod_copies_per_ml for p in pts]
        assert lods == sorted(lods)

    def test_achieves_target_flag(self) -> None:
        target = 100.0
        pts = extraction_sweep(
            target,
            sample_volumes_ul=(200.0, 1000.0, 2000.0),
            efficiencies=(0.30, 0.70),
        )
        for p in pts:
            if p.achieves_target:
                assert p.lod_copies_per_ml <= target
            else:
                assert p.lod_copies_per_ml > target

    def test_effective_vol_formula(self) -> None:
        # V_eff = V_sample * eta * (V_rxn / V_eluate)
        pts = extraction_sweep(
            100.0,
            sample_volumes_ul=(1000.0,),
            efficiencies=(0.50,),
            eluate_volume_ul=100.0,
            reaction_input_ul=5.0,
        )
        assert len(pts) == 1
        p = pts[0]
        expected_v_eff = 1000.0 * 0.50 * (5.0 / 100.0)  # = 25.0 uL
        assert abs(p.effective_vol_ul - expected_v_eff) < 1e-9

    def test_lod_formula_consistency(self) -> None:
        # The forward LOD formula should agree with the achieved LOD in sweep
        from lamp_forge.lod import ExtractionParams, estimate_lod

        pts = extraction_sweep(
            200.0,
            sample_volumes_ul=(1000.0,),
            efficiencies=(0.50,),
            eluate_volume_ul=100.0,
            reaction_input_ul=5.0,
        )
        p = pts[0]
        ext = ExtractionParams(
            sample_volume_ul=1000.0,
            extraction_efficiency=0.50,
            eluate_volume_ul=100.0,
            reaction_input_ul=5.0,
        )
        fwd = estimate_lod(ext, detection_probability=0.95)
        assert abs(p.lod_copies_per_ml - fwd.lod_copies_per_ml) < 1e-6

    def test_invalid_eluate_volume(self) -> None:
        with pytest.raises(ValueError, match="eluate_volume_ul"):
            extraction_sweep(100.0, eluate_volume_ul=0.0)

    def test_invalid_reaction_input_exceeds_eluate(self) -> None:
        with pytest.raises(ValueError):
            extraction_sweep(100.0, eluate_volume_ul=50.0, reaction_input_ul=60.0)

    def test_invalid_negative_sample_volume(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            extraction_sweep(100.0, sample_volumes_ul=(-100.0,))

    def test_invalid_efficiency_zero(self) -> None:
        with pytest.raises(ValueError, match=r"in \(0, 1\]"):
            extraction_sweep(100.0, efficiencies=(0.0,))

    def test_invalid_efficiency_over_one(self) -> None:
        with pytest.raises(ValueError, match=r"in \(0, 1\]"):
            extraction_sweep(100.0, efficiencies=(1.5,))

    def test_empty_sample_volumes(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            extraction_sweep(100.0, sample_volumes_ul=())

    def test_empty_efficiencies(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            extraction_sweep(100.0, efficiencies=())

    def test_default_parameters_present(self) -> None:
        # Default sweep: 4 volumes × 5 efficiencies = 20 points
        pts = extraction_sweep(100.0)
        assert len(pts) == 20  # 4 volumes * 5 efficiencies

    def test_all_fail_for_impossible_target(self) -> None:
        # Require LOD of 0.001 copies/mL -- nothing achieves this
        pts = extraction_sweep(0.001)
        assert not any(p.achieves_target for p in pts)

    def test_all_pass_for_easy_target(self) -> None:
        # Require LOD of 1e9 copies/mL -- everything achieves this
        pts = extraction_sweep(1e9)
        assert all(p.achieves_target for p in pts)

    def test_frozen_dataclass(self) -> None:
        pts = extraction_sweep(100.0, sample_volumes_ul=(200.0,), efficiencies=(0.50,))
        p = pts[0]
        with pytest.raises((AttributeError, TypeError)):
            p.achieves_target = not p.achieves_target  # type: ignore[misc]

    def test_point_fields(self) -> None:
        pts = extraction_sweep(
            100.0,
            sample_volumes_ul=(400.0,),
            efficiencies=(0.60,),
            eluate_volume_ul=80.0,
            reaction_input_ul=4.0,
        )
        p = pts[0]
        assert p.sample_volume_ul == 400.0
        assert p.extraction_efficiency == 0.60
        assert p.eluate_volume_ul == 80.0
        assert p.reaction_input_ul == 4.0


# ---------------------------------------------------------------------------
# write_sweep_csv
# ---------------------------------------------------------------------------


class TestWriteSweepCsv:
    def test_csv_header_and_rows(self, tmp_path: Path) -> None:
        pts = extraction_sweep(
            100.0,
            sample_volumes_ul=(200.0, 1000.0),
            efficiencies=(0.50,),
        )
        out = tmp_path / "sweep.csv"
        write_sweep_csv(pts, 100.0, out)

        assert out.exists()
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        assert len(rows) == 2
        expected_headers = {
            "sample_volume_ul",
            "extraction_efficiency",
            "eluate_volume_ul",
            "reaction_input_ul",
            "effective_vol_ul",
            "lod_copies_per_ml",
            "target_lod_copies_per_ml",
            "achieves_target",
        }
        assert expected_headers == set(rows[0].keys())

    def test_achieves_target_column_values(self, tmp_path: Path) -> None:
        pts = extraction_sweep(
            100.0,
            sample_volumes_ul=(200.0, 2000.0),
            efficiencies=(0.30,),  # 200*0.3*0.05=3 uL -> ~999 cp/mL FAIL
            # 2000*0.3*0.05=30 uL -> ~100 cp/mL PASS
        )
        out = tmp_path / "sweep.csv"
        write_sweep_csv(pts, 100.0, out)
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        # Sorted by LOD ascending: first row is 2000 uL (lower LOD = PASS)
        achieves = [r["achieves_target"] for r in rows]
        assert "yes" in achieves
        assert "no" in achieves

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        pts = extraction_sweep(100.0, sample_volumes_ul=(200.0,), efficiencies=(0.50,))
        out = tmp_path / "nested" / "dir" / "sweep.csv"
        write_sweep_csv(pts, 100.0, out)
        assert out.exists()

    def test_target_lod_written_to_every_row(self, tmp_path: Path) -> None:
        target = 250.0
        pts = extraction_sweep(target, sample_volumes_ul=(200.0,), efficiencies=(0.50,))
        out = tmp_path / "sweep.csv"
        write_sweep_csv(pts, target, out)
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        for row in rows:
            assert float(row["target_lod_copies_per_ml"]) == pytest.approx(target)
