"""Tests for the time-to-positive (TTP) estimator.

All tests are pure Python -- no external binaries required.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from lamp_forge.ttp import (
    TtpParams,
    TtpPoint,
    TtpPreset,
    TtpTable,
    ttp_at_copies,
    ttp_from_sample_concentration,
    ttp_table,
    write_ttp_csv,
)

# ---------------------------------------------------------------------------
# TtpParams construction
# ---------------------------------------------------------------------------


class TestTtpParams:
    def test_defaults(self) -> None:
        p = TtpParams()
        assert p.ttp_one_copy_min == 55.0
        assert p.slope_min_per_decade == 6.0
        assert p.device_window_min == 60.0
        assert p.min_ttp_min == 10.0

    def test_invalid_ttp_one_copy_zero(self) -> None:
        with pytest.raises(ValueError, match="ttp_one_copy_min"):
            TtpParams(ttp_one_copy_min=0.0)

    def test_invalid_ttp_one_copy_negative(self) -> None:
        with pytest.raises(ValueError, match="ttp_one_copy_min"):
            TtpParams(ttp_one_copy_min=-1.0)

    def test_invalid_slope_zero(self) -> None:
        with pytest.raises(ValueError, match="slope_min_per_decade"):
            TtpParams(slope_min_per_decade=0.0)

    def test_invalid_device_window_zero(self) -> None:
        with pytest.raises(ValueError, match="device_window_min"):
            TtpParams(device_window_min=0.0)

    def test_invalid_min_ttp_negative(self) -> None:
        with pytest.raises(ValueError, match="min_ttp_min"):
            TtpParams(min_ttp_min=-1.0)

    def test_from_preset_dna_lamp(self) -> None:
        p = TtpParams.from_preset(TtpPreset.DNA_LAMP)
        assert p.ttp_one_copy_min == 55.0
        assert p.slope_min_per_decade == 6.0
        assert p.device_window_min == 60.0
        assert p.min_ttp_min == 10.0

    def test_from_preset_rt_lamp(self) -> None:
        p = TtpParams.from_preset(TtpPreset.RT_LAMP)
        assert p.ttp_one_copy_min == 60.0
        assert p.min_ttp_min == 12.0

    def test_from_preset_fast_lamp(self) -> None:
        p = TtpParams.from_preset(TtpPreset.FAST_LAMP)
        assert p.ttp_one_copy_min == 50.0
        assert p.min_ttp_min == 8.0

    def test_from_preset_custom_device_window(self) -> None:
        p = TtpParams.from_preset(TtpPreset.DNA_LAMP, device_window_min=45.0)
        assert p.device_window_min == 45.0

    def test_frozen(self) -> None:
        p = TtpParams()
        with pytest.raises(AttributeError):
            p.ttp_one_copy_min = 99.0  # type: ignore[misc]


# ---------------------------------------------------------------------------
# ttp_at_copies
# ---------------------------------------------------------------------------


class TestTtpAtCopies:
    def test_one_copy_returns_ttp_one(self) -> None:
        params = TtpParams(ttp_one_copy_min=55.0, slope_min_per_decade=6.0, min_ttp_min=0.0)
        pt = ttp_at_copies(1.0, params)
        assert pt.copies_per_reaction == 1.0
        assert pt.ttp_minutes == pytest.approx(55.0, abs=0.01)

    def test_ten_copies_reduces_by_slope(self) -> None:
        params = TtpParams(ttp_one_copy_min=55.0, slope_min_per_decade=6.0, min_ttp_min=0.0)
        pt = ttp_at_copies(10.0, params)
        assert pt.ttp_minutes == pytest.approx(55.0 - 6.0, abs=0.01)

    def test_hundred_copies_two_decades(self) -> None:
        params = TtpParams(ttp_one_copy_min=55.0, slope_min_per_decade=6.0, min_ttp_min=0.0)
        pt = ttp_at_copies(100.0, params)
        assert pt.ttp_minutes == pytest.approx(55.0 - 12.0, abs=0.01)

    def test_min_ttp_floor_applied(self) -> None:
        params = TtpParams(ttp_one_copy_min=20.0, slope_min_per_decade=6.0, min_ttp_min=15.0)
        pt = ttp_at_copies(1e6, params)
        assert pt.ttp_minutes == pytest.approx(15.0, abs=0.01)

    def test_in_window_true_when_ttp_within_limit(self) -> None:
        params = TtpParams(ttp_one_copy_min=40.0, slope_min_per_decade=6.0, device_window_min=60.0)
        pt = ttp_at_copies(10.0, params)
        assert pt.in_device_window is True

    def test_in_window_false_when_ttp_exceeds_limit(self) -> None:
        params = TtpParams(ttp_one_copy_min=70.0, slope_min_per_decade=6.0, device_window_min=60.0)
        pt = ttp_at_copies(1.0, params)
        assert pt.in_device_window is False

    def test_zero_copies_raises(self) -> None:
        with pytest.raises(ValueError, match="copies_per_reaction"):
            ttp_at_copies(0.0, TtpParams())

    def test_negative_copies_raises(self) -> None:
        with pytest.raises(ValueError, match="copies_per_reaction"):
            ttp_at_copies(-5.0, TtpParams())

    def test_fractional_copies_treated_as_one(self) -> None:
        params = TtpParams(ttp_one_copy_min=55.0, slope_min_per_decade=6.0, min_ttp_min=0.0)
        pt_half = ttp_at_copies(0.5, params)
        pt_one = ttp_at_copies(1.0, params)
        assert pt_half.ttp_minutes == pytest.approx(pt_one.ttp_minutes, abs=0.01)

    def test_returns_ttp_point(self) -> None:
        pt = ttp_at_copies(100.0, TtpParams())
        assert isinstance(pt, TtpPoint)

    def test_ttp_monotonically_decreasing(self) -> None:
        params = TtpParams(min_ttp_min=0.0)
        prev_ttp = math.inf
        for copies in [1, 10, 100, 1000, 10000]:
            pt = ttp_at_copies(float(copies), params)
            assert pt.ttp_minutes < prev_ttp
            prev_ttp = pt.ttp_minutes


# ---------------------------------------------------------------------------
# ttp_table
# ---------------------------------------------------------------------------


class TestTtpTable:
    def test_returns_ttp_table(self) -> None:
        result = ttp_table(TtpParams())
        assert isinstance(result, TtpTable)

    def test_default_has_13_points(self) -> None:
        result = ttp_table(TtpParams())
        assert len(result.points) == 13

    def test_points_ordered_ascending(self) -> None:
        result = ttp_table(TtpParams())
        copies = [pt.copies_per_reaction for pt in result.points]
        assert copies == sorted(copies)

    def test_custom_n_points(self) -> None:
        result = ttp_table(TtpParams(), n_points=5)
        assert len(result.points) == 5

    def test_min_detectable_copies_set_correctly(self) -> None:
        params = TtpParams(ttp_one_copy_min=70.0, slope_min_per_decade=6.0, device_window_min=60.0)
        result = ttp_table(params, copies_range=(1.0, 1e8))
        if result.min_detectable_copies is not None:
            pt = ttp_at_copies(result.min_detectable_copies, params)
            assert pt.in_device_window is True
            assert result.points[0].in_device_window is False

    def test_min_detectable_copies_none_when_all_too_slow(self) -> None:
        params = TtpParams(ttp_one_copy_min=200.0, device_window_min=10.0, min_ttp_min=0.0)
        result = ttp_table(params, copies_range=(1.0, 100.0), n_points=5)
        assert result.min_detectable_copies is None

    def test_min_detectable_copies_set_when_all_in_window(self) -> None:
        params = TtpParams(ttp_one_copy_min=10.0, device_window_min=60.0, min_ttp_min=0.0)
        result = ttp_table(params, copies_range=(1.0, 100.0), n_points=5)
        assert result.min_detectable_copies is not None
        assert result.min_detectable_copies == pytest.approx(1.0, rel=0.01)

    def test_invalid_range_lo_gt_hi(self) -> None:
        with pytest.raises(ValueError, match="copies_range"):
            ttp_table(TtpParams(), copies_range=(100.0, 1.0))

    def test_invalid_range_zero_lo(self) -> None:
        with pytest.raises(ValueError, match="copies_range"):
            ttp_table(TtpParams(), copies_range=(0.0, 100.0))

    def test_invalid_n_points_one(self) -> None:
        with pytest.raises(ValueError, match="n_points"):
            ttp_table(TtpParams(), n_points=1)

    def test_single_point_range(self) -> None:
        result = ttp_table(TtpParams(), copies_range=(10.0, 10.0), n_points=2)
        assert len(result.points) == 2
        for pt in result.points:
            assert pt.copies_per_reaction == pytest.approx(10.0, rel=0.01)


# ---------------------------------------------------------------------------
# ttp_from_sample_concentration
# ---------------------------------------------------------------------------


class TestTtpFromSampleConcentration:
    def _make_extraction(self) -> object:
        from lamp_forge.lod import ExtractionParams

        return ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=2.5,
        )

    def test_basic_call_returns_ttp_point(self) -> None:
        ext = self._make_extraction()
        from lamp_forge.lod import ExtractionParams

        assert isinstance(ext, ExtractionParams)
        pt = ttp_from_sample_concentration(1e5, ext, TtpParams())  # type: ignore[arg-type]
        assert isinstance(pt, TtpPoint)

    def test_higher_concentration_gives_shorter_ttp(self) -> None:
        ext = self._make_extraction()
        pt_low = ttp_from_sample_concentration(100.0, ext, TtpParams())  # type: ignore[arg-type]
        pt_high = ttp_from_sample_concentration(1e5, ext, TtpParams())  # type: ignore[arg-type]
        assert pt_high.ttp_minutes < pt_low.ttp_minutes

    def test_zero_concentration_raises(self) -> None:
        ext = self._make_extraction()
        with pytest.raises(ValueError, match="conc_copies_per_ml"):
            ttp_from_sample_concentration(0.0, ext, TtpParams())  # type: ignore[arg-type]

    def test_negative_concentration_raises(self) -> None:
        ext = self._make_extraction()
        with pytest.raises(ValueError, match="conc_copies_per_ml"):
            ttp_from_sample_concentration(-1.0, ext, TtpParams())  # type: ignore[arg-type]

    def test_copies_in_rxn_consistent_with_extraction(self) -> None:
        from lamp_forge.lod import ExtractionParams

        ext = ExtractionParams(
            sample_volume_ul=1000.0,
            extraction_efficiency=1.0,
            eluate_volume_ul=100.0,
            reaction_input_ul=5.0,
        )
        # conv = 1000 * 1.0 * 5/100 = 50 uL effective sample
        # 1000 copies/mL = 1 copy/uL -> 50 copies in rxn
        pt = ttp_from_sample_concentration(1000.0, ext, TtpParams(min_ttp_min=0.0))
        expected_ttp = 55.0 - 6.0 * math.log10(50.0)
        assert pt.ttp_minutes == pytest.approx(expected_ttp, abs=0.1)


# ---------------------------------------------------------------------------
# write_ttp_csv
# ---------------------------------------------------------------------------


class TestWriteTtpCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        table = ttp_table(TtpParams(), n_points=5)
        out = tmp_path / "ttp.csv"
        write_ttp_csv(table, out)
        assert out.exists()

    def test_csv_has_header_row(self, tmp_path: Path) -> None:
        table = ttp_table(TtpParams(), n_points=3)
        out = tmp_path / "ttp.csv"
        write_ttp_csv(table, out)
        with out.open(newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert "copies_per_reaction" in header
        assert "ttp_minutes" in header
        assert "in_device_window" in header

    def test_csv_row_count_matches_points(self, tmp_path: Path) -> None:
        table = ttp_table(TtpParams(), n_points=7)
        out = tmp_path / "ttp.csv"
        write_ttp_csv(table, out)
        with out.open(newline="") as fh:
            rows = list(csv.reader(fh))
        # 1 header + n_points data rows
        assert len(rows) == 8

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        table = ttp_table(TtpParams(), n_points=2)
        out = tmp_path / "subdir" / "deep" / "ttp.csv"
        write_ttp_csv(table, out)
        assert out.exists()

    def test_in_device_window_written_as_lowercase_bool(self, tmp_path: Path) -> None:
        table = ttp_table(TtpParams(), n_points=3)
        out = tmp_path / "ttp.csv"
        write_ttp_csv(table, out)
        with out.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        for row in rows:
            assert row["in_device_window"] in ("true", "false")


# ---------------------------------------------------------------------------
# Preset coverage -- all presets must be in the lookup dict
# ---------------------------------------------------------------------------


class TestPresetCoverage:
    def test_all_presets_build_params(self) -> None:
        for preset in TtpPreset:
            params = TtpParams.from_preset(preset)
            assert isinstance(params, TtpParams)

    def test_rt_lamp_ttp_slower_than_dna_lamp_at_low_copies(self) -> None:
        dna = TtpParams.from_preset(TtpPreset.DNA_LAMP)
        rt = TtpParams.from_preset(TtpPreset.RT_LAMP)
        pt_dna = ttp_at_copies(1.0, dna)
        pt_rt = ttp_at_copies(1.0, rt)
        assert pt_rt.ttp_minutes > pt_dna.ttp_minutes

    def test_fast_lamp_ttp_faster_than_dna_lamp(self) -> None:
        dna = TtpParams.from_preset(TtpPreset.DNA_LAMP)
        fast = TtpParams.from_preset(TtpPreset.FAST_LAMP)
        pt_dna = ttp_at_copies(1.0, dna)
        pt_fast = ttp_at_copies(1.0, fast)
        assert pt_fast.ttp_minutes < pt_dna.ttp_minutes
