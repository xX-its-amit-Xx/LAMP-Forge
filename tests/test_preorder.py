"""Tests for the combined pre-order readiness check (lamp_forge.preorder).

All tests are pure Python -- no external binaries required.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import pytest

from lamp_forge.lod import ExtractionParams
from lamp_forge.preorder import (
    PreorderResult,
    PreorderStatus,
    run_preorder,
    write_preorder_csv,
)
from lamp_forge.rt_lamp import RtLampParams, TargetNucleicAcid
from lamp_forge.ttp import TtpParams, TtpPreset

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_extraction(
    sample_volume_ul: float = 1000.0,
    extraction_efficiency: float = 0.50,
    eluate_volume_ul: float = 100.0,
    reaction_input_ul: float = 5.0,
) -> ExtractionParams:
    return ExtractionParams(
        sample_volume_ul=sample_volume_ul,
        extraction_efficiency=extraction_efficiency,
        eluate_volume_ul=eluate_volume_ul,
        reaction_input_ul=reaction_input_ul,
    )


def _make_primer_set(set_id: str, tm: float = 64.0) -> dict[str, Any]:
    """Minimal primer set dict as returned by primer_sets.json."""
    primer = {"role": "F3", "tm": tm, "sequence": "ACGTACGT"}
    return {
        "set_id": set_id,
        "f3": {**primer, "role": "F3"},
        "b3": {**primer, "role": "B3"},
        "fip": {**primer, "role": "FIP"},
        "bip": {**primer, "role": "BIP"},
        "lf": {**primer, "role": "LF"},
        "lb": {**primer, "role": "LB"},
    }


def _make_sets(n: int = 5, tm: float = 64.0) -> list[dict[str, Any]]:
    return [_make_primer_set(f"region_01_set_{i:03d}", tm=tm) for i in range(1, n + 1)]


def _dna_lamp_params(device_window: float = 60.0) -> TtpParams:
    return TtpParams.from_preset(TtpPreset.DNA_LAMP, device_window_min=device_window)


def _rt_lamp_params(device_window: float = 60.0) -> TtpParams:
    return TtpParams.from_preset(TtpPreset.RT_LAMP, device_window_min=device_window)


# ---------------------------------------------------------------------------
# LOD_95 accuracy
# ---------------------------------------------------------------------------


class TestLod95:
    def test_lod_95_copies_per_rxn_approx_poisson(self) -> None:
        extraction = _make_extraction()
        result = run_preorder(_make_sets(), extraction, _dna_lamp_params())
        expected_lod_rxn = -math.log(1.0 - 0.95)
        assert result.lod_95_copies_per_rxn == pytest.approx(expected_lod_rxn, rel=1e-6)

    def test_lod_95_copies_per_ml_consistent_with_extraction(self) -> None:
        extraction = _make_extraction(
            sample_volume_ul=1000.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=100.0,
            reaction_input_ul=5.0,
        )
        result = run_preorder(_make_sets(), extraction, _dna_lamp_params())
        # effective_ul = 1000 * 0.5 * 5/100 = 25 uL
        # LOD_ml = LOD_rxn / 25 * 1000
        expected_ml = result.lod_95_copies_per_rxn / 25.0 * 1000.0
        assert result.lod_95_copies_per_ml == pytest.approx(expected_ml, rel=1e-6)

    def test_larger_sample_gives_lower_lod_ml(self) -> None:
        small = run_preorder(
            _make_sets(), _make_extraction(sample_volume_ul=200.0), _dna_lamp_params()
        )
        large = run_preorder(
            _make_sets(), _make_extraction(sample_volume_ul=1000.0), _dna_lamp_params()
        )
        assert large.lod_95_copies_per_ml < small.lod_95_copies_per_ml


# ---------------------------------------------------------------------------
# TTP at LOD_95
# ---------------------------------------------------------------------------


class TestTtpAtLod:
    def test_ttp_at_lod_positive(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        assert result.ttp_at_lod_min > 0.0

    def test_ttp_less_than_device_window_for_standard_params(self) -> None:
        # 1000 uL * 50% * 5/100 = 25 uL effective; LOD_95 ~ 3 copies/rxn
        # TTP(3 copies, dna-lamp) = 55 - 6*log10(3) ~ 52 min < 60 min
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params(60.0))
        assert result.ttp_in_window is True

    def test_very_small_extraction_yields_no_go(self) -> None:
        # tiny extraction -> LOD_95 ~= 1 copy/rxn -> TTP ~= 55 min
        # set device window to 30 min to force NO_GO
        tiny = ExtractionParams(
            sample_volume_ul=1.0,
            extraction_efficiency=0.01,
            eluate_volume_ul=100.0,
            reaction_input_ul=0.1,
        )
        result = run_preorder(_make_sets(), tiny, _dna_lamp_params(device_window=30.0))
        assert result.status is PreorderStatus.NO_GO
        assert result.ttp_in_window is False

    def test_ttp_at_lod_monotonic_with_sample_volume(self) -> None:
        volumes = [200.0, 500.0, 1000.0, 2000.0]
        ttps = []
        for vol in volumes:
            ext = _make_extraction(sample_volume_ul=vol)
            r = run_preorder(_make_sets(), ext, _dna_lamp_params())
            ttps.append(r.ttp_at_lod_min)
        # More sample = more copies/rxn = lower TTP
        for i in range(len(ttps) - 1):
            assert ttps[i] >= ttps[i + 1]


# ---------------------------------------------------------------------------
# Status logic
# ---------------------------------------------------------------------------


class TestStatusLogic:
    def test_go_status_dna_target_good_extraction(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        assert result.status is PreorderStatus.GO

    def test_no_go_status_when_window_too_short(self) -> None:
        # Very short device window to force NO_GO
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params(device_window=1.0))
        assert result.status is PreorderStatus.NO_GO

    def test_go_status_rna_all_primers_optimised(self) -> None:
        # tm=64 > 63 degC RT floor -- all sets pass RT check.
        # Use a generous 80-min window so TTP (~57 min) is well inside
        # and does not trigger the near-limit WARNING path.
        result = run_preorder(
            _make_sets(tm=64.0),
            _make_extraction(),
            TtpParams.from_preset(TtpPreset.RT_LAMP, device_window_min=80.0),
            target_na=TargetNucleicAcid.RNA,
            rt_params=RtLampParams(rt_min_tm=63.0),
        )
        assert result.status is PreorderStatus.GO
        assert result.n_sets_rt_ok == result.n_sets_assessed

    def test_warning_status_rna_some_primers_below_rt_floor(self) -> None:
        # tm=60 < 63 degC -- all core primers are below RT floor
        result = run_preorder(
            _make_sets(tm=60.0),
            _make_extraction(),
            _rt_lamp_params(),
            target_na=TargetNucleicAcid.RNA,
            rt_params=RtLampParams(rt_min_tm=63.0),
        )
        assert result.status is PreorderStatus.WARNING
        assert result.n_sets_rt_ok == 0

    def test_warning_status_ttp_near_window_boundary(self) -> None:
        # LOD_95 copies/rxn ~ 2.996.
        # TTP(3 copies, dna-lamp) = 55 - 6*log10(3) ~ 52.14 min.
        # For TTP/window > 0.95 (near-limit), need window < 52.14/0.95 ~ 54.9.
        # Use device_window=54 -> 52.14/54 ~ 0.966 > 0.95 => WARNING.
        result = run_preorder(
            _make_sets(), _make_extraction(), _dna_lamp_params(device_window=54.0)
        )
        assert result.ttp_in_window is True
        assert result.status is PreorderStatus.WARNING

    def test_reasons_not_empty(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        assert len(result.reasons) > 0

    def test_no_go_has_reason_about_device_window(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params(device_window=1.0))
        assert any("device window" in r.lower() for r in result.reasons)


# ---------------------------------------------------------------------------
# Target NA handling
# ---------------------------------------------------------------------------


class TestTargetNaHandling:
    def test_dna_target_n_sets_rt_ok_is_minus_one(self) -> None:
        result = run_preorder(
            _make_sets(),
            _make_extraction(),
            _dna_lamp_params(),
            target_na=TargetNucleicAcid.DNA,
        )
        assert result.n_sets_rt_ok == -1

    def test_rna_target_n_sets_rt_ok_non_negative(self) -> None:
        result = run_preorder(
            _make_sets(tm=64.0),
            _make_extraction(),
            _rt_lamp_params(),
            target_na=TargetNucleicAcid.RNA,
        )
        assert result.n_sets_rt_ok >= 0

    def test_target_na_in_result(self) -> None:
        result_dna = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        assert result_dna.target_na == "dna"

        result_rna = run_preorder(
            _make_sets(), _make_extraction(), _rt_lamp_params(), target_na=TargetNucleicAcid.RNA
        )
        assert result_rna.target_na == "rna"


# ---------------------------------------------------------------------------
# top_n slicing
# ---------------------------------------------------------------------------


class TestTopN:
    def test_top_n_limits_sets_assessed(self) -> None:
        sets = _make_sets(n=10, tm=64.0)
        result = run_preorder(
            sets,
            _make_extraction(),
            _rt_lamp_params(),
            target_na=TargetNucleicAcid.RNA,
            top_n=3,
        )
        assert result.n_sets_assessed == 3

    def test_top_n_none_uses_all(self) -> None:
        sets = _make_sets(n=7, tm=64.0)
        result = run_preorder(
            sets,
            _make_extraction(),
            _rt_lamp_params(),
            target_na=TargetNucleicAcid.RNA,
            top_n=None,
        )
        assert result.n_sets_assessed == 7


# ---------------------------------------------------------------------------
# PreorderResult structure
# ---------------------------------------------------------------------------


class TestPreorderResult:
    def test_returns_preorder_result(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        assert isinstance(result, PreorderResult)

    def test_device_window_stored(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params(45.0))
        assert result.device_window_min == pytest.approx(45.0)

    def test_frozen(self) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        with pytest.raises(AttributeError):
            result.status = PreorderStatus.NO_GO  # type: ignore[misc]


# ---------------------------------------------------------------------------
# write_preorder_csv
# ---------------------------------------------------------------------------


class TestWritePreorderCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        out = tmp_path / "preorder.csv"
        write_preorder_csv(result, out)
        assert out.exists()

    def test_csv_has_header(self, tmp_path: Path) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        out = tmp_path / "preorder.csv"
        write_preorder_csv(result, out)
        with out.open(newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert header == ["field", "value"]

    def test_csv_contains_status_row(self, tmp_path: Path) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        out = tmp_path / "preorder.csv"
        write_preorder_csv(result, out)
        with out.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        fields = {row["field"] for row in rows}
        assert "status" in fields

    def test_csv_contains_all_expected_fields(self, tmp_path: Path) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        out = tmp_path / "preorder.csv"
        write_preorder_csv(result, out)
        with out.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        fields = {row["field"] for row in rows}
        for expected in (
            "target_na",
            "n_sets_assessed",
            "lod_95_copies_per_rxn",
            "lod_95_copies_per_ml",
            "ttp_at_lod_95_min",
            "device_window_min",
            "ttp_in_window",
            "status",
        ):
            assert expected in fields, f"Missing field: {expected}"

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        out = tmp_path / "nested" / "dir" / "preorder.csv"
        write_preorder_csv(result, out)
        assert out.exists()

    def test_ttp_in_window_written_as_lowercase_bool(self, tmp_path: Path) -> None:
        result = run_preorder(_make_sets(), _make_extraction(), _dna_lamp_params())
        out = tmp_path / "preorder.csv"
        write_preorder_csv(result, out)
        with out.open(newline="") as fh:
            rows = {row["field"]: row["value"] for row in csv.DictReader(fh)}
        assert rows["ttp_in_window"] in ("true", "false")
