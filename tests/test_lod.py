"""Tests for the LOD / copies-per-reaction estimator.

All tests are pure Python (no external binaries). The Poisson math is
deterministic, so exact floating-point comparisons with a tolerance are safe.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lamp_forge.lod import (
    ExtractionParams,
    LodEstimate,
    detection_probability_at,
    estimate_lod,
    lod_table,
    poisson_lod,
    write_lod_csv,
)

# ---------------------------------------------------------------------------
# poisson_lod
# ---------------------------------------------------------------------------


class TestPoissonLod:
    def test_lod95_approx_3(self) -> None:
        # -ln(0.05) = 2.99573…
        result = poisson_lod(0.95)
        assert abs(result - 2.9957) < 1e-3

    def test_lod99_approx_4_6(self) -> None:
        # -ln(0.01) = 4.60517…
        result = poisson_lod(0.99)
        assert abs(result - 4.6052) < 1e-3

    def test_lod999_approx_6_9(self) -> None:
        result = poisson_lod(0.999)
        assert abs(result - 6.9078) < 1e-3

    def test_roundtrip_with_detection_probability_at(self) -> None:
        for p in (0.50, 0.90, 0.95, 0.99):
            lam = poisson_lod(p)
            recovered = detection_probability_at(lam)
            assert abs(recovered - p) < 1e-12

    def test_raises_below_zero(self) -> None:
        with pytest.raises(ValueError, match="detection_probability"):
            poisson_lod(0.0)

    def test_raises_at_one(self) -> None:
        with pytest.raises(ValueError, match="detection_probability"):
            poisson_lod(1.0)

    def test_raises_above_one(self) -> None:
        with pytest.raises(ValueError, match="detection_probability"):
            poisson_lod(1.5)

    def test_monotonically_increasing(self) -> None:
        probs = [0.50, 0.80, 0.90, 0.95, 0.99, 0.999]
        lods = [poisson_lod(p) for p in probs]
        assert lods == sorted(lods)


# ---------------------------------------------------------------------------
# detection_probability_at
# ---------------------------------------------------------------------------


class TestDetectionProbabilityAt:
    def test_zero_copies_gives_zero_probability(self) -> None:
        assert detection_probability_at(0.0) == 0.0

    def test_large_copy_number_approaches_one(self) -> None:
        assert detection_probability_at(100.0) > 0.9999

    def test_raises_on_negative(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            detection_probability_at(-1.0)


# ---------------------------------------------------------------------------
# ExtractionParams validation
# ---------------------------------------------------------------------------


class TestExtractionParams:
    def test_valid_construction(self) -> None:
        p = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        assert p.sample_volume_ul == 200.0

    def test_conversion_factor(self) -> None:
        # 200 uL sample * 0.5 * (5/50) = 10 uL effective sample in reaction
        p = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        assert abs(p.copies_per_rxn_per_copy_per_ul - 10.0) < 1e-10

    def test_rejects_zero_sample_volume(self) -> None:
        with pytest.raises(ValueError, match="sample_volume_ul"):
            ExtractionParams(
                sample_volume_ul=0.0,
                extraction_efficiency=0.5,
                eluate_volume_ul=50.0,
                reaction_input_ul=5.0,
            )

    def test_rejects_negative_sample_volume(self) -> None:
        with pytest.raises(ValueError, match="sample_volume_ul"):
            ExtractionParams(
                sample_volume_ul=-10.0,
                extraction_efficiency=0.5,
                eluate_volume_ul=50.0,
                reaction_input_ul=5.0,
            )

    def test_rejects_efficiency_above_one(self) -> None:
        with pytest.raises(ValueError, match="extraction_efficiency"):
            ExtractionParams(
                sample_volume_ul=200.0,
                extraction_efficiency=1.1,
                eluate_volume_ul=50.0,
                reaction_input_ul=5.0,
            )

    def test_rejects_efficiency_of_zero(self) -> None:
        with pytest.raises(ValueError, match="extraction_efficiency"):
            ExtractionParams(
                sample_volume_ul=200.0,
                extraction_efficiency=0.0,
                eluate_volume_ul=50.0,
                reaction_input_ul=5.0,
            )

    def test_efficiency_of_one_is_valid(self) -> None:
        p = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=1.0,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        assert p.extraction_efficiency == 1.0

    def test_rejects_reaction_input_exceeding_eluate(self) -> None:
        with pytest.raises(ValueError, match="reaction_input_ul"):
            ExtractionParams(
                sample_volume_ul=200.0,
                extraction_efficiency=0.5,
                eluate_volume_ul=50.0,
                reaction_input_ul=60.0,
            )

    def test_reaction_input_equal_to_eluate_is_valid(self) -> None:
        p = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=50.0,
        )
        assert p.reaction_input_ul == 50.0


# ---------------------------------------------------------------------------
# estimate_lod
# ---------------------------------------------------------------------------


class TestEstimateLod:
    def _params(self) -> ExtractionParams:
        # 200 uL sample, 50% efficiency, 50 uL eluate, 5 uL to reaction
        # conversion factor = 200 * 0.5 * 5/50 = 10 uL
        return ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )

    def test_returns_lod_estimate(self) -> None:
        result = estimate_lod(self._params())
        assert isinstance(result, LodEstimate)

    def test_lod95_copies_per_reaction(self) -> None:
        result = estimate_lod(self._params(), detection_probability=0.95)
        # lambda_LOD = -ln(0.05) ~ 2.9957
        assert abs(result.lod_copies_per_reaction - 2.9957) < 1e-3

    def test_lod95_copies_per_ml(self) -> None:
        # conversion = 10 uL; lod_per_ul = 2.9957/10 = 0.29957; per mL = 299.57
        result = estimate_lod(self._params(), detection_probability=0.95)
        assert abs(result.lod_copies_per_ml - 299.57) < 1.0

    def test_lod_per_ul_and_per_ml_consistent(self) -> None:
        result = estimate_lod(self._params())
        assert abs(result.lod_copies_per_ml - result.lod_copies_per_ul * 1000.0) < 1e-9

    def test_genome_eq_per_ml_alias(self) -> None:
        result = estimate_lod(self._params())
        assert result.lod_genome_eq_per_ml == result.lod_copies_per_ml

    def test_detection_probability_stored(self) -> None:
        result = estimate_lod(self._params(), detection_probability=0.99)
        assert result.detection_probability == 0.99

    def test_higher_probability_gives_higher_lod(self) -> None:
        params = self._params()
        lod95 = estimate_lod(params, 0.95)
        lod99 = estimate_lod(params, 0.99)
        assert lod99.lod_copies_per_ml > lod95.lod_copies_per_ml

    def test_larger_effective_volume_gives_lower_lod(self) -> None:
        # More sample in the reaction => lower LOD in original sample.
        small_sample = ExtractionParams(
            sample_volume_ul=50.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        large_sample = ExtractionParams(
            sample_volume_ul=500.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        lod_small = estimate_lod(small_sample, 0.95)
        lod_large = estimate_lod(large_sample, 0.95)
        assert lod_large.lod_copies_per_ml < lod_small.lod_copies_per_ml

    def test_higher_extraction_efficiency_gives_lower_lod(self) -> None:
        low_eff = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.3,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        high_eff = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.8,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        lod_low = estimate_lod(low_eff, 0.95)
        lod_high = estimate_lod(high_eff, 0.95)
        assert lod_high.lod_copies_per_ml < lod_low.lod_copies_per_ml

    def test_raises_on_bad_probability(self) -> None:
        with pytest.raises(ValueError):
            estimate_lod(self._params(), detection_probability=1.0)


# ---------------------------------------------------------------------------
# lod_table
# ---------------------------------------------------------------------------


class TestLodTable:
    def test_returns_one_estimate_per_probability(self) -> None:
        params = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        probs = (0.90, 0.95, 0.99)
        table = lod_table(params, probs)
        assert len(table) == 3

    def test_probabilities_match_input_order(self) -> None:
        params = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        probs = (0.99, 0.90, 0.95)
        table = lod_table(params, probs)
        assert [e.detection_probability for e in table] == list(probs)

    def test_default_probabilities(self) -> None:
        params = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        table = lod_table(params)
        assert len(table) == 4  # default: 0.90, 0.95, 0.99, 0.999


# ---------------------------------------------------------------------------
# write_lod_csv
# ---------------------------------------------------------------------------


class TestWriteLodCsv:
    def test_creates_file_and_writes_header(self, tmp_path: Path) -> None:
        params = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        estimates = lod_table(params)
        out = tmp_path / "lod.csv"
        write_lod_csv(estimates, out)
        assert out.exists()
        rows = list(csv.reader(out.open()))
        assert rows[0][0] == "detection_probability"
        assert len(rows) == 1 + 4  # header + 4 default probabilities

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        params = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        out = tmp_path / "sub" / "dir" / "lod.csv"
        write_lod_csv(lod_table(params), out)
        assert out.exists()

    def test_lod_copies_per_ml_column_is_numeric(self, tmp_path: Path) -> None:
        params = ExtractionParams(
            sample_volume_ul=200.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=50.0,
            reaction_input_ul=5.0,
        )
        out = tmp_path / "lod.csv"
        write_lod_csv(lod_table(params), out)
        rows = list(csv.DictReader(out.open()))
        for row in rows:
            val = float(row["lod_copies_per_ml_sample"])
            assert val > 0.0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestLodCli:
    def test_basic_invocation(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "lod",
                "--sample-volume", "200",
                "--eluate-volume", "50",
                "--reaction-input", "5",
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert "copies/mL" in result.output or "LOD" in result.output or "lambda" in result.output

    def test_output_contains_table_rows(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "lod",
                "--sample-volume", "200",
                "--eluate-volume", "50",
                "--reaction-input", "5",
                "--probability", "0.95",
                "--probability", "0.99",
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        # Two data rows (one per --probability) should appear.
        lines = [ln for ln in result.output.splitlines() if ln.strip()]
        # At least two lines should contain the probabilities.
        numeric_lines = [ln for ln in lines if "0.950" in ln or "0.990" in ln]
        assert len(numeric_lines) >= 2

    def test_writes_csv_when_requested(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        out = tmp_path / "lod_out.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "lod",
                "--sample-volume", "200",
                "--eluate-volume", "50",
                "--reaction-input", "5",
                "--out-csv", str(out),
            ],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_rejects_bad_efficiency(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "lod",
                "--sample-volume", "200",
                "--eluate-volume", "50",
                "--reaction-input", "5",
                "--efficiency", "1.5",
            ],
            obj={},
        )
        assert result.exit_code != 0

    def test_rejects_reaction_input_exceeding_eluate(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "lod",
                "--sample-volume", "200",
                "--eluate-volume", "50",
                "--reaction-input", "100",
            ],
            obj={},
        )
        assert result.exit_code != 0
