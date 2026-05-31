"""Tests for the surveillance sample-size planner (lamp_forge.sampling).

All tests are pure Python -- no external binaries required.  The
underlying mathematics is deterministic, so exact numeric comparisons
with a small tolerance are safe.

Key identity used for validation:
  Given n samples from an infinite population with effective prevalence
  p_eff = prevalence * sensitivity, the probability of at least one
  positive is:
      P = 1 - (1 - p_eff)^n

  The formula inverts this to find n:
      n = ceil[ log(1-confidence) / log(1 - p_eff) ]
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lamp_forge.sampling import (
    SamplingEstimate,
    SamplingParams,
    _achieved_confidence_infinite,
    _infinite_n,
    sample_size,
    sampling_plan,
    write_sampling_csv,
)

# ---------------------------------------------------------------------------
# SamplingParams validation
# ---------------------------------------------------------------------------


class TestSamplingParamsValidation:
    def test_valid_construction(self) -> None:
        p = SamplingParams(
            target_prevalence=0.05,
            confidence=0.95,
            test_sensitivity=0.95,
        )
        assert p.target_prevalence == 0.05

    def test_defaults(self) -> None:
        p = SamplingParams(target_prevalence=0.10, confidence=0.95)
        assert p.test_sensitivity == 0.95
        assert p.population_size is None

    def test_rejects_prevalence_zero(self) -> None:
        with pytest.raises(ValueError, match="target_prevalence"):
            SamplingParams(target_prevalence=0.0, confidence=0.95)

    def test_rejects_prevalence_one(self) -> None:
        with pytest.raises(ValueError, match="target_prevalence"):
            SamplingParams(target_prevalence=1.0, confidence=0.95)

    def test_rejects_confidence_zero(self) -> None:
        with pytest.raises(ValueError, match="confidence"):
            SamplingParams(target_prevalence=0.05, confidence=0.0)

    def test_rejects_confidence_one(self) -> None:
        with pytest.raises(ValueError, match="confidence"):
            SamplingParams(target_prevalence=0.05, confidence=1.0)

    def test_rejects_sensitivity_zero(self) -> None:
        with pytest.raises(ValueError, match="test_sensitivity"):
            SamplingParams(target_prevalence=0.05, confidence=0.95, test_sensitivity=0.0)

    def test_sensitivity_one_is_valid(self) -> None:
        p = SamplingParams(target_prevalence=0.05, confidence=0.95, test_sensitivity=1.0)
        assert p.test_sensitivity == 1.0

    def test_rejects_negative_population_size(self) -> None:
        with pytest.raises(ValueError, match="population_size"):
            SamplingParams(
                target_prevalence=0.05,
                confidence=0.95,
                population_size=-1,
            )

    def test_rejects_zero_population_size(self) -> None:
        with pytest.raises(ValueError, match="population_size"):
            SamplingParams(
                target_prevalence=0.05,
                confidence=0.95,
                population_size=0,
            )

    def test_positive_population_size_is_valid(self) -> None:
        p = SamplingParams(target_prevalence=0.05, confidence=0.95, population_size=100)
        assert p.population_size == 100


# ---------------------------------------------------------------------------
# _infinite_n (internal helper)
# ---------------------------------------------------------------------------


class TestInfiniteN:
    """n = ceil(log(1-C) / log(1 - p*Se))."""

    def test_known_value_5pct_95pct(self) -> None:
        # p=0.05, Se=1.0, C=0.95: n = ceil(log(0.05)/log(0.95))
        # log(0.05)/log(0.95) = -2.9957 / -0.0513 ~ 58.4 -> 59
        n = _infinite_n(prevalence=0.05, confidence=0.95, sensitivity=1.0)
        assert n == 59

    def test_known_value_10pct_95pct(self) -> None:
        # p=0.10, Se=1.0, C=0.95
        # log(0.05)/log(0.90) ~ 28.4 -> 29
        n = _infinite_n(prevalence=0.10, confidence=0.95, sensitivity=1.0)
        assert n == 29

    def test_higher_prevalence_needs_fewer_samples(self) -> None:
        n_low = _infinite_n(prevalence=0.02, confidence=0.95, sensitivity=0.95)
        n_high = _infinite_n(prevalence=0.20, confidence=0.95, sensitivity=0.95)
        assert n_high < n_low

    def test_higher_confidence_needs_more_samples(self) -> None:
        n_95 = _infinite_n(prevalence=0.05, confidence=0.95, sensitivity=0.95)
        n_99 = _infinite_n(prevalence=0.05, confidence=0.99, sensitivity=0.95)
        assert n_99 > n_95

    def test_lower_sensitivity_needs_more_samples(self) -> None:
        n_high_se = _infinite_n(prevalence=0.10, confidence=0.95, sensitivity=0.99)
        n_low_se = _infinite_n(prevalence=0.10, confidence=0.95, sensitivity=0.70)
        assert n_low_se > n_high_se

    def test_result_is_positive_integer(self) -> None:
        n = _infinite_n(prevalence=0.05, confidence=0.95, sensitivity=0.95)
        assert isinstance(n, int)
        assert n > 0


# ---------------------------------------------------------------------------
# _achieved_confidence_infinite (internal helper)
# ---------------------------------------------------------------------------


class TestAchievedConfidenceInfinite:
    def test_zero_samples_gives_zero(self) -> None:
        assert _achieved_confidence_infinite(0, 0.05, 0.95) == 0.0

    def test_large_n_approaches_one(self) -> None:
        p = _achieved_confidence_infinite(10000, 0.05, 0.95)
        assert p > 0.9999

    def test_roundtrip_with_infinite_n(self) -> None:
        n = _infinite_n(prevalence=0.05, confidence=0.95, sensitivity=0.95)
        achieved = _achieved_confidence_infinite(n, 0.05, 0.95)
        assert achieved >= 0.95


# ---------------------------------------------------------------------------
# sample_size (infinite population)
# ---------------------------------------------------------------------------


class TestSampleSizeInfinite:
    def test_returns_sampling_estimate(self) -> None:
        params = SamplingParams(target_prevalence=0.05, confidence=0.95)
        result = sample_size(params)
        assert isinstance(result, SamplingEstimate)

    def test_n_samples_is_positive_integer(self) -> None:
        params = SamplingParams(target_prevalence=0.10, confidence=0.95)
        result = sample_size(params)
        assert isinstance(result.n_samples, int)
        assert result.n_samples > 0

    def test_achieved_confidence_meets_requirement(self) -> None:
        params = SamplingParams(target_prevalence=0.05, confidence=0.95)
        result = sample_size(params)
        assert result.achieved_confidence >= 0.95

    def test_achieved_not_much_above_requirement(self) -> None:
        # Ceiling rounding should not overshoot by more than one step
        params = SamplingParams(target_prevalence=0.05, confidence=0.95)
        result = sample_size(params)
        assert result.achieved_confidence <= 0.98

    def test_no_finite_population_correction(self) -> None:
        params = SamplingParams(target_prevalence=0.10, confidence=0.95)
        result = sample_size(params)
        assert result.finite_population_corrected is False
        assert result.population_size is None

    def test_high_prevalence_needs_fewer_samples(self) -> None:
        r_low = sample_size(SamplingParams(target_prevalence=0.02, confidence=0.95))
        r_high = sample_size(SamplingParams(target_prevalence=0.20, confidence=0.95))
        assert r_high.n_samples < r_low.n_samples

    def test_perfect_sensitivity_gives_minimum_n(self) -> None:
        r_perfect = sample_size(
            SamplingParams(target_prevalence=0.10, confidence=0.95, test_sensitivity=1.0)
        )
        r_imperfect = sample_size(
            SamplingParams(target_prevalence=0.10, confidence=0.95, test_sensitivity=0.80)
        )
        assert r_perfect.n_samples <= r_imperfect.n_samples

    def test_stored_params_match_input(self) -> None:
        params = SamplingParams(
            target_prevalence=0.07,
            confidence=0.99,
            test_sensitivity=0.90,
        )
        result = sample_size(params)
        assert result.target_prevalence == 0.07
        assert result.confidence == 0.99
        assert result.test_sensitivity == 0.90


# ---------------------------------------------------------------------------
# sample_size (finite population)
# ---------------------------------------------------------------------------


class TestSampleSizeFinite:
    def test_finite_gives_fewer_or_equal_samples(self) -> None:
        params_inf = SamplingParams(target_prevalence=0.10, confidence=0.95)
        params_fin = SamplingParams(target_prevalence=0.10, confidence=0.95, population_size=50)
        r_inf = sample_size(params_inf)
        r_fin = sample_size(params_fin)
        assert r_fin.n_samples <= r_inf.n_samples

    def test_whole_population_when_n_exceeds_N(self) -> None:
        # Very low prevalence + small population -> n_infinite >> N -> test all
        params = SamplingParams(
            target_prevalence=0.001,
            confidence=0.95,
            population_size=20,
        )
        result = sample_size(params)
        assert result.n_samples == 20

    def test_population_size_stored(self) -> None:
        params = SamplingParams(target_prevalence=0.10, confidence=0.95, population_size=200)
        result = sample_size(params)
        assert result.population_size == 200

    def test_fpc_flag_set_when_population_is_small(self) -> None:
        params = SamplingParams(target_prevalence=0.10, confidence=0.95, population_size=50)
        result = sample_size(params)
        # For N=50 and p=0.10, n_infinite ~29; FPC reduces it
        assert result.finite_population_corrected is True

    def test_achieved_confidence_with_finite_population(self) -> None:
        params = SamplingParams(target_prevalence=0.10, confidence=0.95, population_size=100)
        result = sample_size(params)
        assert result.achieved_confidence >= 0.95


# ---------------------------------------------------------------------------
# sampling_plan (table builder)
# ---------------------------------------------------------------------------


class TestSamplingPlan:
    def test_returns_one_estimate_per_prevalence(self) -> None:
        prevs = (0.05, 0.10, 0.20)
        table = sampling_plan(prevs)
        assert len(table) == 3

    def test_prevalences_match_input_order(self) -> None:
        prevs = (0.20, 0.05, 0.10)
        table = sampling_plan(prevs)
        assert [e.target_prevalence for e in table] == list(prevs)

    def test_all_confidence_values_meet_requirement(self) -> None:
        table = sampling_plan((0.05, 0.10, 0.20), confidence=0.95)
        for e in table:
            assert e.achieved_confidence >= 0.95

    def test_consistent_sensitivity_across_table(self) -> None:
        table = sampling_plan((0.05, 0.10), confidence=0.95, sensitivity=0.80)
        for e in table:
            assert e.test_sensitivity == 0.80

    def test_monotone_n_decreasing_with_prevalence(self) -> None:
        prevs = (0.01, 0.05, 0.10, 0.20, 0.50)
        table = sampling_plan(prevs)
        ns = [e.n_samples for e in table]
        assert ns == sorted(ns, reverse=True)

    def test_finite_population_forwarded(self) -> None:
        table = sampling_plan((0.05, 0.10), population_size=50)
        for e in table:
            assert e.population_size == 50


# ---------------------------------------------------------------------------
# write_sampling_csv
# ---------------------------------------------------------------------------


class TestWriteSamplingCsv:
    def test_creates_file_and_header(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.05, 0.10, 0.20))
        out = tmp_path / "plan.csv"
        write_sampling_csv(estimates, out)
        assert out.exists()
        rows = list(csv.reader(out.open()))
        assert rows[0][0] == "target_prevalence_pct"
        assert len(rows) == 1 + 3  # header + 3 data rows

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.10,))
        out = tmp_path / "sub" / "dir" / "plan.csv"
        write_sampling_csv(estimates, out)
        assert out.exists()

    def test_n_samples_column_is_integer_string(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.05,))
        out = tmp_path / "plan.csv"
        write_sampling_csv(estimates, out)
        rows = list(csv.DictReader(out.open()))
        assert int(rows[0]["n_samples"]) > 0

    def test_prevalence_pct_matches_input(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.05, 0.10))
        out = tmp_path / "plan.csv"
        write_sampling_csv(estimates, out)
        rows = list(csv.DictReader(out.open()))
        assert abs(float(rows[0]["target_prevalence_pct"]) - 5.0) < 0.01
        assert abs(float(rows[1]["target_prevalence_pct"]) - 10.0) < 0.01

    def test_achieved_confidence_numeric(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.10,))
        out = tmp_path / "plan.csv"
        write_sampling_csv(estimates, out)
        rows = list(csv.DictReader(out.open()))
        val = float(rows[0]["achieved_confidence_pct"])
        assert val >= 95.0

    def test_population_size_empty_for_infinite(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.10,))
        out = tmp_path / "plan.csv"
        write_sampling_csv(estimates, out)
        rows = list(csv.DictReader(out.open()))
        assert rows[0]["population_size"] == ""

    def test_population_size_present_for_finite(self, tmp_path: Path) -> None:
        estimates = sampling_plan((0.10,), population_size=100)
        out = tmp_path / "plan.csv"
        write_sampling_csv(estimates, out)
        rows = list(csv.DictReader(out.open()))
        assert rows[0]["population_size"] == "100"


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------


class TestSamplingPlanCli:
    def _run(self, args: list[str]) -> object:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        return runner.invoke(cli, ["sampling-plan", *args], obj={})

    def test_basic_invocation(self) -> None:
        result = self._run(["--prevalence", "0.10"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "n_samples" in result.output  # type: ignore[union-attr]

    def test_multiple_prevalences(self) -> None:
        result = self._run(["--prevalence", "0.05", "--prevalence", "0.10", "--prevalence", "0.20"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        # Three data rows should appear (each row starts with a decimal prevalence)
        import re

        lines = [
            ln
            for ln in result.output.splitlines()  # type: ignore[union-attr]
            if ln.strip() and re.match(r"\s*\d+\.\d+\s+\d+", ln)
        ]
        assert len(lines) >= 3

    def test_output_contains_achieved_confidence(self) -> None:
        result = self._run(["--prevalence", "0.10"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "95" in result.output  # type: ignore[union-attr]

    def test_writes_csv_when_requested(self, tmp_path: Path) -> None:
        out = tmp_path / "sampling_plan.csv"
        result = self._run(["--prevalence", "0.10", "--out-csv", str(out)])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_rejects_invalid_prevalence(self) -> None:
        result = self._run(["--prevalence", "1.5"])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_rejects_prevalence_zero(self) -> None:
        result = self._run(["--prevalence", "0.0"])
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_custom_confidence(self) -> None:
        result = self._run(["--prevalence", "0.10", "--confidence", "0.99"])
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "99" in result.output  # type: ignore[union-attr]

    def test_finite_population_shows_fpc_column(self) -> None:
        result = self._run(
            [
                "--prevalence",
                "0.10",
                "--population-size",
                "50",
            ]
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "FPC" in result.output  # type: ignore[union-attr]


# ---------------------------------------------------------------------------
# BioVind vertical spot-checks
# ---------------------------------------------------------------------------


class TestBioVindVerticals:
    """Spot-checks against known textbook values for the three BioVind verticals.

    Reference values computed with the standard formula:
        n = ceil(log(1-C) / log(1 - p*Se))
    """

    def test_oilfield_5pct_srb_detection(self) -> None:
        """Oil & gas: detect SRB presence if >= 5% of wells are positive (Se=0.95)."""
        params = SamplingParams(
            target_prevalence=0.05,
            confidence=0.95,
            test_sensitivity=0.95,
        )
        result = sample_size(params)
        # n = ceil(log(0.05)/log(1-0.0475)) ~ ceil(62.2) = 63
        assert 55 <= result.n_samples <= 70

    def test_farm_asfv_10pct_prevalence(self) -> None:
        """Farm biosecurity: detect ASFV if >= 10% of animals in herd are positive."""
        params = SamplingParams(
            target_prevalence=0.10,
            confidence=0.95,
            test_sensitivity=0.95,
        )
        result = sample_size(params)
        # n = ceil(log(0.05)/log(1-0.095)) ~ ceil(30.6) = 31
        assert 28 <= result.n_samples <= 35

    def test_poc_mtb_1pct_prevalence(self) -> None:
        """Human POC: detect MTB if >= 1% of patients are positive (Se=0.95)."""
        params = SamplingParams(
            target_prevalence=0.01,
            confidence=0.95,
            test_sensitivity=0.95,
        )
        result = sample_size(params)
        # n = ceil(log(0.05)/log(1-0.0095)) ~ ceil(313.6) = 314
        assert 300 <= result.n_samples <= 330

    def test_farm_finite_herd_gives_fewer_samples_than_infinite(self) -> None:
        """Finite herd (50 pigs) needs fewer tests than infinite-population formula."""
        r_inf = sample_size(SamplingParams(target_prevalence=0.10, confidence=0.95))
        r_fin = sample_size(
            SamplingParams(target_prevalence=0.10, confidence=0.95, population_size=50)
        )
        assert r_fin.n_samples <= r_inf.n_samples
