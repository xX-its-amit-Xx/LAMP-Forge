"""Tests for validation_plan module."""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

from lamp_forge.validation_plan import (
    ValidationPlan,
    ValidationSizeRow,
    build_validation_plan,
    n_expected_performance,
    n_zero_failure,
    wilson_lower,
    write_plan_csv,
)

# ---------------------------------------------------------------------------
# wilson_lower
# ---------------------------------------------------------------------------


def test_wilson_lower_all_successes_approx_p_min() -> None:
    """wilson_lower(n, n) should be close to n/(n+z^2)."""
    z = 1.959963985  # 95% two-sided
    n = 100
    expected = n / (n + z**2)
    got = wilson_lower(n, n, confidence=0.95)
    assert abs(got - expected) < 1e-6


def test_wilson_lower_zero_successes() -> None:
    """With k=0, lower bound is 0."""
    assert wilson_lower(0, 50, 0.95) == 0.0


def test_wilson_lower_monotonic_in_k() -> None:
    """Wilson lower bound should be non-decreasing as k increases."""
    n = 30
    lowers = [wilson_lower(k, n, 0.95) for k in range(n + 1)]
    for a, b in itertools.pairwise(lowers):
        assert b >= a - 1e-12


def test_wilson_lower_invalid_n() -> None:
    with pytest.raises(ValueError, match="n must be"):
        wilson_lower(0, 0, 0.95)


def test_wilson_lower_k_out_of_range() -> None:
    with pytest.raises(ValueError, match="k must be in"):
        wilson_lower(5, 3, 0.95)


# ---------------------------------------------------------------------------
# n_zero_failure
# ---------------------------------------------------------------------------


def test_n_zero_failure_sn_95() -> None:
    """For Sn>=0.95 at 95% CI, zero-failure n should be ~73."""
    n = n_zero_failure(0.95, 0.95)
    assert n >= 73
    # The lower CI bound at n must be >= 0.95
    assert wilson_lower(n, n, 0.95) >= 0.95 - 1e-10


def test_n_zero_failure_sp_99() -> None:
    """For Sp>=0.99 at 95% CI, zero-failure n should be ~381."""
    n = n_zero_failure(0.99, 0.95)
    assert n >= 380
    assert wilson_lower(n, n, 0.95) >= 0.99 - 1e-10


def test_n_zero_failure_lower_bound_tight() -> None:
    """n-1 must have lower CI below p_min."""
    p_min = 0.95
    n = n_zero_failure(p_min, 0.95)
    if n > 1:
        assert wilson_lower(n - 1, n - 1, 0.95) < p_min + 1e-10


def test_n_zero_failure_invalid_p_min() -> None:
    with pytest.raises(ValueError):
        n_zero_failure(0.0, 0.95)

    with pytest.raises(ValueError):
        n_zero_failure(1.0, 0.95)


# ---------------------------------------------------------------------------
# n_expected_performance
# ---------------------------------------------------------------------------


def test_n_expected_performance_returns_tuple() -> None:
    n, power = n_expected_performance(0.95, 0.95, 0.97, 0.80)
    assert isinstance(n, int)
    assert isinstance(power, float)
    assert n >= 1
    assert power >= 0.80 - 1e-6


def test_n_expected_performance_achieves_power() -> None:
    """Achieved power at returned n must be >= requested power."""
    nep, achieved = n_expected_performance(0.95, 0.95, 0.97, 0.80)
    assert nep >= 1
    assert achieved >= 0.80 - 1e-6


def test_n_expected_performance_high_power_needs_more() -> None:
    """Higher required power demands more specimens."""
    n80, _ = n_expected_performance(0.95, 0.95, 0.97, 0.80)
    n90, _ = n_expected_performance(0.95, 0.95, 0.97, 0.90)
    assert n90 >= n80


def test_n_expected_performance_invalid_args() -> None:
    with pytest.raises(ValueError):
        n_expected_performance(0.0, 0.95, 0.97, 0.80)

    with pytest.raises(ValueError):
        n_expected_performance(0.95, 0.95, 0.0, 0.80)

    with pytest.raises(ValueError):
        n_expected_performance(0.95, 0.95, 0.97, 0.0)


# ---------------------------------------------------------------------------
# build_validation_plan
# ---------------------------------------------------------------------------


def test_build_validation_plan_default() -> None:
    plan = build_validation_plan()
    assert isinstance(plan, ValidationPlan)
    assert len(plan.rows) == 2
    sn_row = next(r for r in plan.rows if r.metric == "sensitivity")
    sp_row = next(r for r in plan.rows if r.metric == "specificity")

    assert sn_row.p_min == pytest.approx(0.95)
    assert sp_row.p_min == pytest.approx(0.99)

    # Zero-failure lower bounds must meet specifications
    assert sn_row.wilson_lower_zf >= 0.95 - 1e-10
    assert sp_row.wilson_lower_zf >= 0.99 - 1e-10

    # n_positive / n_negative match the zero-failure values
    assert plan.n_positive_required == sn_row.n_zero_failure
    assert plan.n_negative_required == sp_row.n_zero_failure


def test_build_validation_plan_custom_spec() -> None:
    plan = build_validation_plan(sn_min=0.90, sp_min=0.95, confidence=0.90)
    sn_row = next(r for r in plan.rows if r.metric == "sensitivity")
    sp_row = next(r for r in plan.rows if r.metric == "specificity")
    assert wilson_lower(sn_row.n_zero_failure, sn_row.n_zero_failure, 0.90) >= 0.90 - 1e-10
    assert wilson_lower(sp_row.n_zero_failure, sp_row.n_zero_failure, 0.90) >= 0.95 - 1e-10


def test_build_validation_plan_rows_are_dataclasses() -> None:
    plan = build_validation_plan()
    for row in plan.rows:
        assert isinstance(row, ValidationSizeRow)
        assert row.achieved_power >= row.power - 1e-6


# ---------------------------------------------------------------------------
# write_plan_csv
# ---------------------------------------------------------------------------


def test_write_plan_csv(tmp_path: Path) -> None:
    plan = build_validation_plan()
    out = tmp_path / "plan.csv"
    write_plan_csv(plan, out)
    text = out.read_text()
    assert "sensitivity" in text
    assert "specificity" in text
    assert "n_zero_failure" in text
    assert "n_expected" in text

    import csv

    with out.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    for row in rows:
        assert int(row["n_zero_failure"]) > 0
        assert float(row["wilson_lower_zf"]) >= float(row["p_min"]) - 1e-6
