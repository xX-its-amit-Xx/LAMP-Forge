"""Tests for the diagnostic accuracy metrics module.

All tests are pure Python (no external binaries, no network).  The statistical
formulae are deterministic so exact floating-point comparisons with tolerance
are safe.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from lamp_forge.accuracy import (
    AccuracyMetrics,
    ContingencyTable,
    accuracy_table,
    compute_accuracy,
    wilson_ci,
    write_accuracy_csv,
)

# ---------------------------------------------------------------------------
# ContingencyTable validation
# ---------------------------------------------------------------------------


class TestContingencyTable:
    def test_valid_table(self) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        assert t.n_positive_ref == 50
        assert t.n_negative_ref == 50
        assert t.n_total == 100

    def test_perfect_assay(self) -> None:
        t = ContingencyTable(tp=50, fp=0, tn=50, fn=0)
        assert t.n_positive_ref == 50
        assert t.n_negative_ref == 50

    def test_rejects_negative_count(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            ContingencyTable(tp=-1, fp=0, tn=50, fn=0)

    def test_rejects_zero_positive_ref(self) -> None:
        with pytest.raises(ValueError, match="positive reference"):
            ContingencyTable(tp=0, fp=5, tn=50, fn=0)

    def test_rejects_zero_negative_ref(self) -> None:
        with pytest.raises(ValueError, match="negative reference"):
            ContingencyTable(tp=50, fp=0, tn=0, fn=0)

    def test_minimum_valid_table(self) -> None:
        t = ContingencyTable(tp=1, fp=0, tn=1, fn=0)
        assert t.n_total == 2


# ---------------------------------------------------------------------------
# wilson_ci
# ---------------------------------------------------------------------------


class TestWilsonCI:
    def test_50_of_100(self) -> None:
        ci = wilson_ci(50, 100)
        assert abs(ci.point_estimate - 0.50) < 1e-12
        assert 0.40 < ci.lower < 0.50
        assert 0.50 < ci.upper < 0.60

    def test_perfect_sensitivity(self) -> None:
        ci = wilson_ci(50, 50)
        assert ci.point_estimate == 1.0
        assert ci.lower > 0.0  # Wilson CI does not collapse to [1, 1] for k==n
        assert ci.upper == 1.0

    def test_zero_proportion(self) -> None:
        ci = wilson_ci(0, 50)
        assert ci.point_estimate == 0.0
        assert ci.lower == 0.0
        assert ci.upper > 0.0

    def test_95_ci_known_values(self) -> None:
        # For 47/50, 95% Wilson CI approximately [0.827, 0.988]
        ci = wilson_ci(47, 50, confidence=0.95)
        assert abs(ci.point_estimate - 47 / 50) < 1e-12
        assert 0.80 < ci.lower < 0.87
        assert 0.97 < ci.upper <= 1.0

    def test_confidence_stored(self) -> None:
        ci = wilson_ci(20, 30, confidence=0.99)
        assert ci.confidence == 0.99

    def test_rejects_n_zero(self) -> None:
        with pytest.raises(ValueError, match="n must be"):
            wilson_ci(0, 0)

    def test_lower_upper_ordering(self) -> None:
        for k, n in [(1, 10), (5, 10), (9, 10), (0, 20), (20, 20)]:
            ci = wilson_ci(k, n)
            assert ci.lower <= ci.point_estimate
            assert ci.point_estimate <= ci.upper

    def test_bounds_in_unit_interval(self) -> None:
        for k, n in [(0, 5), (5, 5), (1, 2), (3, 4)]:
            ci = wilson_ci(k, n)
            assert 0.0 <= ci.lower <= 1.0
            assert 0.0 <= ci.upper <= 1.0


# ---------------------------------------------------------------------------
# compute_accuracy
# ---------------------------------------------------------------------------


class TestComputeAccuracy:
    def _typical_table(self) -> ContingencyTable:
        return ContingencyTable(tp=47, fp=2, tn=48, fn=3)

    def test_returns_accuracy_metrics(self) -> None:
        m = compute_accuracy(self._typical_table())
        assert isinstance(m, AccuracyMetrics)

    def test_sensitivity_point_estimate(self) -> None:
        m = compute_accuracy(self._typical_table())
        # Sn = 47 / (47 + 3) = 0.94
        assert abs(m.sensitivity.point_estimate - 47 / 50) < 1e-12

    def test_specificity_point_estimate(self) -> None:
        m = compute_accuracy(self._typical_table())
        # Sp = 48 / (48 + 2) = 0.96
        assert abs(m.specificity.point_estimate - 48 / 50) < 1e-12

    def test_lr_positive(self) -> None:
        m = compute_accuracy(self._typical_table())
        sn = 47 / 50
        sp = 48 / 50
        expected_lr_pos = sn / (1 - sp)
        assert abs(m.lr_positive - expected_lr_pos) < 1e-9

    def test_lr_negative(self) -> None:
        m = compute_accuracy(self._typical_table())
        sn = 47 / 50
        sp = 48 / 50
        expected_lr_neg = (1 - sn) / sp
        assert abs(m.lr_negative - expected_lr_neg) < 1e-9

    def test_youden_j(self) -> None:
        m = compute_accuracy(self._typical_table())
        sn = 47 / 50
        sp = 48 / 50
        assert abs(m.youden_j - (sn + sp - 1)) < 1e-12

    def test_ppv_at_50pct_prevalence(self) -> None:
        # At prevalence = 0.5 and balanced TP/FP:
        # PPV = (Sn * 0.5) / (Sn * 0.5 + (1-Sp) * 0.5) = Sn / (Sn + 1 - Sp)
        m = compute_accuracy(self._typical_table(), prevalence=0.5)
        sn = 47 / 50
        sp = 48 / 50
        expected = (sn * 0.5) / (sn * 0.5 + (1 - sp) * 0.5)
        assert abs(m.ppv - expected) < 1e-9

    def test_npv_at_10pct_prevalence(self) -> None:
        m = compute_accuracy(self._typical_table(), prevalence=0.10)
        sn = 47 / 50
        sp = 48 / 50
        prev = 0.10
        expected = (sp * (1 - prev)) / (sp * (1 - prev) + (1 - sn) * prev)
        assert abs(m.npv - expected) < 1e-9

    def test_prevalence_stored(self) -> None:
        m = compute_accuracy(self._typical_table(), prevalence=0.07)
        assert m.prevalence == 0.07

    def test_perfect_assay_lr_positive_infinity(self) -> None:
        t = ContingencyTable(tp=50, fp=0, tn=50, fn=0)
        m = compute_accuracy(t)
        assert math.isinf(m.lr_positive)

    def test_perfect_assay_youden_j_one(self) -> None:
        t = ContingencyTable(tp=50, fp=0, tn=50, fn=0)
        m = compute_accuracy(t)
        assert abs(m.youden_j - 1.0) < 1e-12

    def test_rejects_invalid_prevalence(self) -> None:
        t = self._typical_table()
        with pytest.raises(ValueError, match="prevalence"):
            compute_accuracy(t, prevalence=0.0)
        with pytest.raises(ValueError, match="prevalence"):
            compute_accuracy(t, prevalence=1.0)

    def test_ppv_increases_with_prevalence(self) -> None:
        t = self._typical_table()
        ppvs = [compute_accuracy(t, prevalence=p).ppv for p in (0.01, 0.05, 0.10, 0.50)]
        assert ppvs == sorted(ppvs)

    def test_npv_decreases_with_prevalence(self) -> None:
        t = self._typical_table()
        npvs = [compute_accuracy(t, prevalence=p).npv for p in (0.01, 0.05, 0.10, 0.50)]
        assert npvs == sorted(npvs, reverse=True)

    def test_min_panel_one_specimen_each(self) -> None:
        t = ContingencyTable(tp=1, fp=0, tn=1, fn=0)
        m = compute_accuracy(t, prevalence=0.10)
        assert m.sensitivity.point_estimate == 1.0
        assert m.specificity.point_estimate == 1.0


# ---------------------------------------------------------------------------
# accuracy_table
# ---------------------------------------------------------------------------


class TestAccuracyTable:
    def test_returns_one_row_per_prevalence(self) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=(0.05, 0.10, 0.20))
        assert len(rows) == 3

    def test_sensitivity_same_across_rows(self) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=(0.01, 0.05, 0.10))
        sn_vals = [r.sensitivity.point_estimate for r in rows]
        assert len(set(sn_vals)) == 1  # all equal

    def test_ppv_varies_across_rows(self) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=(0.01, 0.05, 0.20))
        ppvs = [r.ppv for r in rows]
        assert ppvs[0] < ppvs[1] < ppvs[2]

    def test_empty_prevalences_returns_empty(self) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=())
        assert rows == []


# ---------------------------------------------------------------------------
# write_accuracy_csv
# ---------------------------------------------------------------------------


class TestWriteAccuracyCsv:
    def test_creates_csv_file(self, tmp_path: Path) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=(0.05, 0.10))
        out = tmp_path / "accuracy.csv"
        write_accuracy_csv(rows, out)
        assert out.exists()

    def test_csv_has_header_and_two_rows(self, tmp_path: Path) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=(0.05, 0.10))
        out = tmp_path / "accuracy.csv"
        write_accuracy_csv(rows, out)
        with out.open(encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            data_rows = list(reader)
        assert len(data_rows) == 2

    def test_csv_columns_present(self, tmp_path: Path) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t)
        out = tmp_path / "accuracy.csv"
        write_accuracy_csv(rows, out)
        with out.open(encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            fieldnames = reader.fieldnames or []
        expected = {
            "prevalence",
            "tp",
            "fp",
            "tn",
            "fn",
            "sensitivity",
            "specificity",
            "ppv",
            "npv",
            "youden_j",
            "lr_positive",
            "lr_negative",
        }
        assert expected.issubset(set(fieldnames))

    def test_csv_sensitivity_value_correct(self, tmp_path: Path) -> None:
        t = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
        rows = accuracy_table(t, prevalences=(0.10,))
        out = tmp_path / "accuracy.csv"
        write_accuracy_csv(rows, out)
        with out.open(encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            row = next(reader)
        assert abs(float(row["sensitivity"]) - 47 / 50) < 1e-4

    def test_csv_inf_lr_positive_string(self, tmp_path: Path) -> None:
        t = ContingencyTable(tp=50, fp=0, tn=50, fn=0)
        rows = accuracy_table(t, prevalences=(0.10,))
        out = tmp_path / "accuracy.csv"
        write_accuracy_csv(rows, out)
        with out.open(encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            row = next(reader)
        assert row["lr_positive"] == "inf"


# ---------------------------------------------------------------------------
# CLI integration (click CliRunner)
# ---------------------------------------------------------------------------


class TestAccuracyCli:
    def test_basic_invocation(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["accuracy", "--tp", "47", "--fp", "2", "--tn", "48", "--fn", "3"],
        )
        assert result.exit_code == 0, result.output
        assert "Sensitivity" in result.output
        assert "Specificity" in result.output
        assert "LR+" in result.output
        assert "Youden" in result.output

    def test_biovind_criteria_met_message(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        # Perfect assay -- all criteria met
        result = runner.invoke(
            cli,
            ["accuracy", "--tp", "50", "--fp", "0", "--tn", "50", "--fn", "0"],
        )
        assert result.exit_code == 0, result.output
        assert "BioVind acceptance criteria met" in result.output

    def test_biovind_criteria_not_met_message(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        # Poor assay: low sensitivity and specificity
        result = runner.invoke(
            cli,
            ["accuracy", "--tp", "10", "--fp", "10", "--tn", "10", "--fn", "10"],
        )
        assert result.exit_code == 0, result.output
        assert "NOT fully met" in result.output

    def test_csv_output_created(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        out_csv = tmp_path / "acc.csv"
        result = runner.invoke(
            cli,
            [
                "accuracy",
                "--tp",
                "47",
                "--fp",
                "2",
                "--tn",
                "48",
                "--fn",
                "3",
                "--out-csv",
                str(out_csv),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out_csv.exists()

    def test_invalid_tp_negative(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["accuracy", "--tp", "-1", "--fp", "0", "--tn", "50", "--fn", "0"],
        )
        assert result.exit_code != 0

    def test_custom_prevalence(self) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "accuracy",
                "--tp",
                "47",
                "--fp",
                "2",
                "--tn",
                "48",
                "--fn",
                "3",
                "--prevalence",
                "0.03",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "3.00" in result.output
