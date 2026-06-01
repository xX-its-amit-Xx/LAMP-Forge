"""CLI smoke tests for BioVind risk-assessment commands.

Covers the four risk commands that have module-level tests but no
CLI-level tests: brucella-risk, csfv-risk, swine-enteric-risk,
poultry-risk.

All tests use Click's CliRunner -- no external binaries required.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from lamp_forge.cli import cli

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _run(*args: str) -> object:
    runner = CliRunner()
    return runner.invoke(cli, list(args), obj={})


# ===========================================================================
# brucella-risk
# ===========================================================================


class TestBrucellaRiskCli:
    """CLI tests for ``lamp-forge brucella-risk``."""

    def _invoke(self, *args: str) -> object:
        return _run("brucella-risk", *args)

    def test_negative_exits_0(self) -> None:
        result = self._invoke("--no-is711")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_negative_shows_negative(self) -> None:
        result = self._invoke("--no-is711")
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_positive_cattle_shows_high(self) -> None:
        result = self._invoke("--is711", "--context", "cattle")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_positive_human_shows_critical(self) -> None:
        result = self._invoke("--is711", "--context", "human")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_positive_environmental_shows_low(self) -> None:
        result = self._invoke("--is711", "--context", "environmental")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "LOW" in result.output  # type: ignore[union-attr]

    def test_positive_small_ruminant_shows_high(self) -> None:
        result = self._invoke("--is711", "--context", "small_ruminant")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_output_contains_is711_label(self) -> None:
        result = self._invoke("--is711")
        assert "IS711" in result.output  # type: ignore[union-attr]

    def test_human_positive_shows_immediate_report(self) -> None:
        result = self._invoke("--is711", "--context", "human")
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE" in output or "immediate" in output.lower()

    def test_out_json_written(self, tmp_path: Path) -> None:
        out = tmp_path / "brucella.json"
        result = self._invoke("--is711", "--out-json", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "alert_level" in data

    def test_out_csv_written(self, tmp_path: Path) -> None:
        out = tmp_path / "brucella.csv"
        result = self._invoke("--is711", "--out-csv", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_positive(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"is711_positive": True, "context": "swine"}), encoding="utf-8")
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_input_json_negative(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"is711_positive": False}), encoding="utf-8")
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]


# ===========================================================================
# csfv-risk
# ===========================================================================


class TestCsfvRiskCli:
    """CLI tests for ``lamp-forge csfv-risk``."""

    def _invoke(self, *args: str) -> object:
        return _run("csfv-risk", *args)

    def test_negative_exits_0(self) -> None:
        result = self._invoke("--no-ns5b")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_negative_shows_negative(self) -> None:
        result = self._invoke("--no-ns5b")
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_positive_shows_critical(self) -> None:
        result = self._invoke("--ns5b")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_positive_on_farm_shows_immediate_report(self) -> None:
        result = self._invoke("--ns5b", "--context", "on_farm")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE" in output or "immediate" in output.lower()

    def test_positive_border_post_shows_critical(self) -> None:
        result = self._invoke("--ns5b", "--context", "border_post")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_positive_market_shows_critical(self) -> None:
        result = self._invoke("--ns5b", "--context", "market")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_out_json_written(self, tmp_path: Path) -> None:
        out = tmp_path / "csfv.json"
        result = self._invoke("--ns5b", "--out-json", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "alert_level" in data
        assert data["alert_level"] == "CRITICAL"

    def test_out_csv_written(self, tmp_path: Path) -> None:
        out = tmp_path / "csfv.csv"
        result = self._invoke("--ns5b", "--out-csv", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_positive(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"ns5b_positive": True}), encoding="utf-8")
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_input_json_negative(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"ns5b_positive": False}), encoding="utf-8")
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_output_contains_ns5b_label(self) -> None:
        result = self._invoke("--ns5b")
        assert "NS5B" in result.output or "CSFV" in result.output  # type: ignore[union-attr]


# ===========================================================================
# swine-enteric-risk
# ===========================================================================


class TestSwineEntericRiskCli:
    """CLI tests for ``lamp-forge swine-enteric-risk``."""

    def _invoke(self, *args: str) -> object:
        return _run("swine-enteric-risk", *args)

    def test_all_negative_exits_0(self) -> None:
        result = self._invoke()
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_all_negative_shows_negative(self) -> None:
        result = self._invoke()
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_pedv_shows_critical(self) -> None:
        result = self._invoke("--pedv")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_pdcov_shows_high(self) -> None:
        result = self._invoke("--pdcov")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_rota_a_only_shows_moderate(self) -> None:
        result = self._invoke("--rota-a")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]

    def test_pedv_overrides_pdcov_to_critical(self) -> None:
        result = self._invoke("--pedv", "--pdcov")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_output_shows_target_table(self) -> None:
        result = self._invoke("--pedv", "--rota-a")
        output = result.output  # type: ignore[union-attr]
        assert "PEDV" in output
        assert "ROTA_A" in output or "Rota" in output or "rota" in output.lower()
        assert "[+]" in output

    def test_pedv_shows_quarantine_message(self) -> None:
        result = self._invoke("--pedv")
        output = result.output  # type: ignore[union-attr]
        assert "quarantine" in output.lower() or "lockdown" in output.lower()

    def test_out_json_written(self, tmp_path: Path) -> None:
        out = tmp_path / "swine_enteric.json"
        result = self._invoke("--pedv", "--out-json", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "alert_level" in data
        assert data["alert_level"] == "CRITICAL"

    def test_out_csv_written(self, tmp_path: Path) -> None:
        out = tmp_path / "swine_enteric.csv"
        result = self._invoke("--pdcov", "--out-csv", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_pedv(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"pedv": True, "pdcov": False, "rota_a": False}))
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_input_json_all_negative(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"pedv": False, "pdcov": False, "rota_a": False}))
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]


# ===========================================================================
# poultry-risk
# ===========================================================================


class TestPoultryRiskCli:
    """CLI tests for ``lamp-forge poultry-risk``."""

    def _invoke(self, *args: str) -> object:
        return _run("poultry-risk", *args)

    def test_all_negative_exits_0(self) -> None:
        result = self._invoke()
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_all_negative_shows_negative(self) -> None:
        result = self._invoke()
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    def test_aiv_shows_critical(self) -> None:
        result = self._invoke("--aiv")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_ndv_alone_shows_high(self) -> None:
        result = self._invoke("--ndv")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "HIGH" in result.output  # type: ignore[union-attr]

    def test_ibdv_alone_shows_moderate(self) -> None:
        result = self._invoke("--ibdv")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "MODERATE" in result.output  # type: ignore[union-attr]

    def test_ibv_alone_shows_low(self) -> None:
        result = self._invoke("--ibv")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "LOW" in result.output  # type: ignore[union-attr]

    def test_aiv_overrides_to_critical(self) -> None:
        result = self._invoke("--aiv", "--ndv", "--ibdv", "--ibv")
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_output_shows_target_table(self) -> None:
        result = self._invoke("--aiv", "--ibv")
        output = result.output  # type: ignore[union-attr]
        assert "AIV" in output
        assert "IBV" in output
        assert "[+]" in output

    def test_aiv_shows_immediate_report(self) -> None:
        result = self._invoke("--aiv")
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE" in output or "immediate" in output.lower()

    def test_ndv_shows_immediate_report(self) -> None:
        result = self._invoke("--ndv")
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE" in output or "immediate" in output.lower()

    def test_out_json_written(self, tmp_path: Path) -> None:
        out = tmp_path / "poultry.json"
        result = self._invoke("--aiv", "--out-json", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "alert_level" in data
        assert data["alert_level"] == "CRITICAL"

    def test_out_csv_written(self, tmp_path: Path) -> None:
        out = tmp_path / "poultry.csv"
        result = self._invoke("--ndv", "--out-csv", str(out))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_input_json_aiv(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"aiv": True, "ndv": False, "ibdv": False, "ibv": False}))
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "CRITICAL" in result.output  # type: ignore[union-attr]

    def test_input_json_all_negative(self, tmp_path: Path) -> None:
        flags = tmp_path / "flags.json"
        flags.write_text(json.dumps({"aiv": False, "ndv": False, "ibdv": False, "ibv": False}))
        result = self._invoke("--input-json", str(flags))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "NEGATIVE" in result.output  # type: ignore[union-attr]

    @pytest.mark.parametrize("flag", ["--aiv", "--ndv"])
    def test_notifiable_targets_shows_immediate_report(self, flag: str) -> None:
        result = self._invoke(flag)
        output = result.output  # type: ignore[union-attr]
        assert "IMMEDIATE" in output or "notify" in output.lower()
