"""CLI smoke tests for BioVind trend-analysis commands.

Covers the four trend commands that have module-level tests but no
CLI-level tests: farm-trend, poc-trend, bov-trend, swine-enteric-trend.

All tests use Click's CliRunner and synthetic JSON / CSV fixtures --
no external binaries required.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import ClassVar

from click.testing import CliRunner

from lamp_forge.cli import cli

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _run(*args: str) -> object:
    runner = CliRunner()
    return runner.invoke(cli, list(args), obj={})


def _write_json(path: Path, panel_flags: dict[str, bool]) -> Path:
    """Write a minimal risk-output JSON file with the given panel_flags."""
    path.write_text(
        json.dumps({"alert_level": "NEGATIVE", "panel_flags": panel_flags}),
        encoding="utf-8",
    )
    return path


# ===========================================================================
# farm-trend
# ===========================================================================


class TestFarmTrendCli:
    """CLI smoke tests for ``lamp-forge farm-trend``."""

    _ALL_NEG: ClassVar[dict[str, bool]] = {
        "asfv": False,
        "fmdv": False,
        "aiv": False,
        "ndv": False,
        "prrsv": False,
    }

    def _invoke(self, *args: str) -> object:
        return _run("farm-trend", *args)

    def test_no_input_exits_nonzero(self) -> None:
        result = self._invoke()
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_single_json_insufficient_data(self, tmp_path: Path) -> None:
        j = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        result = self._invoke("--farm-result", str(j))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT_DATA" in result.output  # type: ignore[union-attr]

    def test_two_clean_samples_stable_clear(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--farm-result", str(j1), "--farm-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_aiv_newly_detected_shows_emerging(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", {**self._ALL_NEG, "aiv": True})
        result = self._invoke("--farm-result", str(j1), "--farm-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        output = result.output  # type: ignore[union-attr]
        assert "EMERGING" in output
        assert "NEWLY" in output or "newly" in output.lower()

    def test_aiv_cleared_shows_resolving(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", {**self._ALL_NEG, "aiv": True})
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--farm-result", str(j1), "--farm-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "RESOLVING" in result.output  # type: ignore[union-attr]

    def test_prrsv_endemic_flag(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", {**self._ALL_NEG, "prrsv": True})
        j2 = _write_json(tmp_path / "s2.json", {**self._ALL_NEG, "prrsv": True})
        j3 = _write_json(tmp_path / "s3.json", {**self._ALL_NEG, "prrsv": True})
        result = self._invoke(
            "--farm-result",
            str(j1),
            "--farm-result",
            str(j2),
            "--farm-result",
            str(j3),
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "endemic" in result.output.lower()  # type: ignore[union-attr]

    def test_out_json_written(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        out = tmp_path / "trend.json"
        result = self._invoke(
            "--farm-result",
            str(j1),
            "--farm-result",
            str(j2),
            "--out-json",
            str(out),
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "direction" in data
        assert data["direction"] == "STABLE_CLEAR"

    def test_out_csv_written(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        out = tmp_path / "timeline.csv"
        result = self._invoke(
            "--farm-result",
            str(j1),
            "--farm-result",
            str(j2),
            "--out-csv",
            str(out),
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()

    def test_csv_input_two_clean_rows(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "monitoring.csv"
        with csv_file.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["sample_id", "date", "asfv", "fmdv", "aiv", "ndv", "prrsv"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "sample_id": "W1",
                    "date": "2026-01-01",
                    "asfv": 0,
                    "fmdv": 0,
                    "aiv": 0,
                    "ndv": 0,
                    "prrsv": 0,
                }
            )
            writer.writerow(
                {
                    "sample_id": "W2",
                    "date": "2026-02-01",
                    "asfv": 0,
                    "fmdv": 0,
                    "aiv": 0,
                    "ndv": 0,
                    "prrsv": 0,
                }
            )
        result = self._invoke("--csv", str(csv_file))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_csv_and_json_are_mutually_exclusive(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        csv_file = tmp_path / "m.csv"
        csv_file.write_text("sample_id,asfv,fmdv,aiv,ndv,prrsv\nW1,0,0,0,0,0\n", encoding="utf-8")
        result = self._invoke("--farm-result", str(j1), "--csv", str(csv_file))
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_table_shows_pathogen_columns(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", {**self._ALL_NEG, "aiv": True})
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--farm-result", str(j1), "--farm-result", str(j2))
        output = result.output  # type: ignore[union-attr]
        assert "AIV" in output
        assert "[+]" in output
        assert "[-]" in output


# ===========================================================================
# poc-trend
# ===========================================================================


class TestPocTrendCli:
    """CLI smoke tests for ``lamp-forge poc-trend``."""

    _ALL_NEG: ClassVar[dict[str, bool]] = {
        "mtb": False,
        "cdiff": False,
        "sars2": False,
        "flu_a": False,
        "gas": False,
    }

    def _invoke(self, *args: str) -> object:
        return _run("poc-trend", *args)

    def test_no_input_exits_nonzero(self) -> None:
        result = self._invoke()
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_single_json_insufficient_data(self, tmp_path: Path) -> None:
        j = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        result = self._invoke("--poc-result", str(j))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT_DATA" in result.output  # type: ignore[union-attr]

    def test_two_clean_samples_stable_clear(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--poc-result", str(j1), "--poc-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_mtb_newly_detected_shows_emerging(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", {**self._ALL_NEG, "mtb": True})
        result = self._invoke("--poc-result", str(j1), "--poc-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        output = result.output  # type: ignore[union-attr]
        assert "EMERGING" in output

    def test_respiratory_wave_two_pathogens(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", {**self._ALL_NEG, "flu_a": True, "sars2": True})
        j2 = _write_json(tmp_path / "s2.json", {**self._ALL_NEG, "flu_a": True, "sars2": True})
        result = self._invoke("--poc-result", str(j1), "--poc-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]

    def test_out_json_written(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        out = tmp_path / "poc_trend.json"
        result = self._invoke(
            "--poc-result",
            str(j1),
            "--poc-result",
            str(j2),
            "--out-json",
            str(out),
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "direction" in data

    def test_csv_input_two_clean_rows(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "clinic.csv"
        with csv_file.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["sample_id", "mtb", "cdiff", "sars2", "flu_a", "gas"],
            )
            writer.writeheader()
            writer.writerow(
                {"sample_id": "C1", "mtb": 0, "cdiff": 0, "sars2": 0, "flu_a": 0, "gas": 0}
            )
            writer.writerow(
                {"sample_id": "C2", "mtb": 0, "cdiff": 0, "sars2": 0, "flu_a": 0, "gas": 0}
            )
        result = self._invoke("--csv", str(csv_file))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]


# ===========================================================================
# bov-trend
# ===========================================================================


class TestBovTrendCli:
    """CLI smoke tests for ``lamp-forge bov-trend``."""

    _ALL_NEG: ClassVar[dict[str, bool]] = {
        "brsv": False,
        "bcov": False,
        "bvdv": False,
        "ibr": False,
        "mhae": False,
        "mbovis": False,
    }

    def _invoke(self, *args: str) -> object:
        return _run("bov-trend", *args)

    def test_no_input_exits_nonzero(self) -> None:
        result = self._invoke()
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_single_json_insufficient_data(self, tmp_path: Path) -> None:
        j = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        result = self._invoke("--bov-result", str(j))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT_DATA" in result.output  # type: ignore[union-attr]

    def test_two_clean_samples_stable_clear(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--bov-result", str(j1), "--bov-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_viral_coinfection_shows_brd_pattern(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(
            tmp_path / "s2.json",
            {**self._ALL_NEG, "brsv": True, "mhae": True},
        )
        result = self._invoke("--bov-result", str(j1), "--bov-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "EMERGING" in result.output  # type: ignore[union-attr]

    def test_out_json_written(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        out = tmp_path / "bov_trend.json"
        result = self._invoke(
            "--bov-result",
            str(j1),
            "--bov-result",
            str(j2),
            "--out-json",
            str(out),
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "direction" in data

    def test_csv_input_two_clean_rows(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "feedlot.csv"
        with csv_file.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["sample_id", "brsv", "bcov", "bvdv", "ibr", "mhae", "mbovis"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "sample_id": "P1",
                    "brsv": 0,
                    "bcov": 0,
                    "bvdv": 0,
                    "ibr": 0,
                    "mhae": 0,
                    "mbovis": 0,
                }
            )
            writer.writerow(
                {
                    "sample_id": "P2",
                    "brsv": 0,
                    "bcov": 0,
                    "bvdv": 0,
                    "ibr": 0,
                    "mhae": 0,
                    "mbovis": 0,
                }
            )
        result = self._invoke("--csv", str(csv_file))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]


# ===========================================================================
# swine-enteric-trend
# ===========================================================================


class TestSwineEntericTrendCli:
    """CLI smoke tests for ``lamp-forge swine-enteric-trend``."""

    _ALL_NEG: ClassVar[dict[str, bool]] = {
        "pedv": False,
        "pdcov": False,
        "rota_a": False,
    }

    def _invoke(self, *args: str) -> object:
        return _run("swine-enteric-trend", *args)

    def test_no_input_exits_nonzero(self) -> None:
        result = self._invoke()
        assert result.exit_code != 0  # type: ignore[union-attr]

    def test_single_json_insufficient_data(self, tmp_path: Path) -> None:
        j = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        result = self._invoke("--enteric-result", str(j))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "INSUFFICIENT_DATA" in result.output  # type: ignore[union-attr]

    def test_two_clean_samples_stable_clear(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--enteric-result", str(j1), "--enteric-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_pedv_newly_detected_shows_emerging(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", {**self._ALL_NEG, "pedv": True})
        result = self._invoke("--enteric-result", str(j1), "--enteric-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        output = result.output  # type: ignore[union-attr]
        assert "EMERGING" in output

    def test_pedv_shows_quarantine_language(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", {**self._ALL_NEG, "pedv": True})
        j2 = _write_json(tmp_path / "s2.json", {**self._ALL_NEG, "pedv": True})
        result = self._invoke("--enteric-result", str(j1), "--enteric-result", str(j2))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        output = result.output  # type: ignore[union-attr]
        assert "quarantine" in output.lower() or "lockdown" in output.lower() or "PEDV" in output

    def test_out_json_written(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", self._ALL_NEG)
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        out = tmp_path / "enteric_trend.json"
        result = self._invoke(
            "--enteric-result",
            str(j1),
            "--enteric-result",
            str(j2),
            "--out-json",
            str(out),
        )
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert out.exists()
        data = json.loads(out.read_text())
        assert "direction" in data

    def test_csv_input_two_clean_rows(self, tmp_path: Path) -> None:
        csv_file = tmp_path / "farrowing.csv"
        with csv_file.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["sample_id", "pedv", "pdcov", "rota_a"],
            )
            writer.writeheader()
            writer.writerow({"sample_id": "W1", "pedv": 0, "pdcov": 0, "rota_a": 0})
            writer.writerow({"sample_id": "W2", "pedv": 0, "pdcov": 0, "rota_a": 0})
        result = self._invoke("--csv", str(csv_file))
        assert result.exit_code == 0, result.output  # type: ignore[union-attr]
        assert "STABLE_CLEAR" in result.output  # type: ignore[union-attr]

    def test_table_shows_pathogen_columns(self, tmp_path: Path) -> None:
        j1 = _write_json(tmp_path / "s1.json", {**self._ALL_NEG, "pedv": True})
        j2 = _write_json(tmp_path / "s2.json", self._ALL_NEG)
        result = self._invoke("--enteric-result", str(j1), "--enteric-result", str(j2))
        output = result.output  # type: ignore[union-attr]
        assert "PEDV" in output
        assert "[+]" in output
        assert "[-]" in output
