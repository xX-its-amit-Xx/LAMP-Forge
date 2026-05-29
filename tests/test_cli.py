"""CLI smoke tests using Click's test runner."""

from __future__ import annotations

import pytest
import yaml
from click.testing import CliRunner

from lamp_forge.cli import cli


def test_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "LAMP-Forge" in result.output


def test_version() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["version"])
    assert result.exit_code == 0
    assert "lamp-forge" in result.output


def test_run_with_bad_config_exits_with_error(tmp_path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("not_a_real_config: true")
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--config", str(bad)])
    assert result.exit_code != 0
    assert "Config error" in result.output or "Missing required field" in result.output


def test_run_missing_config_file() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--config", "/no/such/file.yaml"])
    assert result.exit_code != 0
