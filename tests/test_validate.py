"""Tests for the config validate command and check_config_warnings helper.

All tests are pure Python -- no external binaries required.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from click.testing import CliRunner

from lamp_forge.cli import cli
from lamp_forge.config import check_config_warnings, load_config

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_config(path: Path, data: dict) -> Path:  # type: ignore[type-arg]
    path.write_text(yaml.safe_dump(data))
    return path


def _base_config(tmp_path: Path) -> dict:  # type: ignore[type-arg]
    """Return a minimal valid config dict."""
    return {
        "target": {
            "name": "test_target",
            "taxon_id": 1773,
            "gene": "rpoB",
            "max_sequences": 20,
            "email": "user@test.edu",
        },
        "off_targets": {
            "fasta_dir": str(tmp_path / "off_targets"),
            "min_identity_threshold": 0.80,
            "min_coverage_threshold": 0.80,
        },
        "conservation": {
            "window_size": 30,
            "entropy_threshold": 0.25,
            "min_region_length": 250,
        },
        "primer": {
            "tm_min": 60.0,
            "tm_max": 65.0,
            "tm_match_tolerance": 2.0,
            "gc_min": 40.0,
            "gc_max": 65.0,
            "hairpin_dg_threshold": -2.0,
            "dimer_dg_threshold": -5.0,
            "amplicon_size": {"f2_b2_min": 120, "f2_b2_max": 160},
        },
        "output": {
            "dir": "results",
            "top_n": 10,
            "generate_html": True,
            "generate_csv": True,
        },
    }


# ---------------------------------------------------------------------------
# check_config_warnings -- directory check
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsOffTargetDir:
    def test_no_warning_when_dir_exists(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        off_dir = tmp_path / "off_targets"
        off_dir.mkdir()
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=True)
        assert not any("[off_targets] fasta_dir" in w for w in warnings)

    def test_warns_when_dir_missing(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        # off_targets dir is NOT created -> should warn
        warnings = check_config_warnings(cfg, check_dirs=True)
        assert any("[off_targets]" in w and "does not exist" in w for w in warnings)

    def test_no_dir_warning_when_check_dirs_false(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("does not exist" in w for w in warnings)


# ---------------------------------------------------------------------------
# check_config_warnings -- entropy threshold
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsEntropy:
    def test_no_warning_for_normal_entropy(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[conservation]" in w for w in warnings)

    def test_warns_when_entropy_too_low(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.05
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[conservation]" in w and "tight" in w for w in warnings)

    def test_warns_when_entropy_too_high(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.55
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[conservation]" in w and "loose" in w for w in warnings)

    def test_boundary_exactly_008_no_warning(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.08
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[conservation]" in w and "tight" in w for w in warnings)

    def test_boundary_exactly_050_no_warning(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.50
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[conservation]" in w and "loose" in w for w in warnings)


# ---------------------------------------------------------------------------
# check_config_warnings -- max_sequences
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsMaxSequences:
    def test_no_warning_for_sufficient_sequences(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[target]" in w and "max_sequences" in w for w in warnings)

    def test_warns_when_max_sequences_too_low(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["target"]["max_sequences"] = 3
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[target]" in w and "max_sequences" in w for w in warnings)

    def test_boundary_5_no_warning(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["target"]["max_sequences"] = 5
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[target]" in w and "max_sequences" in w for w in warnings)


# ---------------------------------------------------------------------------
# check_config_warnings -- Tm range
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsTm:
    def test_no_warning_for_normal_tm(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[primer]" in w and "degC" in w for w in warnings)

    def test_warns_when_tm_min_too_low(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["primer"]["tm_min"] = 55.0
        d["primer"]["tm_max"] = 60.0
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[primer]" in w and "tm_min" in w for w in warnings)

    def test_warns_when_tm_max_too_high(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["primer"]["tm_min"] = 68.0
        d["primer"]["tm_max"] = 75.0
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[primer]" in w and "tm_max" in w for w in warnings)

    def test_boundary_tm_min_58_no_warning(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["primer"]["tm_min"] = 58.0
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[primer]" in w and "tm_min" in w for w in warnings)

    def test_boundary_tm_max_72_no_warning(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["primer"]["tm_min"] = 68.0
        d["primer"]["tm_max"] = 72.0
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[primer]" in w and "tm_max" in w for w in warnings)


# ---------------------------------------------------------------------------
# check_config_warnings -- top_n
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsTopN:
    def test_no_warning_for_normal_top_n(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[output]" in w for w in warnings)

    def test_warns_when_top_n_very_large(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["output"]["top_n"] = 100
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[output]" in w and "top_n" in w for w in warnings)

    def test_boundary_top_n_50_no_warning(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["output"]["top_n"] = 50
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[output]" in w and "top_n" in w for w in warnings)


# ---------------------------------------------------------------------------
# check_config_warnings -- identity threshold
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsIdentity:
    def test_no_warning_for_normal_identity(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert not any("[off_targets]" in w and "identity" in w for w in warnings)

    def test_warns_when_identity_threshold_very_low(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["off_targets"]["min_identity_threshold"] = 0.60
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        assert any("[off_targets]" in w and "identity" in w for w in warnings)

    def test_boundary_070_warns(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["off_targets"]["min_identity_threshold"] = 0.70
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=False)
        # exactly 0.70 is the boundary -- not triggered (< 0.70)
        assert not any("[off_targets]" in w and "identity" in w for w in warnings)

    def test_clean_config_returns_empty(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        off_dir = tmp_path / "off_targets"
        off_dir.mkdir()
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        warnings = check_config_warnings(cfg, check_dirs=True)
        assert warnings == []


# ---------------------------------------------------------------------------
# check_config_warnings -- returns list type
# ---------------------------------------------------------------------------


class TestCheckConfigWarningsReturnType:
    def test_returns_list(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        result = check_config_warnings(cfg, check_dirs=False)
        assert isinstance(result, list)

    def test_warnings_are_strings(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.60
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        result = check_config_warnings(cfg, check_dirs=False)
        assert all(isinstance(w, str) for w in result)

    def test_multiple_problems_all_reported(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.60
        d["target"]["max_sequences"] = 2
        cfg = load_config(_write_config(tmp_path / "cfg.yaml", d))
        result = check_config_warnings(cfg, check_dirs=False)
        assert len(result) >= 2


# ---------------------------------------------------------------------------
# CLI: lamp-forge validate
# ---------------------------------------------------------------------------


class TestValidateCli:
    def _write_valid(self, tmp_path: Path) -> Path:
        d = _base_config(tmp_path)
        off_dir = tmp_path / "off_targets"
        off_dir.mkdir()
        return _write_config(tmp_path / "cfg.yaml", d)

    def test_exits_zero_for_valid_config(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert result.exit_code == 0, result.output

    def test_output_contains_ok_marker(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "Config OK" in result.output

    def test_output_contains_target_name(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "test_target" in result.output

    def test_output_contains_tm_range(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "60.0" in result.output
        assert "65.0" in result.output

    def test_output_contains_off_target_section(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "Off-target" in result.output

    def test_exits_two_for_invalid_config(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.yaml"
        bad.write_text("not_a_real_config: true")
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", "--config", str(bad)], obj={})
        assert result.exit_code == 2

    def test_invalid_config_shows_error(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.yaml"
        bad.write_text("not_a_real_config: true")
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", "--config", str(bad)], obj={})
        assert "Config error" in result.output or "Missing required" in result.output

    def test_soft_warning_shown_for_missing_dir(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg_path = _write_config(tmp_path / "cfg.yaml", d)
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        # No dir created -> warning expected, but exit 0
        assert result.exit_code == 0
        assert "does not exist" in result.output

    def test_no_dir_warning_with_flag(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        cfg_path = _write_config(tmp_path / "cfg.yaml", d)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["validate", "--config", str(cfg_path), "--no-check-dirs"],
            obj={},
        )
        assert result.exit_code == 0
        assert "does not exist" not in result.output

    def test_no_warnings_line_shown_when_clean(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "No warnings" in result.output

    def test_warnings_count_shown(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.60
        d["target"]["max_sequences"] = 2
        cfg_path = _write_config(tmp_path / "cfg.yaml", d)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["validate", "--config", str(cfg_path), "--no-check-dirs"],
            obj={},
        )
        assert result.exit_code == 0
        assert "Warnings" in result.output

    def test_accessions_only_config_shows_accession_count(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        del d["target"]["taxon_id"]
        d["target"]["accessions"] = ["OQ518516.1", "OR351212.1"]
        cfg_path = _write_config(tmp_path / "cfg.yaml", d)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["validate", "--config", str(cfg_path), "--no-check-dirs"],
            obj={},
        )
        assert result.exit_code == 0
        assert "accession" in result.output

    def test_validate_help_text(self) -> None:
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", "--help"], obj={})
        assert result.exit_code == 0
        assert "validate" in result.output.lower() or "Validate" in result.output

    def test_html_and_csv_flags_shown(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "HTML=yes" in result.output
        assert "CSV=yes" in result.output

    def test_entropy_warning_shown_when_high(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["conservation"]["entropy_threshold"] = 0.60
        cfg_path = _write_config(tmp_path / "cfg.yaml", d)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["validate", "--config", str(cfg_path), "--no-check-dirs"],
            obj={},
        )
        assert result.exit_code == 0
        assert "[conservation]" in result.output

    def test_gene_shown_in_source_line(self, tmp_path: Path) -> None:
        runner = CliRunner()
        cfg_path = self._write_valid(tmp_path)
        result = runner.invoke(cli, ["validate", "--config", str(cfg_path)], obj={})
        assert "rpoB" in result.output

    def test_taxon_and_accessions_both_shown(self, tmp_path: Path) -> None:
        d = _base_config(tmp_path)
        d["target"]["accessions"] = ["OQ518516.1"]
        cfg_path = _write_config(tmp_path / "cfg.yaml", d)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["validate", "--config", str(cfg_path), "--no-check-dirs"],
            obj={},
        )
        assert result.exit_code == 0
        assert "1773" in result.output
        assert "accession" in result.output
