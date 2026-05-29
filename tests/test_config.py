"""Tests for YAML config loading and validation."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lamp_forge.config import ConfigError, load_config


def _write(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.fixture
def minimal_config_dict() -> dict:
    return {
        "target": {
            "name": "test",
            "taxon_id": 1773,
            "gene": "rpoB",
            "max_sequences": 10,
            "email": "user@test.edu",
        },
        "off_targets": {
            "fasta_dir": "off_targets",
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


class TestLoadConfig:
    def test_minimal_loads(self, tmp_path: Path, minimal_config_dict: dict) -> None:
        cfg = load_config(_write(tmp_path / "c.yaml", minimal_config_dict))
        assert cfg.target_name == "test"
        assert cfg.taxon_id == 1773
        assert cfg.tm_min == 60.0
        assert cfg.f2_b2_max == 160

    def test_missing_file(self) -> None:
        with pytest.raises(ConfigError, match="not found"):
            load_config(Path("/nope/never/exists.yaml"))

    def test_missing_section(self, tmp_path: Path, minimal_config_dict: dict) -> None:
        del minimal_config_dict["primer"]
        with pytest.raises(ConfigError, match="primer"):
            load_config(_write(tmp_path / "c.yaml", minimal_config_dict))

    def test_tm_min_above_max_rejected(
        self, tmp_path: Path, minimal_config_dict: dict
    ) -> None:
        minimal_config_dict["primer"]["tm_min"] = 70.0
        with pytest.raises(ConfigError, match="tm_min"):
            load_config(_write(tmp_path / "c.yaml", minimal_config_dict))

    def test_no_target_id_or_accessions_rejected(
        self, tmp_path: Path, minimal_config_dict: dict
    ) -> None:
        minimal_config_dict["target"].pop("taxon_id")
        with pytest.raises(ConfigError, match="taxon_id"):
            load_config(_write(tmp_path / "c.yaml", minimal_config_dict))

    def test_email_required_unless_env(
        self,
        tmp_path: Path,
        minimal_config_dict: dict,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        minimal_config_dict["target"].pop("email")
        monkeypatch.delenv("NCBI_EMAIL", raising=False)
        with pytest.raises(ConfigError, match="email"):
            load_config(_write(tmp_path / "c.yaml", minimal_config_dict))

    def test_email_from_env(
        self,
        tmp_path: Path,
        minimal_config_dict: dict,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        minimal_config_dict["target"].pop("email")
        monkeypatch.setenv("NCBI_EMAIL", "envuser@test.edu")
        cfg = load_config(_write(tmp_path / "c.yaml", minimal_config_dict))
        assert cfg.email == "envuser@test.edu"

    def test_region_too_short_for_amplicon_rejected(
        self, tmp_path: Path, minimal_config_dict: dict
    ) -> None:
        minimal_config_dict["conservation"]["min_region_length"] = 50
        with pytest.raises(ConfigError, match="min_region_length"):
            load_config(_write(tmp_path / "c.yaml", minimal_config_dict))

    def test_f2_b2_too_narrow_rejected(
        self, tmp_path: Path, minimal_config_dict: dict
    ) -> None:
        minimal_config_dict["primer"]["amplicon_size"]["f2_b2_min"] = 50
        with pytest.raises(ConfigError, match="LAMP geometric floor"):
            load_config(_write(tmp_path / "c.yaml", minimal_config_dict))
