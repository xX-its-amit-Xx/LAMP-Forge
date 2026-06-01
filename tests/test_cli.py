"""CLI smoke tests using Click's test runner."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

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


# ---------------------------------------------------------------------------
# `lamp-forge cartridge` end-to-end smoke test
# ---------------------------------------------------------------------------


def _primer(role: str, seq: str, tm: float = 64.0) -> dict[str, Any]:
    return {
        "role": role,
        "sequence": seq,
        "start": 0,
        "end": len(seq),
        "strand": "+",
        "tm": tm,
        "gc_percent": 50.0,
        "hairpin_dg": -1.0,
        "homodimer_dg": -2.0,
        "inner_subseq": None,
        "outer_subseq": None,
    }


def _write_sets_json(path: Path, set_id: str, seqs: dict[str, str]) -> None:
    entry: dict[str, Any] = {
        "set_id": set_id,
        "region_id": "R000",
        "composite_score": 0.9,
    }
    for field in ("f3", "b3", "fip", "bip", "lf", "lb"):
        entry[field] = _primer(field.upper(), seqs[field]) if field in seqs else None
    path.write_text(
        json.dumps(
            {
                "version": "0.1.0",
                "generated_at": "2026-01-01T00:00:00Z",
                "primer_sets": [entry],
            }
        ),
        encoding="utf-8",
    )


def _two_channel_fixture(tmp_path: Path) -> tuple[Path, Path]:
    a = tmp_path / "rna.json"
    b = tmp_path / "dna.json"
    _write_sets_json(
        a,
        "RNA1_set_001",
        {
            "f3": "GACCTGACATTCGGATCACG",
            "b3": "TTGCAGTCGAACTGGTACGA",
            "fip": "CAGTTCGACTGCAAGGTACCGATCGAATGTCAGGTC",
            "bip": "GTACCAGTTCGACTGCAATGCGTGATCCGAATGTCA",
        },
    )
    _write_sets_json(
        b,
        "DNA1_set_001",
        {
            "f3": "AAACCCGGGTTTACGTAACG",
            "b3": "CCTTGGAATTGCCAATGCGT",
            "fip": "TTAGGCCAATTCCGGAAACCCGGGAATTAGGCCAAT",
            "bip": "GGGTTAACCGGAATTAGGCCAATTAGGGTTAACCGG",
        },
    )
    return a, b


def _write_cartridge_yaml(
    path: Path,
    *,
    rna_json: Path,
    dna_json: Path,
    max_channels: int = 4,
) -> None:
    payload = {
        "cartridge_id": "cli_test_cart",
        "chemistry": "rt-lamp",
        "max_channels": max_channels,
        "extraction": {
            "sample_volume_ul": 1000.0,
            "extraction_efficiency": 0.50,
            "eluate_volume_ul": 100.0,
            "reaction_input_ul": 5.0,
        },
        "device_window_min": 60.0,
        "pool": {"stock_conc_um": 100.0, "total_pool_volume_ul": 500.0},
        "vendor": "idt",
        "targets": [
            {"label": "RNA1", "sets_json": str(rna_json), "target_na": "rna"},
            {"label": "DNA1", "sets_json": str(dna_json), "target_na": "dna"},
        ],
    }
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


_EXPECTED = (
    "cartridge_manifest.json",
    "panel.json",
    "panel_primers.csv",
    "pool_plan.csv",
    "lod_panel.csv",
    "rt_check.csv",
    "vendor_order.csv",
    "build_report.html",
)


def test_cartridge_command_go_writes_all_artifacts(tmp_path: Path) -> None:
    rna, dna = _two_channel_fixture(tmp_path)
    manifest = tmp_path / "cart.yaml"
    _write_cartridge_yaml(manifest, rna_json=rna, dna_json=dna)
    out_dir = tmp_path / "build"

    runner = CliRunner()
    result = runner.invoke(cli, ["cartridge", str(manifest), "--out-dir", str(out_dir)], obj={})
    assert result.exit_code == 0, result.output
    for name in _EXPECTED:
        assert (out_dir / name).exists(), f"missing {name}"
    payload = json.loads((out_dir / "cartridge_manifest.json").read_text())
    expected_keys = {
        "schema_version",
        "cartridge_id",
        "build_status",
        "manifest_hash",
        "per_target",
        "panel_summary",
        "panel_lod_worst",
    }
    assert expected_keys.issubset(set(payload.keys()))


def test_cartridge_command_no_go_returns_nonzero(tmp_path: Path) -> None:
    rna, dna = _two_channel_fixture(tmp_path)
    manifest = tmp_path / "cart.yaml"
    # max_channels=1 < 2 targets => load_cartridge_manifest raises;
    # the CLI maps that to exit 2.
    _write_cartridge_yaml(manifest, rna_json=rna, dna_json=dna, max_channels=1)
    out_dir = tmp_path / "build"

    runner = CliRunner()
    result = runner.invoke(cli, ["cartridge", str(manifest), "--out-dir", str(out_dir)], obj={})
    assert result.exit_code != 0
    assert "Manifest error" in result.output
