"""Tests for the 6-channel oilfield MIC/souring BioID cartridge manifest.

Verifies that:
- The production YAML is syntactically valid and has the expected top-level
  structure (parsed as raw YAML, without requiring primer_sets.json to exist).
- A synthetic 6-channel all-DNA manifest round-trips through
  load_cartridge_manifest correctly.
- Chemistry consistency: all channels are DNA, chemistry is "lamp".
- Pool stoichiometry: stock_conc_um meets the 6-channel minimum (6 × 44 µM).
- Channel labels match the five MIC guilds plus the 16S control.
- The control channel is labelled with role="control"; all others are "detect".
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from lamp_forge.cartridge import CartridgeManifest, load_cartridge_manifest
from lamp_forge.lod import ExtractionParams
from lamp_forge.pool import PoolingParams

# ---------------------------------------------------------------------------
# Constants matching the production manifest
# ---------------------------------------------------------------------------

_PRODUCTION_YAML = (
    Path(__file__).parent.parent / "config" / "cartridges" / "cartridge_mic_oilfield_6ch.yaml"
)

# Canonical MIC-guild channel labels in panel order (detect first, control last)
_EXPECTED_LABELS = ("SRB_dsrB", "MCR_mcrA", "IRB_omcA", "APB_fthfs", "NRB_narG", "CTRL_16S")

# Minimum pool stock for 6 channels: 6 × 44 µM = 264 µM
_MIN_STOCK_UM = 6 * 44

# ---------------------------------------------------------------------------
# Shared fixture helpers (mirrors test_cartridge.py conventions)
# ---------------------------------------------------------------------------


def _primer_dict(role: str, seq: str, tm: float = 64.0) -> dict[str, Any]:
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


def _write_sets_json(path: Path, set_id: str, *, seqs: dict[str, str], tm: float = 64.0) -> None:
    entry: dict[str, Any] = {
        "set_id": set_id,
        "region_id": "R000",
        "composite_score": 0.9,
    }
    for field in ("f3", "b3", "fip", "bip", "lf", "lb"):
        seq = seqs.get(field)
        entry[field] = _primer_dict(field.upper(), seq, tm=tm) if seq else None
    payload = {
        "version": "0.1.0",
        "generated_at": "2026-01-01T00:00:00Z",
        "primer_sets": [entry],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


# Six distinct non-clashing primer sequence sets (one per MIC-guild channel).
# Sequences are synthetic; only distinctness matters for the cartridge loader.
_CHANNEL_SEQS: list[dict[str, str]] = [
    {  # SRB_dsrB
        "f3": "AACCTTGGCAGTACGATCGT",
        "b3": "TTGGCAATCGGTACGAACCT",
        "fip": "CAGTCGAACTTGGCAATCGGAACCTTGGCAGTACGA",
        "bip": "GTTCGACTGCCAAGGTTAGGTTGGCAATCGGTACGA",
    },
    {  # MCR_mcrA
        "f3": "GCTAGCATGTACCGATCGAC",
        "b3": "ATCGCTAGCATGTACCGTAC",
        "fip": "TACGATCGGCTAGCATGTACGATCGACGCTAGCATG",
        "bip": "CATGCTAGCGATCGTATGCTATCGCTAGCATGTACC",
    },
    {  # IRB_omcA
        "f3": "TTACGGCAATCGTACGAGCT",
        "b3": "AGCTACGGCAATCGTACGAG",
        "fip": "TACGAGCTTACGGCAATCGTTTACGGCAATCGTACG",
        "bip": "CGAGCTAGCTACGGCAATCGAGCTACGGCAATCGTA",
    },
    {  # APB_fthfs
        "f3": "CCATGGACTGCATCGTACGA",
        "b3": "CGTACGACCATGGACTGCAT",
        "fip": "GCATCGTACGACCATGGACTCCATGGACTGCATCGTA",
        "bip": "CTGCATCCGTACGACCATGGCGTACGACCATGGACTG",
    },
    {  # NRB_narG
        "f3": "GGATTCAGTCGCATACGGTA",
        "b3": "CGGTAATCAGTCGCATACGG",
        "fip": "GCATACGGTAATCAGTCGCAGGATTCAGTCGCATACG",
        "bip": "CATCGCATACGGTATCAGTCCGGTAATCAGTCGCATA",
    },
    {  # CTRL_16S
        "f3": "AAACCCGGGTTTACGTAACG",
        "b3": "CCTTGGAATTGCCAATGCGT",
        "fip": "TTAGGCCAATTCCGGAAACCCGGGAATTAGGCCAAT",
        "bip": "GGGTTAACCGGAATTAGGCCAATTAGGGTTAACCGG",
    },
]


def _make_six_channel_fixture(tmp_path: Path) -> list[Path]:
    """Write six distinct primer_sets.json files for a 6-channel DNA cartridge."""
    paths: list[Path] = []
    for label, seqs in zip(_EXPECTED_LABELS, _CHANNEL_SEQS, strict=True):
        p = tmp_path / f"{label.lower()}_sets.json"
        _write_sets_json(p, f"{label}_set_001", seqs=seqs)
        paths.append(p)
    return paths


def _write_mic_manifest_yaml(
    path: Path,
    target_paths: list[Path],
    *,
    cartridge_id: str = "mic_oilfield_6ch_v1",
    chemistry: str = "lamp",
    max_channels: int = 6,
    stock_conc_um: float = 300.0,
    total_pool_volume_ul: float = 500.0,
) -> None:
    """Write a 6-channel MIC oilfield manifest YAML for tests."""
    roles = ["detect", "detect", "detect", "detect", "detect", "control"]
    targets = [
        {
            "label": label,
            "sets_json": str(p),
            "target_na": "dna",
            "role": role,
        }
        for label, p, role in zip(_EXPECTED_LABELS, target_paths, roles, strict=True)
    ]
    payload: dict[str, Any] = {
        "cartridge_id": cartridge_id,
        "chemistry": chemistry,
        "max_channels": max_channels,
        "extraction": {
            "sample_volume_ul": 1000.0,
            "extraction_efficiency": 0.50,
            "eluate_volume_ul": 100.0,
            "reaction_input_ul": 5.0,
        },
        "device_window_min": 60.0,
        "pool": {
            "stock_conc_um": stock_conc_um,
            "total_pool_volume_ul": total_pool_volume_ul,
        },
        "vendor": "idt",
        "targets": targets,
    }
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


# ---------------------------------------------------------------------------
# Production YAML structure tests (raw YAML parse, no file-existence checks)
# ---------------------------------------------------------------------------


class TestProductionYamlStructure:
    """Validate the committed cartridge YAML without running the full pipeline."""

    def test_production_yaml_is_parseable(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert isinstance(raw, dict)

    def test_cartridge_id(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert raw["cartridge_id"] == "mic_oilfield_6ch_v1"

    def test_chemistry_is_lamp(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert raw["chemistry"] == "lamp"

    def test_max_channels_is_six(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert raw["max_channels"] == 6

    def test_exactly_six_targets(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert len(raw["targets"]) == 6

    def test_target_labels_match_mic_guilds(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        labels = tuple(t["label"] for t in raw["targets"])
        assert labels == _EXPECTED_LABELS

    def test_all_targets_are_dna(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        for t in raw["targets"]:
            assert t["target_na"] == "dna", f"Expected DNA for {t['label']}, got {t['target_na']}"

    def test_control_channel_role(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        ctrl = raw["targets"][-1]
        assert ctrl["label"] == "CTRL_16S"
        assert ctrl["role"] == "control"

    def test_detect_channel_roles(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        for t in raw["targets"][:-1]:
            assert t["role"] == "detect", f"Expected 'detect' for {t['label']}"

    def test_pool_stock_meets_minimum(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert raw["pool"]["stock_conc_um"] >= _MIN_STOCK_UM

    def test_extraction_fields_present(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        ext = raw["extraction"]
        for key in (
            "sample_volume_ul",
            "extraction_efficiency",
            "eluate_volume_ul",
            "reaction_input_ul",
        ):
            assert key in ext, f"extraction missing key: {key}"

    def test_vendor_is_idt(self) -> None:
        raw = yaml.safe_load(_PRODUCTION_YAML.read_text(encoding="utf-8"))
        assert raw["vendor"] == "idt"


# ---------------------------------------------------------------------------
# Synthetic manifest load tests (load_cartridge_manifest with real file I/O)
# ---------------------------------------------------------------------------


class TestSyntheticManifestLoad:
    """Exercise load_cartridge_manifest with six synthetic primer_sets.json files."""

    def test_loads_six_channel_dna_manifest(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        assert isinstance(manifest, CartridgeManifest)
        assert manifest.cartridge_id == "mic_oilfield_6ch_v1"
        assert manifest.chemistry == "lamp"
        assert manifest.max_channels == 6
        assert len(manifest.targets) == 6

    def test_all_channels_are_dna(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        for spec in manifest.targets:
            assert spec.target_na == "dna", f"{spec.label} should be DNA"

    def test_channel_labels_in_order(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        loaded_labels = tuple(t.label for t in manifest.targets)
        assert loaded_labels == _EXPECTED_LABELS

    def test_control_channel_role(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        assert manifest.targets[-1].role == "control"
        assert manifest.targets[-1].label == "CTRL_16S"

    def test_detect_channel_roles(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        for spec in manifest.targets[:-1]:
            assert spec.role == "detect", f"{spec.label} should be 'detect'"

    def test_extraction_params(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        assert isinstance(manifest.extraction, ExtractionParams)
        assert manifest.extraction.sample_volume_ul == 1000.0
        assert manifest.extraction.extraction_efficiency == 0.50
        assert manifest.extraction.eluate_volume_ul == 100.0
        assert manifest.extraction.reaction_input_ul == 5.0

    def test_pool_params(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        assert isinstance(manifest.pool, PoolingParams)
        assert manifest.pool.stock_conc_um == 300.0
        assert manifest.pool.total_pool_volume_ul == 500.0

    def test_pool_stock_above_six_channel_minimum(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        assert manifest.pool.stock_conc_um >= _MIN_STOCK_UM

    def test_sets_json_paths_resolved(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths)

        manifest = load_cartridge_manifest(manifest_path)

        for spec in manifest.targets:
            assert spec.sets_json.exists(), f"{spec.label}: sets_json not found"

    def test_lamp_chemistry_rejects_rna_channel(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        # Inject an RNA channel into an otherwise all-DNA lamp manifest
        roles = ["detect"] * 5 + ["control"]
        targets = [
            {
                "label": label,
                "sets_json": str(p),
                "target_na": "rna" if label == "SRB_dsrB" else "dna",
                "role": role,
            }
            for label, p, role in zip(_EXPECTED_LABELS, target_paths, roles, strict=True)
        ]
        payload: dict[str, Any] = {
            "cartridge_id": "bad_cart",
            "chemistry": "lamp",
            "max_channels": 6,
            "extraction": {
                "sample_volume_ul": 1000.0,
                "extraction_efficiency": 0.50,
                "eluate_volume_ul": 100.0,
                "reaction_input_ul": 5.0,
            },
            "device_window_min": 60.0,
            "pool": {"stock_conc_um": 300.0, "total_pool_volume_ul": 500.0},
            "vendor": "idt",
            "targets": targets,
        }
        manifest_path.write_text(yaml.safe_dump(payload), encoding="utf-8")
        with pytest.raises(ValueError, match="chemistry"):
            load_cartridge_manifest(manifest_path)

    def test_channel_cap_enforced(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        _write_mic_manifest_yaml(manifest_path, target_paths, max_channels=4)
        with pytest.raises(ValueError, match="max_channels"):
            load_cartridge_manifest(manifest_path)

    def test_duplicate_label_rejected(self, tmp_path: Path) -> None:
        target_paths = _make_six_channel_fixture(tmp_path)
        manifest_path = tmp_path / "mic_cart.yaml"
        roles = ["detect"] * 5 + ["control"]
        targets = [
            {
                "label": "SRB_dsrB",  # duplicate label
                "sets_json": str(p),
                "target_na": "dna",
                "role": role,
            }
            for p, role in zip(target_paths, roles, strict=True)
        ]
        payload: dict[str, Any] = {
            "cartridge_id": "dup_cart",
            "chemistry": "lamp",
            "max_channels": 6,
            "extraction": {
                "sample_volume_ul": 1000.0,
                "extraction_efficiency": 0.50,
                "eluate_volume_ul": 100.0,
                "reaction_input_ul": 5.0,
            },
            "device_window_min": 60.0,
            "pool": {"stock_conc_um": 300.0, "total_pool_volume_ul": 500.0},
            "vendor": "idt",
            "targets": targets,
        }
        manifest_path.write_text(yaml.safe_dump(payload), encoding="utf-8")
        with pytest.raises(ValueError, match="duplicate"):
            load_cartridge_manifest(manifest_path)
