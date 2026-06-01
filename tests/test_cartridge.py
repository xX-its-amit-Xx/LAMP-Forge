"""Tests for the one-shot BioID cartridge build orchestrator.

The orchestrator only composes existing LAMP-Forge primitives, so these tests
inject a synthetic primer_sets.json fixture (no primer3, no blast, no mafft)
and exercise the manifest schema, the GO/WARNING/NO_GO aggregation, and the
end-to-end run path including manifest-hash determinism.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from lamp_forge.cartridge import (
    CARTRIDGE_SCHEMA_VERSION,
    CartridgeBuildResult,
    CartridgeManifest,
    CartridgeTargetSpec,
    _aggregate_status,
    load_cartridge_manifest,
    render_cartridge_html,
    run_cartridge,
    write_cartridge_manifest_json,
)
from lamp_forge.lod import ExtractionParams
from lamp_forge.pool import PoolingParams

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _primer_dict(role: str, seq: str, tm: float = 64.0) -> dict[str, Any]:
    """Minimal serialised Primer matching report.write_json output."""
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


def _write_sets_json(
    path: Path,
    set_id: str,
    *,
    seqs: dict[str, str],
    tm: float = 64.0,
    score: float = 0.9,
) -> None:
    """Write a small primer_sets.json with one set."""
    entry: dict[str, Any] = {
        "set_id": set_id,
        "region_id": "R000",
        "composite_score": score,
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


def _make_two_target_fixture(
    tmp_path: Path,
    *,
    rna_tm: float = 64.0,
    rna_label: str = "RNA1",
    dna_label: str = "DNA1",
) -> tuple[Path, Path]:
    """Write two non-clashing primer_sets.json files for a 2-channel cartridge."""
    rna = tmp_path / "rna_sets.json"
    dna = tmp_path / "dna_sets.json"
    _write_sets_json(
        rna,
        f"{rna_label}_set_001",
        seqs={
            "f3": "GACCTGACATTCGGATCACG",
            "b3": "TTGCAGTCGAACTGGTACGA",
            "fip": "CAGTTCGACTGCAAGGTACCGATCGAATGTCAGGTC",
            "bip": "GTACCAGTTCGACTGCAATGCGTGATCCGAATGTCA",
        },
        tm=rna_tm,
    )
    _write_sets_json(
        dna,
        f"{dna_label}_set_001",
        seqs={
            "f3": "AAACCCGGGTTTACGTAACG",
            "b3": "CCTTGGAATTGCCAATGCGT",
            "fip": "TTAGGCCAATTCCGGAAACCCGGGAATTAGGCCAAT",
            "bip": "GGGTTAACCGGAATTAGGCCAATTAGGGTTAACCGG",
        },
        tm=64.0,
    )
    return rna, dna


def _write_manifest_yaml(
    path: Path,
    *,
    cartridge_id: str = "test_cart",
    chemistry: str = "rt-lamp",
    max_channels: int = 4,
    sample_volume_ul: float = 1000.0,
    extraction_efficiency: float = 0.50,
    eluate_volume_ul: float = 100.0,
    reaction_input_ul: float = 5.0,
    device_window_min: float = 60.0,
    stock_conc_um: float = 100.0,
    total_pool_volume_ul: float = 500.0,
    vendor: str = "idt",
    targets: list[dict[str, Any]] | None = None,
    extra: dict[str, Any] | None = None,
) -> None:
    """Write a cartridge manifest YAML for tests."""
    payload: dict[str, Any] = {
        "cartridge_id": cartridge_id,
        "chemistry": chemistry,
        "max_channels": max_channels,
        "extraction": {
            "sample_volume_ul": sample_volume_ul,
            "extraction_efficiency": extraction_efficiency,
            "eluate_volume_ul": eluate_volume_ul,
            "reaction_input_ul": reaction_input_ul,
        },
        "device_window_min": device_window_min,
        "pool": {
            "stock_conc_um": stock_conc_um,
            "total_pool_volume_ul": total_pool_volume_ul,
        },
        "vendor": vendor,
        "targets": targets if targets is not None else [],
    }
    if extra:
        payload.update(extra)
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------


class TestManifestLoad:
    def test_loads_minimal_two_channel_yaml(self, tmp_path: Path) -> None:
        rna, dna = _make_two_target_fixture(tmp_path)
        manifest_path = tmp_path / "cart.yaml"
        _write_manifest_yaml(
            manifest_path,
            targets=[
                {"label": "RNA1", "sets_json": str(rna), "target_na": "rna"},
                {"label": "DNA1", "sets_json": str(dna), "target_na": "dna"},
            ],
        )
        manifest = load_cartridge_manifest(manifest_path)
        assert isinstance(manifest, CartridgeManifest)
        assert manifest.cartridge_id == "test_cart"
        assert manifest.chemistry == "rt-lamp"
        assert len(manifest.targets) == 2
        assert manifest.targets[0].target_na == "rna"
        assert manifest.targets[1].target_na == "dna"
        assert isinstance(manifest.extraction, ExtractionParams)
        assert isinstance(manifest.pool, PoolingParams)

    def test_channel_cap_enforced(self, tmp_path: Path) -> None:
        rna, dna = _make_two_target_fixture(tmp_path)
        manifest_path = tmp_path / "cart.yaml"
        _write_manifest_yaml(
            manifest_path,
            max_channels=1,  # < 2 targets => violation
            targets=[
                {"label": "RNA1", "sets_json": str(rna), "target_na": "rna"},
                {"label": "DNA1", "sets_json": str(dna), "target_na": "dna"},
            ],
        )
        with pytest.raises(ValueError, match="max_channels"):
            load_cartridge_manifest(manifest_path)

    def test_rna_dna_mix_is_supported(self, tmp_path: Path) -> None:
        rna, dna = _make_two_target_fixture(tmp_path)
        manifest_path = tmp_path / "cart.yaml"
        _write_manifest_yaml(
            manifest_path,
            chemistry="rt-lamp",
            targets=[
                {"label": "RNA1", "sets_json": str(rna), "target_na": "rna"},
                {"label": "DNA1", "sets_json": str(dna), "target_na": "dna"},
            ],
        )
        manifest = load_cartridge_manifest(manifest_path)
        nas = {t.target_na for t in manifest.targets}
        assert nas == {"rna", "dna"}

    def test_lamp_chemistry_rejects_rna_target(self, tmp_path: Path) -> None:
        rna, dna = _make_two_target_fixture(tmp_path)
        manifest_path = tmp_path / "cart.yaml"
        _write_manifest_yaml(
            manifest_path,
            chemistry="lamp",  # inconsistent with an rna target
            targets=[
                {"label": "RNA1", "sets_json": str(rna), "target_na": "rna"},
                {"label": "DNA1", "sets_json": str(dna), "target_na": "dna"},
            ],
        )
        with pytest.raises(ValueError, match="rt-lamp"):
            load_cartridge_manifest(manifest_path)

    def test_missing_top_level_field_raises(self, tmp_path: Path) -> None:
        manifest_path = tmp_path / "cart.yaml"
        # Missing 'targets' / 'extraction' / etc.
        manifest_path.write_text(
            yaml.safe_dump({"cartridge_id": "x", "chemistry": "lamp"}),
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="missing required field"):
            load_cartridge_manifest(manifest_path)

    def test_schema_version_constant_present(self) -> None:
        assert isinstance(CARTRIDGE_SCHEMA_VERSION, str)
        assert CARTRIDGE_SCHEMA_VERSION

    def test_duplicate_label_rejected(self, tmp_path: Path) -> None:
        rna, dna = _make_two_target_fixture(tmp_path)
        manifest_path = tmp_path / "cart.yaml"
        _write_manifest_yaml(
            manifest_path,
            targets=[
                {"label": "X", "sets_json": str(rna), "target_na": "rna"},
                {"label": "X", "sets_json": str(dna), "target_na": "dna"},
            ],
        )
        with pytest.raises(ValueError, match="duplicate"):
            load_cartridge_manifest(manifest_path)


# ---------------------------------------------------------------------------
# Aggregation logic
# ---------------------------------------------------------------------------


def _pt(
    label: str,
    *,
    preorder_status: str = "GO",
    preorder_reasons: tuple[str, ...] = (),
) -> dict[str, Any]:
    return {
        "label": label,
        "preorder_status": preorder_status,
        "preorder_reasons": list(preorder_reasons),
    }


_CLEAN_PANEL: dict[str, Any] = {
    "worst_inter_assay_dg": -3.0,
    "is_clean": True,
    "flagged_cross_dimers": [],
}
_DIRTY_PANEL: dict[str, Any] = {
    "worst_inter_assay_dg": -9.5,
    "is_clean": False,
    "flagged_cross_dimers": [{"pair": "A/x x B/y", "dg_kcal_mol": -9.5}],
}


class TestAggregateStatus:
    def test_all_go_returns_go(self) -> None:
        status, reasons = _aggregate_status([_pt("A"), _pt("B")], _CLEAN_PANEL, rt_results=[])
        assert status == "GO"
        assert any("pass preorder" in r for r in reasons)

    def test_any_target_no_go_returns_no_go(self) -> None:
        status, reasons = _aggregate_status(
            [_pt("A"), _pt("B", preorder_status="NO_GO", preorder_reasons=("TTP > window",))],
            _CLEAN_PANEL,
            rt_results=[],
        )
        assert status == "NO_GO"
        assert any("'B'" in r for r in reasons)

    def test_dirty_panel_forces_no_go(self) -> None:
        status, reasons = _aggregate_status([_pt("A"), _pt("B")], _DIRTY_PANEL, rt_results=[])
        assert status == "NO_GO"
        assert any("not clean" in r for r in reasons)

    def test_rt_zero_optimised_forces_no_go(self) -> None:
        status, reasons = _aggregate_status(
            [_pt("R")],
            _CLEAN_PANEL,
            rt_results=[{"label": "R", "n_sets": 5, "n_sets_rt_ok": 0}],
        )
        assert status == "NO_GO"
        assert any("RT-LAMP Tm check" in r for r in reasons)

    def test_warning_when_some_rt_sets_pass(self) -> None:
        status, reasons = _aggregate_status(
            [_pt("R")],
            _CLEAN_PANEL,
            rt_results=[{"label": "R", "n_sets": 5, "n_sets_rt_ok": 3}],
        )
        assert status == "WARNING"
        assert any("sub-optimal" in r for r in reasons)

    def test_warning_propagates_from_per_target(self) -> None:
        status, reasons = _aggregate_status(
            [_pt("A"), _pt("B", preorder_status="WARNING", preorder_reasons=("near limit",))],
            _CLEAN_PANEL,
            rt_results=[],
        )
        assert status == "WARNING"
        assert any("near limit" in r for r in reasons)

    def test_dna_only_skips_rt_results(self) -> None:
        # rt_results empty -> RT branch is a no-op.
        status, _ = _aggregate_status([_pt("A"), _pt("B")], _CLEAN_PANEL, rt_results=[])
        assert status == "GO"


# ---------------------------------------------------------------------------
# End-to-end run_cartridge against a 2-target fixture
# ---------------------------------------------------------------------------


_EXPECTED_ARTIFACTS: tuple[str, ...] = (
    "cartridge_manifest.json",
    "panel.json",
    "panel_primers.csv",
    "pool_plan.csv",
    "lod_panel.csv",
    "rt_check.csv",
    "vendor_order.csv",
    "build_report.html",
)


def _build_manifest(tmp_path: Path) -> CartridgeManifest:
    """In-memory CartridgeManifest pointing at two synthetic sets."""
    rna, dna = _make_two_target_fixture(tmp_path)
    return CartridgeManifest(
        cartridge_id="cart_test",
        chemistry="rt-lamp",
        max_channels=4,
        extraction=ExtractionParams(
            sample_volume_ul=1000.0,
            extraction_efficiency=0.5,
            eluate_volume_ul=100.0,
            reaction_input_ul=5.0,
        ),
        device_window_min=60.0,
        pool=PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0),
        vendor="idt",
        targets=(
            CartridgeTargetSpec(label="RNA1", sets_json=rna, target_na="rna", role="detect"),
            CartridgeTargetSpec(label="DNA1", sets_json=dna, target_na="dna", role="detect"),
        ),
    )


class TestRunCartridge:
    def test_writes_every_artifact(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        out_dir = tmp_path / "build"
        result = run_cartridge(manifest, out_dir, write_html=True)
        assert isinstance(result, CartridgeBuildResult)
        for name in _EXPECTED_ARTIFACTS:
            assert (out_dir / name).exists(), f"missing artifact: {name}"
        # Per-target preorder CSVs.
        assert (out_dir / "preorder_RNA1.csv").exists()
        assert (out_dir / "preorder_DNA1.csv").exists()

    def test_manifest_hash_is_deterministic_across_runs(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        r1 = run_cartridge(manifest, tmp_path / "b1", write_html=False)
        r2 = run_cartridge(manifest, tmp_path / "b2", write_html=False)
        assert r1.manifest_hash == r2.manifest_hash

    def test_manifest_hash_changes_when_primer_sequence_changes(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        r1 = run_cartridge(manifest, tmp_path / "b1", write_html=False)

        # Mutate one of the underlying primer_sets.json files (change one base).
        rna_path = manifest.targets[0].sets_json
        data = json.loads(rna_path.read_text(encoding="utf-8"))
        data["primer_sets"][0]["f3"]["sequence"] = (
            data["primer_sets"][0]["f3"]["sequence"][:-1] + "T"
        )
        rna_path.write_text(json.dumps(data), encoding="utf-8")

        r2 = run_cartridge(manifest, tmp_path / "b2", write_html=False)
        assert r1.manifest_hash != r2.manifest_hash

    def test_cartridge_manifest_json_has_expected_keys(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        out_dir = tmp_path / "build"
        run_cartridge(manifest, out_dir, write_html=False)
        payload = json.loads((out_dir / "cartridge_manifest.json").read_text())
        expected = {
            "schema_version",
            "generated_at",
            "lamp_forge_version",
            "cartridge_id",
            "chemistry",
            "max_channels",
            "n_targets",
            "device_window_min",
            "extraction",
            "pool_params",
            "vendor",
            "per_target",
            "panel_summary",
            "panel_lod_worst",
            "build_status",
            "reasons",
            "manifest_hash",
        }
        assert expected.issubset(set(payload.keys()))


# ---------------------------------------------------------------------------
# Panel-LoD worst row
# ---------------------------------------------------------------------------


class TestPanelLodWorstRow:
    def test_worst_row_is_max_across_channels(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        out_dir = tmp_path / "build"
        result = run_cartridge(manifest, out_dir, write_html=False)
        # Both channels share the same extraction chain, so LoDs are equal.
        # The worst-row label must match one of the channel labels.
        worst_label = result.panel_lod_worst.get("label", "")
        assert worst_label in {"RNA1", "DNA1"}

        # The lod_panel.csv should also carry a row prefixed "WORST(".
        csv_text = (out_dir / "lod_panel.csv").read_text(encoding="utf-8")
        assert "WORST(" in csv_text


# ---------------------------------------------------------------------------
# Manifest JSON serialiser + HTML renderer
# ---------------------------------------------------------------------------


class TestSerialisers:
    def test_write_cartridge_manifest_json_creates_file(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        result = run_cartridge(manifest, tmp_path / "b", write_html=False)
        out = tmp_path / "out" / "cartridge_manifest.json"
        write_cartridge_manifest_json(result, out)
        assert out.exists()
        payload = json.loads(out.read_text())
        assert payload["cartridge_id"] == manifest.cartridge_id

    def test_render_cartridge_html_writes_file(self, tmp_path: Path) -> None:
        manifest = _build_manifest(tmp_path)
        result = run_cartridge(manifest, tmp_path / "b", write_html=False)
        out = tmp_path / "report.html"
        render_cartridge_html(result, out)
        text = out.read_text(encoding="utf-8")
        assert manifest.cartridge_id in text
        assert "Build status" in text
