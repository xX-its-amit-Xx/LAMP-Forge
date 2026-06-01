"""One-shot BioID cartridge build orchestrator.

A BioVind-style portable diagnostic platform ships per-deployment as a
multi-channel **cartridge**: an SKU that bundles a fixed set of LAMP/RT-LAMP
assays (channels) sharing a single sample chamber, a single extraction step,
and a single device run-time window.  Spinning up a new cartridge SKU is
not just designing one assay -- it is composing the entire LAMP-Forge
pipeline across N targets and reconciling the per-target pass/fail signals
into a single cartridge-level **GO / WARNING / NO_GO** verdict.

Today users chain seven commands by hand (``panel``, ``pool``, ``lod``,
``preorder``, ``rt-check``, ``vendor-export``) and reconcile the outputs in
a spreadsheet.  This module collapses that workflow into one call::

    from lamp_forge.cartridge import load_cartridge_manifest, run_cartridge
    manifest = load_cartridge_manifest(Path("config/cartridge_resp_18ch.yaml"))
    result = run_cartridge(manifest, Path("results/resp_18ch_v1"))
    if result.build_status == "NO_GO":
        sys.exit(1)  # block vendor order submission

Every primitive used here is already part of LAMP-Forge -- this module owns
no new science, only the orchestration and the cartridge-level aggregation
logic that decides when a panel can be ordered.

Cartridge-level status rollup:

* **GO** -- every target is preorder-GO, panel ``is_clean=True``, every RNA
  channel passes the RT-LAMP Tm check.
* **WARNING** -- panel is clean, but at least one target is preorder-WARNING
  (e.g. TTP near the device-window limit, RT-LAMP sub-optimal on some sets).
* **NO_GO** -- any per-target NO_GO, panel ``is_clean=False`` (flagged
  cross-dimers), or RNA channel with zero RT-optimised sets.

Each :class:`CartridgeBuildResult` is signed with a SHA-256 ``manifest_hash``
over the canonical JSON of the manifest plus all per-target outputs, so the
``cartridge_manifest.json`` artifact is reproducible and tamper-evident.

References:
    Notomi T et al. (2000) Nucleic Acids Res 28(12):e63.
    Tanner NA et al. (2012) BioTechniques 53(2):81-89.
    NEB WarmStart(R) LAMP Kit manual (E1700).
"""

from __future__ import annotations

import base64
import hashlib
import io
import json
import logging
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal

import yaml

from lamp_forge.lod import ExtractionParams, estimate_lod
from lamp_forge.panel import PanelResult, run_panel
from lamp_forge.pool import PoolingParams, build_pool_plan, write_pool_csv
from lamp_forge.preorder import run_preorder
from lamp_forge.rt_lamp import (
    RtLampCheckResult,
    RtLampParams,
    TargetNucleicAcid,
    check_primer_sets_for_rt_lamp,
    write_rt_check_csv,
)
from lamp_forge.ttp import TtpParams, TtpPreset
from lamp_forge.vendor_export import VendorFormat, rows_from_panel_json, write_vendor_csv

logger = logging.getLogger(__name__)

# Manifest schema version. Bump when the on-disk JSON shape changes.
CARTRIDGE_SCHEMA_VERSION: str = "1.0"

# Chemistry enum (string-literal): "lamp" for DNA-only cartridges,
# "rt-lamp" for cartridges with at least one RNA channel.
Chemistry = Literal["lamp", "rt-lamp"]
TargetNa = Literal["dna", "rna"]
TargetRole = Literal["detect", "control"]
Vendor = Literal["idt", "twist"]
BuildStatus = Literal["GO", "WARNING", "NO_GO"]


# ---------------------------------------------------------------------------
# Manifest dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CartridgeTargetSpec:
    """One channel on the cartridge.

    Attributes:
        label: Short channel label (e.g. ``"SARS2_N"``).  Used in artifact
            filenames and the per-target JSON block.
        sets_json: Path to the ``primer_sets.json`` produced for this target
            by a prior ``lamp-forge run``.
        target_na: ``"dna"`` or ``"rna"``.  RNA targets are routed through
            the RT-LAMP Tm check and the rt-lamp TTP preset.
        role: ``"detect"`` for diagnostic channels, ``"control"`` for the
            16S/RNase-P sample-adequacy control channel.  The cartridge-level
            status rollup is unaffected by the role -- the build either fits
            on the cartridge or it does not -- but the role is preserved in
            the manifest so downstream tooling can render the channels in
            order.
    """

    label: str
    sets_json: Path
    target_na: TargetNa
    role: TargetRole = "detect"


@dataclass(frozen=True, slots=True)
class CartridgeManifest:
    """A vetted multi-target BioID cartridge build specification.

    Attributes:
        cartridge_id: Short cartridge SKU identifier (e.g. ``"resp_18ch_v1"``).
        chemistry: ``"lamp"`` if every channel is DNA, ``"rt-lamp"`` if at
            least one channel is RNA.  Validated against the per-channel
            ``target_na`` at load time.
        max_channels: Physical channel cap of the cartridge platform
            (e.g. 18 for a BioVind BioID-18).  ``len(targets)`` must not
            exceed this.
        extraction: Sample-preparation chain shared by every channel.
        device_window_min: Maximum reaction time allowed on the device
            (e.g. 60 min).  Applied uniformly to every channel.
        pool: Pooling parameters (stock concentration, total pool volume)
            shared by every channel.
        vendor: ``"idt"`` or ``"twist"`` -- the vendor order CSV format.
        targets: Ordered tuple of channel specs.
    """

    cartridge_id: str
    chemistry: Chemistry
    max_channels: int
    extraction: ExtractionParams
    device_window_min: float
    pool: PoolingParams
    vendor: Vendor
    targets: tuple[CartridgeTargetSpec, ...]


@dataclass(frozen=True, slots=True)
class CartridgeBuildResult:
    """End-to-end cartridge build record.

    Attributes:
        manifest: The :class:`CartridgeManifest` this result was built from.
        per_target: Tuple of per-channel summary dicts (see
            :func:`_per_target_summary` for the field set).
        panel_summary: Subset of the ``panel.json`` payload (worst
            inter-assay dG, is_clean flag, flagged cross-dimers).
        panel_lod_worst: Worst-case LoD across channels (max
            ``lod_copies_per_ml`` per the rollup convention -- "what fraction
            of organisms in the sample will the *least* sensitive channel
            miss?").
        build_status: Cartridge-level ``"GO"`` / ``"WARNING"`` / ``"NO_GO"``.
        reasons: Human-readable strings explaining the rollup.
        schema_version: Version of the on-disk cartridge_manifest.json shape.
        generated_at: ISO-8601 UTC timestamp.
        manifest_hash: SHA-256 over the canonical JSON of manifest + outputs.
    """

    manifest: CartridgeManifest
    per_target: tuple[dict[str, Any], ...]
    panel_summary: dict[str, Any]
    panel_lod_worst: dict[str, Any]
    build_status: BuildStatus
    reasons: tuple[str, ...]
    schema_version: str
    generated_at: str
    manifest_hash: str


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------


_REQUIRED_KEYS: tuple[str, ...] = (
    "cartridge_id",
    "chemistry",
    "max_channels",
    "extraction",
    "device_window_min",
    "pool",
    "vendor",
    "targets",
)

_REQUIRED_EXTRACTION_KEYS: tuple[str, ...] = (
    "sample_volume_ul",
    "extraction_efficiency",
    "eluate_volume_ul",
    "reaction_input_ul",
)

_REQUIRED_POOL_KEYS: tuple[str, ...] = ("stock_conc_um", "total_pool_volume_ul")

_REQUIRED_TARGET_KEYS: tuple[str, ...] = ("label", "sets_json", "target_na")


def load_cartridge_manifest(path: Path) -> CartridgeManifest:
    """Load and validate a cartridge manifest YAML file.

    The YAML schema mirrors the :class:`CartridgeManifest` dataclass.  All
    required keys must be present; per-target ``sets_json`` paths are
    resolved relative to the manifest file's directory if not absolute.

    Args:
        path: Path to the manifest YAML.

    Returns:
        Validated :class:`CartridgeManifest`.

    Raises:
        ValueError: If any required key is missing, the chemistry/target-na
            combination is inconsistent, the number of targets exceeds
            ``max_channels``, or a per-target ``sets_json`` is missing.
    """
    path = Path(path)
    with path.open(encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)

    if not isinstance(raw, dict):
        raise ValueError(f"{path}: top-level YAML must be a mapping")

    missing = [k for k in _REQUIRED_KEYS if k not in raw]
    if missing:
        raise ValueError(f"{path}: missing required field(s): {', '.join(missing)}")

    chemistry = str(raw["chemistry"]).lower()
    if chemistry not in ("lamp", "rt-lamp"):
        raise ValueError(f"{path}: chemistry must be 'lamp' or 'rt-lamp', got {chemistry!r}")

    max_channels = int(raw["max_channels"])
    if max_channels < 2:
        raise ValueError(f"{path}: max_channels must be >= 2, got {max_channels}")

    device_window_min = float(raw["device_window_min"])
    if device_window_min <= 0.0:
        raise ValueError(f"{path}: device_window_min must be > 0, got {device_window_min}")

    extraction_raw = raw["extraction"]
    if not isinstance(extraction_raw, dict):
        raise ValueError(f"{path}: 'extraction' must be a mapping")
    missing_ext = [k for k in _REQUIRED_EXTRACTION_KEYS if k not in extraction_raw]
    if missing_ext:
        raise ValueError(f"{path}: extraction missing field(s): {', '.join(missing_ext)}")
    extraction = ExtractionParams(
        sample_volume_ul=float(extraction_raw["sample_volume_ul"]),
        extraction_efficiency=float(extraction_raw["extraction_efficiency"]),
        eluate_volume_ul=float(extraction_raw["eluate_volume_ul"]),
        reaction_input_ul=float(extraction_raw["reaction_input_ul"]),
    )

    pool_raw = raw["pool"]
    if not isinstance(pool_raw, dict):
        raise ValueError(f"{path}: 'pool' must be a mapping")
    missing_pool = [k for k in _REQUIRED_POOL_KEYS if k not in pool_raw]
    if missing_pool:
        raise ValueError(f"{path}: pool missing field(s): {', '.join(missing_pool)}")
    pool = PoolingParams(
        stock_conc_um=float(pool_raw["stock_conc_um"]),
        total_pool_volume_ul=float(pool_raw["total_pool_volume_ul"]),
    )

    vendor = str(raw["vendor"]).lower()
    if vendor not in ("idt", "twist"):
        raise ValueError(f"{path}: vendor must be 'idt' or 'twist', got {vendor!r}")

    targets_raw = raw["targets"]
    if not isinstance(targets_raw, list) or not targets_raw:
        raise ValueError(f"{path}: 'targets' must be a non-empty list")

    targets: list[CartridgeTargetSpec] = []
    seen_labels: set[str] = set()
    base_dir = path.parent
    for i, t_raw in enumerate(targets_raw):
        if not isinstance(t_raw, dict):
            raise ValueError(f"{path}: targets[{i}] must be a mapping")
        missing_t = [k for k in _REQUIRED_TARGET_KEYS if k not in t_raw]
        if missing_t:
            raise ValueError(f"{path}: targets[{i}] missing field(s): {', '.join(missing_t)}")
        label = str(t_raw["label"])
        if label in seen_labels:
            raise ValueError(f"{path}: duplicate target label {label!r}")
        seen_labels.add(label)

        sets_json = Path(str(t_raw["sets_json"]))
        if not sets_json.is_absolute():
            sets_json = (base_dir / sets_json).resolve()
        if not sets_json.exists():
            raise ValueError(f"{path}: targets[{i}] sets_json not found: {sets_json}")

        target_na = str(t_raw["target_na"]).lower()
        if target_na not in ("dna", "rna"):
            raise ValueError(
                f"{path}: targets[{i}] target_na must be 'dna' or 'rna', got {target_na!r}"
            )

        role = str(t_raw.get("role", "detect")).lower()
        if role not in ("detect", "control"):
            raise ValueError(
                f"{path}: targets[{i}] role must be 'detect' or 'control', got {role!r}"
            )

        targets.append(
            CartridgeTargetSpec(
                label=label,
                sets_json=sets_json,
                target_na=target_na,  # type: ignore[arg-type]
                role=role,  # type: ignore[arg-type]
            )
        )

    if len(targets) > max_channels:
        raise ValueError(
            f"{path}: n_targets ({len(targets)}) exceeds max_channels ({max_channels})"
        )

    has_rna = any(t.target_na == "rna" for t in targets)
    if chemistry == "lamp" and has_rna:
        raise ValueError(
            f"{path}: chemistry='lamp' is inconsistent with at least one RNA target; "
            "use chemistry='rt-lamp'."
        )

    return CartridgeManifest(
        cartridge_id=str(raw["cartridge_id"]),
        chemistry=chemistry,  # type: ignore[arg-type]
        max_channels=max_channels,
        extraction=extraction,
        device_window_min=device_window_min,
        pool=pool,
        vendor=vendor,  # type: ignore[arg-type]
        targets=tuple(targets),
    )


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def _aggregate_status(
    per_target: Sequence[dict[str, Any]],
    panel_summary: dict[str, Any],
    rt_results: Sequence[dict[str, Any]],
) -> tuple[BuildStatus, tuple[str, ...]]:
    """Roll per-target + panel + RT outputs into a cartridge-level verdict.

    Rules:

    * Any per-target ``preorder_status == "NO_GO"`` -> cartridge ``NO_GO``.
    * Panel ``is_clean == False`` -> cartridge ``NO_GO``.
    * Any RNA channel with ``is_rt_optimized == False`` for *every* set ->
      cartridge ``NO_GO`` (zero usable sets for that channel).
    * Any per-target ``preorder_status == "WARNING"`` (and no NO_GO above)
      -> cartridge ``WARNING``.
    * Any RNA channel with at least one set RT-optimised but not all sets ->
      cartridge ``WARNING``.
    * Otherwise -> cartridge ``GO``.

    Args:
        per_target: Per-channel summary dicts (must include
            ``label`` and ``preorder_status``).
        panel_summary: Output of :func:`_panel_summary_from_result`.
        rt_results: One stacked entry per RNA channel with
            ``label``, ``n_sets``, ``n_sets_rt_ok``.

    Returns:
        ``(status, reasons)`` tuple.
    """
    reasons: list[str] = []
    no_go = False

    for pt in per_target:
        if pt.get("preorder_status") == "NO_GO":
            no_go = True
            label = pt.get("label", "?")
            reasons.append(
                f"Channel {label!r}: per-target preorder NO_GO -- "
                f"{'; '.join(pt.get('preorder_reasons', ())) or 'no reason given'}"
            )

    if not panel_summary.get("is_clean", True):
        no_go = True
        n_flagged = len(panel_summary.get("flagged_cross_dimers", ()))
        worst = panel_summary.get("worst_inter_assay_dg", 0.0)
        reasons.append(
            f"Panel is not clean: {n_flagged} inter-assay heterodimer(s) flagged "
            f"(worst dG {worst:.2f} kcal/mol)."
        )

    for rt in rt_results:
        n_sets = int(rt.get("n_sets", 0))
        n_ok = int(rt.get("n_sets_rt_ok", 0))
        if n_sets > 0 and n_ok == 0:
            no_go = True
            reasons.append(
                f"Channel {rt.get('label', '?')!r}: 0/{n_sets} RNA sets pass the "
                "RT-LAMP Tm check -- no usable set for one-step RT-LAMP."
            )

    if no_go:
        return "NO_GO", tuple(reasons)

    warning = False
    for pt in per_target:
        if pt.get("preorder_status") == "WARNING":
            warning = True
            label = pt.get("label", "?")
            reasons.append(
                f"Channel {label!r}: per-target preorder WARNING -- "
                f"{'; '.join(pt.get('preorder_reasons', ())) or 'review before ordering'}"
            )

    for rt in rt_results:
        n_sets = int(rt.get("n_sets", 0))
        n_ok = int(rt.get("n_sets_rt_ok", 0))
        if n_sets > 0 and 0 < n_ok < n_sets:
            warning = True
            reasons.append(
                f"Channel {rt.get('label', '?')!r}: only {n_ok}/{n_sets} RNA sets "
                "pass the RT-LAMP Tm check; some sets are sub-optimal."
            )

    if warning:
        return "WARNING", tuple(reasons)

    reasons.append(
        f"All {len(per_target)} channel(s) pass preorder, panel is clean, "
        "and every RNA channel is RT-LAMP optimised."
    )
    return "GO", tuple(reasons)


# ---------------------------------------------------------------------------
# Helpers used inside run_cartridge
# ---------------------------------------------------------------------------


def _panel_summary_from_result(result: PanelResult) -> dict[str, Any]:
    """Compact dict summarising a :class:`PanelResult` for the manifest."""
    return {
        "worst_inter_assay_dg": result.worst_dg,
        "is_clean": result.is_clean,
        "search_mode": result.search_mode,
        "n_combinations_evaluated": result.n_combinations_evaluated,
        "dimer_dg_threshold": result.dimer_dg_threshold,
        "flagged_cross_dimers": [
            {
                "pair": cd.label,
                "dg_kcal_mol": cd.dg,
            }
            for cd in result.flagged
        ],
    }


def _composite_score_for(panel_result: PanelResult, label: str) -> float:
    """Look up the chosen set's composite score for one channel."""
    chosen = panel_result.selection.get(label)
    if chosen is None:
        return 0.0
    return float(chosen.composite_score)


def _chosen_set_id_for(panel_result: PanelResult, label: str) -> str:
    chosen = panel_result.selection.get(label)
    if chosen is None:
        return ""
    return str(chosen.set_id)


def _canonical_json(obj: Any) -> bytes:
    """Deterministic JSON encoding for hashing -- sorted keys, no whitespace."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")


def _manifest_canonical_obj(manifest: CartridgeManifest) -> dict[str, Any]:
    """Hashable canonical dict of the input manifest (paths + params + sequences)."""
    target_blocks: list[dict[str, Any]] = []
    for t in manifest.targets:
        # Include the file's content hash too, so changing a primer sequence
        # changes the manifest_hash even if the path string stays the same.
        with t.sets_json.open("rb") as fh:
            content_sha = hashlib.sha256(fh.read()).hexdigest()
        target_blocks.append(
            {
                "label": t.label,
                "sets_json": str(t.sets_json),
                "sets_json_sha256": content_sha,
                "target_na": t.target_na,
                "role": t.role,
            }
        )
    return {
        "cartridge_id": manifest.cartridge_id,
        "chemistry": manifest.chemistry,
        "max_channels": manifest.max_channels,
        "device_window_min": manifest.device_window_min,
        "extraction": {
            "sample_volume_ul": manifest.extraction.sample_volume_ul,
            "extraction_efficiency": manifest.extraction.extraction_efficiency,
            "eluate_volume_ul": manifest.extraction.eluate_volume_ul,
            "reaction_input_ul": manifest.extraction.reaction_input_ul,
        },
        "pool": {
            "stock_conc_um": manifest.pool.stock_conc_um,
            "total_pool_volume_ul": manifest.pool.total_pool_volume_ul,
        },
        "vendor": manifest.vendor,
        "targets": target_blocks,
    }


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------


def run_cartridge(
    manifest: CartridgeManifest,
    out_dir: Path,
    *,
    write_html: bool = True,
) -> CartridgeBuildResult:
    """Run every LAMP-Forge primitive a cartridge build needs, in one call.

    Composes :func:`panel.run_panel`, :func:`pool.build_pool_plan`,
    :func:`preorder.run_preorder`, :func:`rt_lamp.check_primer_sets_for_rt_lamp`,
    :func:`lod.estimate_lod`, and :func:`vendor_export.rows_from_panel_json`
    into one signed cartridge build.

    Args:
        manifest: Loaded cartridge manifest.
        out_dir: Output directory; created if missing.
        write_html: If True, render the cartridge HTML report and the
            ``panel_report.html`` passthrough.  Disable in CI to avoid
            matplotlib usage.

    Returns:
        :class:`CartridgeBuildResult` carrying every output path and the
        cartridge-level GO/WARNING/NO_GO verdict.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Panel selection (one set per channel that minimises cross-dimer)
    # ------------------------------------------------------------------
    set_specs = {t.label: t.sets_json for t in manifest.targets}
    panel_result = run_panel(
        set_specs,
        out_dir,
        top_per_target=5,
        dimer_dg_threshold=-5.0,
        generate_html=write_html,
    )
    panel_summary = _panel_summary_from_result(panel_result)

    # ------------------------------------------------------------------
    # 2. Pool plan over the selected panel
    # ------------------------------------------------------------------
    with (out_dir / "panel.json").open(encoding="utf-8") as fh:
        panel_json_payload = json.load(fh)
    pool_plan = build_pool_plan(panel_json_payload["selection"], manifest.pool)
    write_pool_csv(pool_plan, out_dir / "pool_plan.csv")

    # ------------------------------------------------------------------
    # 3. Per-target preorder + LoD + RT-check (stacked) + vendor export
    # ------------------------------------------------------------------
    per_target: list[dict[str, Any]] = []
    rt_results_for_status: list[dict[str, Any]] = []
    all_rt_check_rows: list[RtLampCheckResult] = []
    lod_panel_rows: list[dict[str, Any]] = []

    rt_params = RtLampParams()

    for tgt in manifest.targets:
        with tgt.sets_json.open(encoding="utf-8") as fh:
            sets_payload = json.load(fh)
        sets_data = sets_payload.get("primer_sets", [])
        if not isinstance(sets_data, list) or not sets_data:
            raise ValueError(
                f"Channel {tgt.label!r}: primer_sets.json has no usable sets ({tgt.sets_json})"
            )

        target_na = TargetNucleicAcid(tgt.target_na)
        preset = TtpPreset.RT_LAMP if target_na is TargetNucleicAcid.RNA else TtpPreset.DNA_LAMP
        ttp_params = TtpParams.from_preset(preset, device_window_min=manifest.device_window_min)

        preorder_result = run_preorder(
            sets_data,
            manifest.extraction,
            ttp_params,
            target_na=target_na,
            rt_params=rt_params,
        )
        write_preorder_csv_path = out_dir / f"preorder_{tgt.label}.csv"
        # Re-use the standard CSV writer.
        from lamp_forge.preorder import write_preorder_csv

        write_preorder_csv(preorder_result, write_preorder_csv_path)

        # Per-channel LoD (also rolled into the panel-LoD summary below).
        lod_est = estimate_lod(manifest.extraction, detection_probability=0.95)

        # RT-check rows for the stacked rt_check.csv. DNA channels still get
        # an entry so the panel-wide file is complete; the helper returns a
        # uniform "RT not required" row for DNA.
        rt_check_results = check_primer_sets_for_rt_lamp(
            sets_data, target_na=target_na, params=rt_params
        )
        all_rt_check_rows.extend(rt_check_results)

        n_sets = len(rt_check_results)
        n_ok = sum(1 for r in rt_check_results if r.is_rt_optimized)
        if target_na is TargetNucleicAcid.RNA:
            rt_results_for_status.append(
                {"label": tgt.label, "n_sets": n_sets, "n_sets_rt_ok": n_ok}
            )

        target_block: dict[str, Any] = {
            "label": tgt.label,
            "role": tgt.role,
            "target_na": tgt.target_na,
            "set_id": _chosen_set_id_for(panel_result, tgt.label),
            "composite_score": _composite_score_for(panel_result, tgt.label),
            "preorder_status": str(preorder_result.status),
            "preorder_reasons": list(preorder_result.reasons),
            "rt_check": {
                "n_sets": n_sets,
                "n_sets_rt_ok": n_ok if target_na is TargetNucleicAcid.RNA else -1,
            },
            "lod": {
                "lod_copies_per_rxn": lod_est.lod_copies_per_reaction,
                "lod_copies_per_ml": lod_est.lod_copies_per_ml,
                "detection_probability": lod_est.detection_probability,
            },
            "ttp_at_lod_min": preorder_result.ttp_at_lod_min,
        }
        per_target.append(target_block)

        lod_panel_rows.append(
            {
                "label": tgt.label,
                "lod_copies_per_rxn": lod_est.lod_copies_per_reaction,
                "lod_copies_per_ml": lod_est.lod_copies_per_ml,
            }
        )

    # ------------------------------------------------------------------
    # 4. Stacked RT-check CSV across every channel
    # ------------------------------------------------------------------
    write_rt_check_csv(all_rt_check_rows, out_dir / "rt_check.csv")

    # ------------------------------------------------------------------
    # 5. Panel-LoD CSV (one row per channel + worst-row)
    # ------------------------------------------------------------------
    panel_lod_worst = _write_panel_lod_csv(lod_panel_rows, out_dir / "lod_panel.csv")

    # ------------------------------------------------------------------
    # 6. Single vendor order CSV across every channel
    # ------------------------------------------------------------------
    vendor_rows = rows_from_panel_json(panel_json_payload["selection"])
    vendor_fmt = VendorFormat(manifest.vendor)
    write_vendor_csv(vendor_rows, out_dir / "vendor_order.csv", vendor_fmt)

    # ------------------------------------------------------------------
    # 7. Cartridge-level status rollup
    # ------------------------------------------------------------------
    build_status, reasons = _aggregate_status(per_target, panel_summary, rt_results_for_status)

    generated_at = datetime.now(UTC).isoformat()

    # ------------------------------------------------------------------
    # 8. Manifest hash -- SHA-256 over canonical (manifest + outputs).
    #     Excludes generated_at and the hash field itself so two back-to-back
    #     runs over the same inputs produce the same hash.
    # ------------------------------------------------------------------
    hash_payload = {
        "manifest": _manifest_canonical_obj(manifest),
        "per_target": per_target,
        "panel_summary": panel_summary,
        "panel_lod_worst": panel_lod_worst,
        "build_status": build_status,
        "reasons": list(reasons),
        "schema_version": CARTRIDGE_SCHEMA_VERSION,
    }
    manifest_hash = hashlib.sha256(_canonical_json(hash_payload)).hexdigest()

    result = CartridgeBuildResult(
        manifest=manifest,
        per_target=tuple(per_target),
        panel_summary=panel_summary,
        panel_lod_worst=panel_lod_worst,
        build_status=build_status,
        reasons=reasons,
        schema_version=CARTRIDGE_SCHEMA_VERSION,
        generated_at=generated_at,
        manifest_hash=manifest_hash,
    )

    write_cartridge_manifest_json(result, out_dir / "cartridge_manifest.json")
    if write_html:
        render_cartridge_html(result, out_dir / "build_report.html")

    logger.info(
        "Cartridge %s built: status=%s (%d target(s), hash=%s)",
        manifest.cartridge_id,
        build_status,
        len(per_target),
        manifest_hash[:12],
    )
    return result


# ---------------------------------------------------------------------------
# Panel-LoD writer (small enough to live here -- not a new science module)
# ---------------------------------------------------------------------------


def _write_panel_lod_csv(rows: list[dict[str, Any]], path: Path) -> dict[str, Any]:
    """Write the per-channel LoD CSV plus an aggregated worst-row.

    The "worst" channel is the one with the highest ``lod_copies_per_ml``
    (i.e. the least sensitive channel -- the cartridge inherits the worst
    LoD across channels because a sample-positive must clear *every* channel
    that interrogates that pathogen).

    Args:
        rows: One ``{label, lod_copies_per_rxn, lod_copies_per_ml}`` per channel.
        path: Output CSV path.

    Returns:
        The worst-row dict, useful for the manifest JSON.
    """
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    worst = max(rows, key=lambda r: float(r["lod_copies_per_ml"])) if rows else {}
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["channel", "lod_copies_per_rxn", "lod_copies_per_ml"])
        for r in rows:
            writer.writerow(
                [
                    r["label"],
                    f"{float(r['lod_copies_per_rxn']):.4f}",
                    f"{float(r['lod_copies_per_ml']):.2f}",
                ]
            )
        if worst:
            writer.writerow(
                [
                    f"WORST({worst['label']})",
                    f"{float(worst['lod_copies_per_rxn']):.4f}",
                    f"{float(worst['lod_copies_per_ml']):.2f}",
                ]
            )
    return dict(worst)


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------


def write_cartridge_manifest_json(result: CartridgeBuildResult, path: Path) -> None:
    """Write the signed cartridge build record to JSON.

    The file is the top-level deliverable: the CI gate, the vendor order
    reviewer, and the wet-lab handoff all consume this single file.

    Args:
        result: :class:`CartridgeBuildResult` from :func:`run_cartridge`.
        path: Output JSON path; parent directories are created.
    """
    from lamp_forge import __version__

    manifest = result.manifest
    payload = {
        "schema_version": result.schema_version,
        "generated_at": result.generated_at,
        "lamp_forge_version": __version__,
        "cartridge_id": manifest.cartridge_id,
        "chemistry": manifest.chemistry,
        "max_channels": manifest.max_channels,
        "n_targets": len(manifest.targets),
        "device_window_min": manifest.device_window_min,
        "extraction": {
            "sample_volume_ul": manifest.extraction.sample_volume_ul,
            "extraction_efficiency": manifest.extraction.extraction_efficiency,
            "eluate_volume_ul": manifest.extraction.eluate_volume_ul,
            "reaction_input_ul": manifest.extraction.reaction_input_ul,
        },
        "pool_params": {
            "stock_conc_um": manifest.pool.stock_conc_um,
            "total_pool_volume_ul": manifest.pool.total_pool_volume_ul,
        },
        "vendor": manifest.vendor,
        "per_target": list(result.per_target),
        "panel_summary": result.panel_summary,
        "panel_lod_worst": result.panel_lod_worst,
        "build_status": result.build_status,
        "reasons": list(result.reasons),
        "manifest_hash": result.manifest_hash,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------


_CARTRIDGE_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>LAMP-Forge cartridge build &mdash; {{ cartridge_id }}</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #1f2937; line-height: 1.5; }
  h1 { border-bottom: 3px solid #0369a1; padding-bottom: 0.3em; }
  h2 { margin-top: 2em; color: #0369a1; }
  .verdict { padding: 0.8em 1em; border-radius: 6px; font-weight: bold; margin: 1em 0; }
  .go { background: #dcfce7; border-left: 4px solid #15803d; color: #14532d; }
  .warn { background: #fef9c3; border-left: 4px solid #ca8a04; color: #713f12; }
  .nogo { background: #fee2e2; border-left: 4px solid #b91c1c; color: #7f1d1d; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 0.88em; }
  th, td { border: 1px solid #cbd5e1; padding: 6px 10px; text-align: left; }
  th { background: #e0f2fe; color: #075985; }
  tr:nth-child(even) { background: #f8fafc; }
  .seq { font-family: ui-monospace, Consolas, monospace; font-size: 0.85em; }
  .plot { margin: 1em 0; text-align: center; }
  .plot img { max-width: 100%; border: 1px solid #e2e8f0; border-radius: 4px; }
  .footer { margin-top: 3em; padding-top: 1em; border-top: 1px solid #e2e8f0;
            font-size: 0.85em; color: #64748b; }
  code { background: #f1f5f9; padding: 0.05em 0.3em; border-radius: 3px; }
</style>
</head>
<body>
<h1>LAMP-Forge cartridge build &mdash; {{ cartridge_id }}</h1>

<div class="verdict {{ verdict_class }}">
  Build status: {{ build_status }}
  &middot; chemistry {{ chemistry }}
  &middot; {{ n_targets }}/{{ max_channels }} channel(s)
  &middot; device window {{ "%.0f"|format(device_window) }} min
</div>

<h2>Reasons</h2>
<ul>
  {% for r in reasons %}
  <li>{{ r }}</li>
  {% endfor %}
</ul>

<h2>Channel table</h2>
<table>
  <thead><tr>
    <th>Channel</th><th>Role</th><th>NA</th><th>Set</th><th>Score</th>
    <th>Preorder</th><th>LoD (cp/mL)</th><th>RT-LAMP OK</th>
  </tr></thead>
  <tbody>
  {% for t in per_target %}
    <tr>
      <td><strong>{{ t.label }}</strong></td>
      <td>{{ t.role }}</td>
      <td>{{ t.target_na }}</td>
      <td class="seq">{{ t.set_id }}</td>
      <td>{{ "%.3f"|format(t.composite_score) }}</td>
      <td>{{ t.preorder_status }}</td>
      <td>{{ "%.1f"|format(t.lod.lod_copies_per_ml) }}</td>
      <td>{{ t.rt_check.n_sets_rt_ok if t.target_na == "rna" else "N/A" }}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>

<h2>Panel cross-reactivity</h2>
<p>Worst inter-assay heterodimer dG:
   <code>{{ "%.2f"|format(panel_summary.worst_inter_assay_dg) }}</code> kcal/mol
   (clean = {{ panel_summary.is_clean }}).</p>
{% if panel_summary.flagged_cross_dimers %}
<table>
  <thead><tr><th>Pair</th><th>dG (kcal/mol)</th></tr></thead>
  <tbody>
  {% for cd in panel_summary.flagged_cross_dimers %}
    <tr><td>{{ cd.pair }}</td><td>{{ "%.2f"|format(cd.dg_kcal_mol) }}</td></tr>
  {% endfor %}
  </tbody>
</table>
{% else %}
<p>No inter-assay heterodimers flagged.</p>
{% endif %}

<h2>Worst-channel LoD</h2>
<p>The cartridge inherits the worst per-channel LoD:
   {% if panel_lod_worst %}
   channel <code>{{ panel_lod_worst.label }}</code> at
   <code>{{ "%.1f"|format(panel_lod_worst.lod_copies_per_ml) }}</code> copies/mL.
   {% else %}
   (no channels)
   {% endif %}
</p>

<div class="plot"><img src="{{ heatmap }}" alt="LoD bar chart"></div>

<div class="footer">
  Generated by LAMP-Forge v{{ version }} at {{ generated_at }}.
  Manifest hash (SHA-256): <code>{{ manifest_hash }}</code>.
  In silico cartridge build &mdash; confirm cartridge performance empirically
  (per-channel LoD, no-template control, panel co-amplification) before
  releasing the SKU for production.
</div>
</body>
</html>
"""


def render_cartridge_html(result: CartridgeBuildResult, path: Path) -> None:
    """Render the cartridge-level HTML build report (channel table + LoD bar).

    Args:
        result: :class:`CartridgeBuildResult` from :func:`run_cartridge`.
        path: Output HTML path; parent directories are created.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from jinja2 import BaseLoader, Environment

    from lamp_forge import __version__

    verdict_class = {"GO": "go", "WARNING": "warn", "NO_GO": "nogo"}[result.build_status]

    labels = [t["label"] for t in result.per_target]
    lods = [float(t["lod"]["lod_copies_per_ml"]) for t in result.per_target]

    fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(labels) + 2), 3.5))
    ax.bar(labels, lods, color="#0369a1")
    ax.set_ylabel("LoD (copies/mL)")
    ax.set_title("Per-channel LoD (lower is better)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    plt.close(fig)
    heatmap = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")

    env = Environment(loader=BaseLoader(), autoescape=True, trim_blocks=True, lstrip_blocks=True)
    html = env.from_string(_CARTRIDGE_HTML_TEMPLATE).render(
        cartridge_id=result.manifest.cartridge_id,
        chemistry=result.manifest.chemistry,
        n_targets=len(result.manifest.targets),
        max_channels=result.manifest.max_channels,
        device_window=result.manifest.device_window_min,
        build_status=result.build_status,
        verdict_class=verdict_class,
        reasons=result.reasons,
        per_target=result.per_target,
        panel_summary=result.panel_summary,
        panel_lod_worst=result.panel_lod_worst,
        heatmap=heatmap,
        version=__version__,
        generated_at=result.generated_at,
        manifest_hash=result.manifest_hash,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
