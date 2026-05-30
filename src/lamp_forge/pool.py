"""Equimolar multiplex primer-pool calculator.

For portable multiplex platforms (e.g. BioVind BioID, up to ~18 targets),
all primers for every assay are combined in a single pre-pooled working mix.
This module calculates the volume of each primer stock needed to achieve the
correct LAMP stoichiometry for every primer in the combined pool.

**LAMP 10x working-stock concentrations per target** (Li et al. 2018 / NEB E1700)::

    FIP / BIP : 16 uM  ->  1.6 uM final in the 1x reaction
    LF  / LB  :  4 uM  ->  0.4 uM final in the 1x reaction
    F3  / B3  :  2 uM  ->  0.2 uM final in the 1x reaction

For a multiplex panel the same concentrations apply *per target*: when the
pool is used at 1/10 volume in the reaction, every assay's primers contribute
the same copy count as they would in a single-target LAMP.

**Important: stock-concentration feasibility.**
All primer volumes must fit inside the desired total pool volume.  For N
targets with P per-role concentrations the minimum required stock is::

    min_stock_uM = sum(target_conc_i for each primer across all targets)
                 = N x (16+16+4+4+2+2)  [6-primer set]
                 = N x 44 uM

At the IDT/Twist standard of 100 uM this is feasible for up to 2 targets
(88 uM < 100 uM).  For 3+ targets request resuspension at a higher
concentration (200 uM, 500 uM) or contact the vendor for a pre-mixed pool.
:func:`build_pool_plan` raises :class:`ValueError` with the minimum required
concentration so the caller can catch and relay this information.

Usage (library)::

    from lamp_forge.pool import PoolingParams, build_pool_plan, write_pool_csv
    import json
    from pathlib import Path

    panel = json.loads(Path("results/souring_panel/panel.json").read_text())
    params = PoolingParams(stock_conc_um=100.0, total_pool_volume_ul=500.0)
    plan = build_pool_plan(panel["selection"], params)
    write_pool_csv(plan, Path("results/souring_panel/pool_sheet.csv"))

References:
    Li J et al. (2018) J Vis Exp 139:e58547. LAMP primer design and pooling.
    NEB WarmStart(R) LAMP Kit protocol (E1700). 63-65 degC one-step protocol.
    Notomi T et al. (2000) Nucleic Acids Res 28(12):e63.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Role-specific concentrations in the 10x working-stock pool (uM).
# FIP/BIP are supplied at 8x molar excess relative to outer primers because
# they are consumed faster (the inner primers initiate most amplification
# cycles once the dumbbell intermediates form).
# ---------------------------------------------------------------------------

_ROLE_CONC_UM: dict[str, float] = {
    "FIP": 16.0,
    "BIP": 16.0,
    "LF": 4.0,
    "LB": 4.0,
    "F3": 2.0,
    "B3": 2.0,
}

_ROLE_ORDER: tuple[str, ...] = ("F3", "B3", "FIP", "BIP", "LF", "LB")

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PoolingParams:
    """Parameters for building a multiplex primer working stock.

    Attributes:
        stock_conc_um: Concentration of each synthesised primer stock tube (uM).
            IDT and Twist ship standard oligos at 100 uM; higher-concentration
            resuspension (200 uM, 500 uM) is available on request and required
            for panels with three or more targets.
        total_pool_volume_ul: Desired total volume of the pooled working mix (uL).
            500 uL is typical for a small lab batch.  The pool is a 10x working
            stock; add 2.5 uL per 25 uL LAMP reaction.
    """

    stock_conc_um: float = 100.0
    total_pool_volume_ul: float = 500.0

    def __post_init__(self) -> None:
        """Validate parameter ranges."""
        if self.stock_conc_um <= 0.0:
            raise ValueError(f"stock_conc_um must be positive, got {self.stock_conc_um!r}")
        if self.total_pool_volume_ul <= 0.0:
            raise ValueError(
                f"total_pool_volume_ul must be positive, got {self.total_pool_volume_ul!r}"
            )


@dataclass(frozen=True, slots=True)
class PoolEntry:
    """One primer's contribution to the pooled working stock.

    Each row corresponds to one pipetting step: transfer ``vol_stock_ul`` uL
    from the labelled synthesis stock tube into the pool vessel.

    Attributes:
        target_label: Short name for the detection target (e.g. ``'SRB_dsrB'``).
        primer_role: LAMP primer role: ``'F3'``, ``'B3'``, ``'FIP'``, ``'BIP'``,
            ``'LF'``, or ``'LB'``.
        primer_name: Unique name for this oligo in the pool.
        sequence: 5'-to-3' primer sequence (DNA, uppercase).
        stock_conc_um: Concentration of the synthesis stock tube (uM).
        target_conc_um: Required concentration of this primer in the pool (uM).
        vol_stock_ul: Microlitres to pipette from the stock tube into the pool.
    """

    target_label: str
    primer_role: str
    primer_name: str
    sequence: str
    stock_conc_um: float
    target_conc_um: float
    vol_stock_ul: float


@dataclass(frozen=True, slots=True)
class PoolPlan:
    """Complete pipetting recipe for a multiplex LAMP primer working mix.

    Attributes:
        entries: One entry per primer, sorted by target then by canonical LAMP
            role order (F3, B3, FIP, BIP, LF, LB).
        total_pool_volume_ul: Total volume of the working pool (uL).
        water_volume_ul: Volume of nuclease-free water to add after all primer
            stocks have been pipetted.  Always >= 0 for a valid plan.
        n_targets: Number of detection targets in the panel.
        n_primers: Total number of oligos in the pool.
    """

    entries: tuple[PoolEntry, ...]
    total_pool_volume_ul: float
    water_volume_ul: float
    n_targets: int
    n_primers: int


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def build_pool_plan(
    selection: list[dict[str, Any]],
    params: PoolingParams | None = None,
) -> PoolPlan:
    """Build a pipetting plan for a multiplex primer working stock.

    Reads the ``selection`` list from a ``panel.json`` produced by
    ``lamp-forge panel`` and returns a :class:`PoolPlan` describing how much
    of each individual primer stock to pipette into a combined working mix.

    Args:
        selection: List of per-target selection dicts, each with fields
            ``target`` (str), ``set_id`` (str), and ``oligos``
            (list of ``{"role": ..., "sequence": ...}``).  This is the
            ``"selection"`` array in ``panel.json``.
        params: Pooling parameters (stock concentration and total pool volume).
            Defaults to :class:`PoolingParams` with 100 uM stock and 500 uL pool.

    Returns:
        :class:`PoolPlan` with per-primer pipetting volumes and the water top-up.

    Raises:
        ValueError: If the stock concentration is insufficient for the panel size,
            i.e. the sum of all primer working-stock concentrations exceeds
            ``params.stock_conc_um``.  The error message states the minimum
            required stock concentration.
        ValueError: If ``selection`` is empty.
    """
    if params is None:
        params = PoolingParams()

    if not selection:
        raise ValueError("selection must not be empty.")

    entries: list[PoolEntry] = []

    for item in selection:
        target = str(item["target"])
        oligos_raw = item.get("oligos", [])
        if not isinstance(oligos_raw, list):
            raise ValueError(
                f"Expected 'oligos' list for target {target!r}, got {type(oligos_raw).__name__}"
            )

        oligos_by_role: dict[str, str] = {}
        for o in oligos_raw:
            if not isinstance(o, dict):
                continue
            role = str(o.get("role", ""))
            seq = str(o.get("sequence", ""))
            if role and seq:
                oligos_by_role[role] = seq.upper()

        for role in _ROLE_ORDER:
            role_seq: str | None = oligos_by_role.get(role)
            if role_seq is None:
                continue
            seq = role_seq
            target_conc = _ROLE_CONC_UM[role]
            vol_stock = target_conc * params.total_pool_volume_ul / params.stock_conc_um
            entries.append(
                PoolEntry(
                    target_label=target,
                    primer_role=role,
                    primer_name=f"{target}_{role}",
                    sequence=seq,
                    stock_conc_um=params.stock_conc_um,
                    target_conc_um=target_conc,
                    vol_stock_ul=vol_stock,
                )
            )

    total_stock_vol = sum(e.vol_stock_ul for e in entries)
    water_volume = params.total_pool_volume_ul - total_stock_vol

    if water_volume < 0.0:
        sum_conc = sum(e.target_conc_um for e in entries)
        raise ValueError(
            f"Stock concentration ({params.stock_conc_um:.0f} uM) is too low for "
            f"this {len(selection)}-target panel. "
            f"The sum of all primer working-stock concentrations is {sum_conc:.0f} uM "
            f"(> {params.stock_conc_um:.0f} uM). "
            f"Request resuspension at >= {int(sum_conc) + 1} uM "
            "(IDT/Twist offer 200 uM and 500 uM options)."
        )

    n_targets = len({e.target_label for e in entries})

    return PoolPlan(
        entries=tuple(entries),
        total_pool_volume_ul=params.total_pool_volume_ul,
        water_volume_ul=water_volume,
        n_targets=n_targets,
        n_primers=len(entries),
    )


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_pool_csv(plan: PoolPlan, path: Path) -> None:
    """Write the pooling plan to a CSV pipetting sheet.

    Each row is one pipetting step.  A final WATER row records the volume of
    nuclease-free water to top up the pool to the desired total volume.

    Args:
        plan: The pooling plan returned by :func:`build_pool_plan`.
        path: Output file path.  Parent directories are created as needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "target",
                "role",
                "primer_name",
                "sequence",
                "stock_conc_um",
                "target_conc_um",
                "vol_stock_ul",
            ]
        )
        for e in plan.entries:
            writer.writerow(
                [
                    e.target_label,
                    e.primer_role,
                    e.primer_name,
                    e.sequence,
                    f"{e.stock_conc_um:.1f}",
                    f"{e.target_conc_um:.1f}",
                    f"{e.vol_stock_ul:.3f}",
                ]
            )
        writer.writerow(
            [
                "WATER",
                "",
                "Nuclease-free water",
                "",
                "",
                "",
                f"{plan.water_volume_ul:.3f}",
            ]
        )
