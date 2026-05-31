r"""Pre-order readiness check for LAMP primer sets.

Before committing to primer synthesis, a BioVind-style portable assay must
clear three independent hurdles:

1. **LOD** -- will the assay detect the pathogen at the expected sample
   concentration, given the extraction chain?
2. **TTP** -- at the LOD copy count, will the reaction turn positive before
   the device time window closes (30-60 min on BioVind BioID)?
3. **RT-LAMP** -- for RNA targets, are all core primers above the 63 degC
   co-activity floor required for one-step reverse transcription?

Running these checks separately (``lamp-forge lod``, ``lamp-forge ttp``,
``lamp-forge rt-check``) is correct but error-prone: users must remember to
run all three and manually reconcile the outputs.  This module chains them
into a single pass that returns a **GO / WARNING / NO-GO** verdict with
machine-readable and human-readable output.

Typical pre-order workflow::

    lamp-forge preorder \
      --input results/prrsv_orf7/primer_sets.json \
      --na-type rna \
      --sample-volume 1000 \
      --efficiency 0.50 \
      --eluate-volume 100 \
      --reaction-input 5 \
      --preset rt-lamp \
      --device-window 60

References:
    Tanner NA et al. (2012) BioTechniques 53(2):81-89.
    Notomi T et al. (2000) Nucleic Acids Res 28(12):e63.
    NEB WarmStart(R) LAMP Kit manual (E1700).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from lamp_forge.lod import ExtractionParams, estimate_lod
from lamp_forge.rt_lamp import (
    RtLampCheckResult,
    RtLampParams,
    TargetNucleicAcid,
    check_primer_sets_for_rt_lamp,
)
from lamp_forge.ttp import TtpParams, TtpPoint, ttp_at_copies

# ---------------------------------------------------------------------------
# Status enum
# ---------------------------------------------------------------------------


class PreorderStatus(StrEnum):
    """Overall pre-order readiness verdict.

    Attributes:
        GO: All checks pass -- the primer set is ready to order for the
            described device/extraction configuration.
        WARNING: TTP is within the device window, but one or more RT-LAMP
            checks failed (RNA target with sub-optimal Tm primers) or TTP is
            within 5 % of the device window limit.  Review before ordering.
        NO_GO: TTP at LOD_95 exceeds the device window.  The assay will not
            read out in time on the target platform.  Optimise the extraction
            chain or tighten the Tm window before ordering.
    """

    GO = "GO"
    WARNING = "WARNING"
    NO_GO = "NO_GO"


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PreorderResult:
    """Combined pre-order readiness result for a primer_sets.json file.

    Attributes:
        target_na: ``"rna"`` or ``"dna"``; determines whether the RT-LAMP
            check is applicable.
        n_sets_assessed: Number of primer sets evaluated (top-N or all).
        n_sets_rt_ok: Number of sets that pass the RT-LAMP Tm check.  Set to
            -1 for DNA targets, where RT is not required.
        lod_95_copies_per_rxn: Lambda_LOD at P(detect)=0.95 (~2.996 copies).
        lod_95_copies_per_ml: Back-calculated LOD_95 in the original sample
            (copies/mL).
        ttp_at_lod_min: Predicted TTP (minutes) at ``lod_95_copies_per_rxn``
            copies using the selected chemistry preset.
        device_window_min: Device run window (minutes).
        ttp_in_window: True when ``ttp_at_lod_min <= device_window_min``.
        status: GO / WARNING / NO_GO.
        reasons: Human-readable strings explaining the status.
    """

    target_na: str
    n_sets_assessed: int
    n_sets_rt_ok: int
    lod_95_copies_per_rxn: float
    lod_95_copies_per_ml: float
    ttp_at_lod_min: float
    device_window_min: float
    ttp_in_window: bool
    status: PreorderStatus
    reasons: tuple[str, ...]


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

_TTP_NEAR_WINDOW_FRACTION: float = 0.05


def run_preorder(
    sets_data: list[dict[str, Any]],
    extraction: ExtractionParams,
    ttp_params: TtpParams,
    *,
    target_na: TargetNucleicAcid = TargetNucleicAcid.DNA,
    rt_params: RtLampParams | None = None,
    top_n: int | None = None,
) -> PreorderResult:
    """Run the combined pre-order readiness check.

    Chains LOD estimation, TTP prediction, and (for RNA targets) RT-LAMP
    compatibility assessment into a single pass.  The LOD_95 copy count per
    reaction is used as the TTP input point -- this is the worst-case scenario
    the device must handle to call a sample positive at 95% probability.

    Args:
        sets_data: The ``"primer_sets"`` list from a ``primer_sets.json`` file.
        extraction: Sample preparation chain parameters (sample volume,
            extraction efficiency, eluate volume, reaction input volume).
        ttp_params: TTP model parameters (preset or custom) including the
            device window against which TTP is judged.
        target_na: Whether the diagnostic target is DNA or RNA; determines
            whether the RT-LAMP Tm check is performed.
        rt_params: RT-LAMP temperature parameters.  Defaults to standard
            one-step RT-LAMP values (63-65 degC core, 60 degC loop).
        top_n: If given, assess only the first ``top_n`` sets from the list.

    Returns:
        :class:`PreorderResult` summarising the combined readiness verdict.
    """
    # --- LOD_95 ---------------------------------------------------------------
    lod_est = estimate_lod(extraction, detection_probability=0.95)
    lod_rxn = lod_est.lod_copies_per_reaction
    lod_ml = lod_est.lod_copies_per_ml

    # --- TTP at LOD_95 --------------------------------------------------------
    ttp_point: TtpPoint = ttp_at_copies(lod_rxn, ttp_params)
    ttp_in_window = ttp_point.in_device_window

    # --- RT-LAMP check (RNA targets only) -------------------------------------
    rt_results: list[RtLampCheckResult] = []
    n_sets_rt_ok = -1  # sentinel for DNA
    if target_na is TargetNucleicAcid.RNA:
        rt_results = check_primer_sets_for_rt_lamp(
            sets_data,
            target_na=target_na,
            params=rt_params,
            top_n=top_n,
        )
        n_sets_rt_ok = sum(1 for r in rt_results if r.is_rt_optimized)

    n_sets = (
        len(rt_results) if rt_results else (len(sets_data[:top_n]) if top_n else len(sets_data))
    )

    # --- Status logic ---------------------------------------------------------
    reasons: list[str] = []

    if not ttp_in_window:
        status = PreorderStatus.NO_GO
        reasons.append(
            f"TTP at LOD_95 ({ttp_point.ttp_minutes:.1f} min) exceeds device window "
            f"({ttp_params.device_window_min:.0f} min) -- assay will not read out in time."
        )
    else:
        near_limit = ttp_point.ttp_minutes / ttp_params.device_window_min > (
            1.0 - _TTP_NEAR_WINDOW_FRACTION
        )
        rt_not_all_ok = target_na is TargetNucleicAcid.RNA and n_sets_rt_ok < n_sets

        if rt_not_all_ok or near_limit:
            status = PreorderStatus.WARNING
            if near_limit:
                reasons.append(
                    f"TTP at LOD_95 ({ttp_point.ttp_minutes:.1f} min) is within "
                    f"{_TTP_NEAR_WINDOW_FRACTION * 100:.0f}% of device window "
                    f"({ttp_params.device_window_min:.0f} min) -- marginal."
                )
            if rt_not_all_ok:
                n_bad = n_sets - n_sets_rt_ok
                rt_min = (rt_params or RtLampParams()).rt_min_tm
                reasons.append(
                    f"{n_bad}/{n_sets} primer set(s) have core primers below the "
                    f"RT-LAMP Tm floor ({rt_min:.1f} degC) -- re-run with "
                    f"primer.tm_min >= {rt_min:.1f} to fix suboptimal sets."
                )
        else:
            status = PreorderStatus.GO
            if target_na is TargetNucleicAcid.RNA:
                reasons.append(
                    f"All {n_sets} set(s) pass RT-LAMP Tm check and TTP is within window."
                )
            else:
                reasons.append(
                    f"TTP at LOD_95 ({ttp_point.ttp_minutes:.1f} min) is within "
                    f"device window ({ttp_params.device_window_min:.0f} min)."
                )

    return PreorderResult(
        target_na=str(target_na),
        n_sets_assessed=n_sets,
        n_sets_rt_ok=n_sets_rt_ok,
        lod_95_copies_per_rxn=lod_rxn,
        lod_95_copies_per_ml=lod_ml,
        ttp_at_lod_min=ttp_point.ttp_minutes,
        device_window_min=ttp_params.device_window_min,
        ttp_in_window=ttp_in_window,
        status=status,
        reasons=tuple(reasons),
    )


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_preorder_csv(result: PreorderResult, path: Path) -> None:
    """Write a pre-order readiness result to a key-value CSV file.

    Args:
        result: :class:`PreorderResult` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str]] = [
        ("target_na", result.target_na),
        ("n_sets_assessed", str(result.n_sets_assessed)),
        ("n_sets_rt_ok", str(result.n_sets_rt_ok)),
        ("lod_95_copies_per_rxn", f"{result.lod_95_copies_per_rxn:.4f}"),
        ("lod_95_copies_per_ml", f"{result.lod_95_copies_per_ml:.2f}"),
        ("ttp_at_lod_95_min", f"{result.ttp_at_lod_min:.2f}"),
        ("device_window_min", f"{result.device_window_min:.1f}"),
        ("ttp_in_window", "true" if result.ttp_in_window else "false"),
        ("status", str(result.status)),
        ("reasons", "; ".join(result.reasons)),
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["field", "value"])
        writer.writerows(rows)
