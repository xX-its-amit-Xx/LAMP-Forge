"""Reverse-LOD calculator: given a target LOD, what extraction chain achieves it?

Complements :mod:`lamp_forge.lod` (forward direction: extraction params -> LOD)
with the backward direction: target LOD -> required extraction parameters.

Use this before finalising a sample-preparation protocol to answer "which
sample volume and extraction efficiency do I need to hit my required detection
threshold?"

The Poisson model is identical to the forward calculator::

    LOD [copies/mL] = -ln(1 - P) / (V_eff [uL] / 1000)

Rearranging for the required effective volume per reaction::

    V_eff_required [uL] = -ln(1 - P) * 1000 / LOD_target [copies/mL]

where V_eff = V_sample * eta * (V_rxn / V_eluate).

Typical field usage (BioVind oil & gas / farm / POC verticals):

* Oil & gas SRB monitoring: ``--target-lod 500`` (threshold of concern for
  souring is ~10^3 cells/mL; 500 copies/mL gives a comfortable safety margin).
* Farm biosecurity PRRSV oral fluid: ``--target-lod 100`` (clinically
  relevant sensitivity for PRRSV ORF7 RT-LAMP from a rope sample).
* Human POC MTB sputum: ``--target-lod 50`` (approaches the clinical
  sensitivity of Xpert MTB/RIF at ~10 GE/mL in optimised preparations).

References:
    Dingle & Basler (2012) Poisson statistics and end-point digital LAMP.
    Bio-Rad dPCR application notes on Poisson correction.
    LAMP-Forge lod.py for the forward LOD calculator.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from lamp_forge.lod import poisson_lod


@dataclass(frozen=True, slots=True)
class ExtractionSweepPoint:
    """One (sample_volume, efficiency) combination in the reverse-LOD sweep.

    Attributes:
        sample_volume_ul: Sample volume input to extraction (uL).
        extraction_efficiency: Extraction efficiency as a fraction in (0, 1].
        eluate_volume_ul: Extraction eluate volume (uL).
        reaction_input_ul: Volume of eluate added to the LAMP reaction (uL).
        effective_vol_ul: Effective sample per reaction =
            V_sample * eta * (V_rxn / V_eluate), in uL.
        lod_copies_per_ml: Achieved LOD at the requested detection probability,
            in copies per mL of original sample.
        achieves_target: True when lod_copies_per_ml <= target_lod_copies_per_ml.
    """

    sample_volume_ul: float
    extraction_efficiency: float
    eluate_volume_ul: float
    reaction_input_ul: float
    effective_vol_ul: float
    lod_copies_per_ml: float
    achieves_target: bool


def required_effective_vol_ul(
    target_lod_copies_per_ml: float,
    detection_probability: float = 0.95,
) -> float:
    """Minimum effective sample volume per reaction to achieve a target LOD.

    This is the direct analytical inverse of the Poisson LOD formula:
    given a desired LOD in copies/mL, return the minimum effective sample
    volume (V_sample * eta * V_rxn / V_eluate) needed in the reaction.

    Args:
        target_lod_copies_per_ml: Required assay LOD in copies/mL of original
            sample.  Must be strictly positive.
        detection_probability: Detection probability at the LOD, e.g. 0.95 for
            LOD_95.  Must be in (0, 1).

    Returns:
        Required effective sample volume per reaction in uL.

    Raises:
        ValueError: If target_lod_copies_per_ml <= 0 or detection_probability
            is outside (0, 1).
    """
    if target_lod_copies_per_ml <= 0.0:
        raise ValueError(
            f"target_lod_copies_per_ml must be positive, got {target_lod_copies_per_ml!r}"
        )
    lambda_lod = poisson_lod(detection_probability)
    return lambda_lod * 1000.0 / target_lod_copies_per_ml


def extraction_sweep(
    target_lod_copies_per_ml: float,
    sample_volumes_ul: tuple[float, ...] = (200.0, 400.0, 1000.0, 2000.0),
    efficiencies: tuple[float, ...] = (0.30, 0.40, 0.50, 0.60, 0.70),
    eluate_volume_ul: float = 100.0,
    reaction_input_ul: float = 5.0,
    detection_probability: float = 0.95,
) -> list[ExtractionSweepPoint]:
    """Sweep sample volumes and efficiencies to find combinations meeting a target LOD.

    For each (sample_volume, efficiency) pair the achieved LOD is computed and
    compared against the target.  Results are sorted by achieved LOD ascending
    so the rows that meet the target appear first.

    Args:
        target_lod_copies_per_ml: Required assay LOD in copies/mL.
        sample_volumes_ul: Sample volumes to evaluate, in uL.  Defaults cover
            the most common field collection volumes (200 uL blood/swab,
            400 uL VTM swab, 1 mL produced water, 2 mL nasal wash).
        efficiencies: Extraction efficiencies to evaluate as fractions.
            Defaults cover the practical range for rapid-extraction kits
            used in point-of-care and field settings (30-70%).
        eluate_volume_ul: Fixed elution volume for all combinations (uL).
        reaction_input_ul: Fixed volume of eluate added to each LAMP reaction.
        detection_probability: Detection probability, default 0.95 (LOD_95).

    Returns:
        List of :class:`ExtractionSweepPoint`, one per (sample_volume, efficiency)
        pair, sorted by achieved LOD ascending (best / lowest LOD first).

    Raises:
        ValueError: On any invalid parameter.
    """
    if eluate_volume_ul <= 0.0:
        raise ValueError(f"eluate_volume_ul must be positive, got {eluate_volume_ul!r}")
    if reaction_input_ul <= 0.0:
        raise ValueError(f"reaction_input_ul must be positive, got {reaction_input_ul!r}")
    if reaction_input_ul > eluate_volume_ul:
        raise ValueError(
            f"reaction_input_ul ({reaction_input_ul}) must not exceed "
            f"eluate_volume_ul ({eluate_volume_ul})"
        )
    if not sample_volumes_ul:
        raise ValueError("sample_volumes_ul must contain at least one value")
    if not efficiencies:
        raise ValueError("efficiencies must contain at least one value")

    lambda_lod = poisson_lod(detection_probability)
    dilution_factor = reaction_input_ul / eluate_volume_ul

    points: list[ExtractionSweepPoint] = []
    for v_sample in sample_volumes_ul:
        if v_sample <= 0.0:
            raise ValueError(f"Each sample_volume_ul must be positive, got {v_sample!r}")
        for eta in efficiencies:
            if not (0.0 < eta <= 1.0):
                raise ValueError(f"Each efficiency must be in (0, 1], got {eta!r}")
            v_eff = v_sample * eta * dilution_factor
            lod_per_ml = lambda_lod * 1000.0 / v_eff
            points.append(
                ExtractionSweepPoint(
                    sample_volume_ul=v_sample,
                    extraction_efficiency=eta,
                    eluate_volume_ul=eluate_volume_ul,
                    reaction_input_ul=reaction_input_ul,
                    effective_vol_ul=v_eff,
                    lod_copies_per_ml=lod_per_ml,
                    achieves_target=(lod_per_ml <= target_lod_copies_per_ml),
                )
            )

    points.sort(key=lambda p: p.lod_copies_per_ml)
    return points


def write_sweep_csv(
    points: list[ExtractionSweepPoint],
    target_lod_copies_per_ml: float,
    path: Path,
) -> None:
    """Write the extraction sweep table to a CSV file.

    Args:
        points: Sweep results from :func:`extraction_sweep`.
        target_lod_copies_per_ml: The target LOD used for the sweep, written
            into every data row for traceability.
        path: Output CSV path; parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "sample_volume_ul",
                "extraction_efficiency",
                "eluate_volume_ul",
                "reaction_input_ul",
                "effective_vol_ul",
                "lod_copies_per_ml",
                "target_lod_copies_per_ml",
                "achieves_target",
            ]
        )
        for p in points:
            writer.writerow(
                [
                    f"{p.sample_volume_ul:.1f}",
                    f"{p.extraction_efficiency:.2f}",
                    f"{p.eluate_volume_ul:.1f}",
                    f"{p.reaction_input_ul:.1f}",
                    f"{p.effective_vol_ul:.3f}",
                    f"{p.lod_copies_per_ml:.2f}",
                    f"{target_lod_copies_per_ml:.2f}",
                    "yes" if p.achieves_target else "no",
                ]
            )
