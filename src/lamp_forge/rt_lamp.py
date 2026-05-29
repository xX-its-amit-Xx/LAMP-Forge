"""RT-LAMP compatibility checker for RNA-target LAMP assays.

Standard LAMP runs at 60-65 degC on DNA templates.  RT-LAMP adds a
reverse-transcription step so the same isothermal reaction can detect RNA
targets — RNA viruses, transcript-level markers, etc.  One-step RT-LAMP
combines a heat-stable reverse transcriptase with a strand-displacing
polymerase (e.g. Bst 2.0 WarmStart) in a single 63-65 degC reaction.

Primers designed for the broad 60-65 degC LAMP window may have Tm values
that are acceptable for DNA-LAMP but suboptimal for one-step RT-LAMP:
reverse transcriptases used in one-step protocols (NEB RTx, AMV-RT) have
reduced processivity below ~63 degC, which can slow cDNA synthesis and
reduce sensitivity.

This module provides a lightweight, primer3-free assessment of whether a
primer set produced by ``lamp-forge run`` is ready for one-step RT-LAMP, or
whether the Tm window should be tightened (``primer.tm_min >= 63.0``) before
ordering.

Key farm-animal RNA virus targets this module is designed to support:

* PRRSV — porcine reproductive and respiratory syndrome virus (ORF7 / ORF6)
* FMDV — foot-and-mouth disease virus (1D/VP1 or 3D RNA polymerase)
* Avian influenza A H5N1 / H7N9 — HA segment or M-segment universal marker
* Classical swine fever virus — RNA pestivirus, differential to ASFV

References:
    Tanner NA et al. (2015) BioTechniques 58(2):59-68.
    Notomi T et al. (2000) Nucleic Acids Res 28(12):e63.
    NEB WarmStart(R) RTAMP Kit manual (E1700): 63-65 degC one-step protocol.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any


class TargetNucleicAcid(StrEnum):
    """Whether the diagnostic target molecule is DNA or RNA.

    LAMP primer design operates on double-stranded DNA.  RNA targets
    (viruses, mRNA markers) must be reverse-transcribed first; ``lamp-forge
    run`` produces primers from the cDNA representation, and this flag tells
    the RT-check step whether a one-step reverse transcription is required.
    """

    DNA = "dna"
    RNA = "rna"


_DEFAULT_RT_MIN: float = 63.0
_DEFAULT_RT_MAX: float = 65.0
_DEFAULT_LOOP_MIN: float = 60.0

# LF/LB bind single-stranded loop structures and tolerate a wider Tm range
# than the core strand-displacing primers.
_LOOP_ROLES: frozenset[str] = frozenset({"LF", "LB"})
_PRIMER_FIELDS: tuple[str, ...] = ("f3", "b3", "fip", "bip", "lf", "lb")


@dataclass(frozen=True, slots=True)
class RtLampParams:
    """Temperature parameters for one-step RT-LAMP optimisation.

    Attributes:
        rt_min_tm: Minimum Tm (degC) for core primers (F3, B3, FIP, BIP).
            Below this value the reverse transcriptase may have reduced
            co-activity with Bst polymerase.  Default: 63.0 degC.
        rt_max_tm: Maximum primer Tm (degC).  Default: 65.0 degC.
        loop_min_tm: Minimum Tm accepted for loop primers (LF, LB).
            Default: 60.0 degC.
    """

    rt_min_tm: float = _DEFAULT_RT_MIN
    rt_max_tm: float = _DEFAULT_RT_MAX
    loop_min_tm: float = _DEFAULT_LOOP_MIN

    def __post_init__(self) -> None:  # noqa: D105
        if self.rt_min_tm >= self.rt_max_tm:
            raise ValueError(
                f"rt_min_tm ({self.rt_min_tm}) must be < rt_max_tm ({self.rt_max_tm})"
            )
        if self.loop_min_tm > self.rt_min_tm:
            raise ValueError(
                f"loop_min_tm ({self.loop_min_tm}) must be <= rt_min_tm ({self.rt_min_tm})"
            )


@dataclass(frozen=True, slots=True)
class RtLampCheckResult:
    """RT-LAMP compatibility assessment for one primer set.

    Attributes:
        set_id: Identifier of the assessed primer set.
        target_na: Whether the target is DNA or RNA.
        n_primers: Total primers in the set (4 without loops, 6 with).
        n_in_rt_range: Primers with Tm in the acceptable range for their role.
        n_core_below_rt_min: Core primers (F3/B3/FIP/BIP) with Tm below
            ``rt_min_tm``.  Zero means the set passes the RT-LAMP Tm check.
        is_rt_optimized: True when all core primers meet the RT Tm floor (RNA
            targets), or when the target is DNA (RT not required).
        warnings: Per-primer diagnostic messages for out-of-range Tms.
        recommended_action: Short guidance string; empty when the set is clean
            or the target is DNA.
    """

    set_id: str
    target_na: TargetNucleicAcid
    n_primers: int
    n_in_rt_range: int
    n_core_below_rt_min: int
    is_rt_optimized: bool
    warnings: tuple[str, ...]
    recommended_action: str


def check_primer_set_for_rt_lamp(
    set_dict: dict[str, Any],
    *,
    target_na: TargetNucleicAcid = TargetNucleicAcid.RNA,
    params: RtLampParams | None = None,
) -> RtLampCheckResult:
    """Assess one primer set (from primer_sets.json) for RT-LAMP readiness.

    Core primers (F3/B3/FIP/BIP) must meet ``params.rt_min_tm``; loop
    primers (LF/LB) only need to reach ``params.loop_min_tm``.  DNA targets
    skip the Tm check entirely — RT is not required.

    Args:
        set_dict: Primer set dict as decoded from a ``primer_sets.json`` file.
        target_na: DNA or RNA target type; see :class:`TargetNucleicAcid`.
        params: RT-LAMP temperature parameters.  Defaults to standard
            one-step RT-LAMP values (63-65 degC core, 60 degC loop).

    Returns:
        :class:`RtLampCheckResult` summarising the set's RT-LAMP suitability.
    """
    if params is None:
        params = RtLampParams()

    set_id = str(set_dict.get("set_id", "unknown"))

    if target_na is TargetNucleicAcid.DNA:
        n_primers = sum(1 for k in _PRIMER_FIELDS if set_dict.get(k))
        return RtLampCheckResult(
            set_id=set_id,
            target_na=target_na,
            n_primers=n_primers,
            n_in_rt_range=n_primers,
            n_core_below_rt_min=0,
            is_rt_optimized=True,
            warnings=(),
            recommended_action="DNA target: RT step not required; standard LAMP applies.",
        )

    warnings_out: list[str] = []
    n_primers = 0
    n_in_rt_range = 0
    n_core_below_rt_min = 0

    for field_key in _PRIMER_FIELDS:
        primer = set_dict.get(field_key)
        if not primer:
            continue
        n_primers += 1
        role = str(primer.get("role", field_key.upper())).upper()
        tm = float(primer.get("tm", 0.0))
        is_loop = role in _LOOP_ROLES
        min_tm = params.loop_min_tm if is_loop else params.rt_min_tm

        if min_tm <= tm <= params.rt_max_tm:
            n_in_rt_range += 1
        elif tm < min_tm:
            if not is_loop:
                n_core_below_rt_min += 1
            warnings_out.append(
                f"{role} Tm {tm:.1f} degC < {min_tm:.1f} degC RT-LAMP minimum"
            )
        else:
            warnings_out.append(
                f"{role} Tm {tm:.1f} degC > {params.rt_max_tm:.1f} degC RT-LAMP maximum"
            )

    is_rt_optimized = n_core_below_rt_min == 0 and n_primers > 0

    if is_rt_optimized:
        recommended_action = ""
    elif n_core_below_rt_min == 1:
        recommended_action = (
            f"Tighten primer.tm_min to >= {params.rt_min_tm:.1f} and re-run; "
            "1 core primer is below the RT-LAMP Tm floor."
        )
    else:
        recommended_action = (
            f"Tighten primer.tm_min to >= {params.rt_min_tm:.1f} and re-run; "
            f"{n_core_below_rt_min} core primers are below the RT-LAMP Tm floor."
        )

    return RtLampCheckResult(
        set_id=set_id,
        target_na=target_na,
        n_primers=n_primers,
        n_in_rt_range=n_in_rt_range,
        n_core_below_rt_min=n_core_below_rt_min,
        is_rt_optimized=is_rt_optimized,
        warnings=tuple(warnings_out),
        recommended_action=recommended_action,
    )


def check_primer_sets_for_rt_lamp(
    sets_data: list[dict[str, Any]],
    *,
    target_na: TargetNucleicAcid = TargetNucleicAcid.RNA,
    params: RtLampParams | None = None,
    top_n: int | None = None,
) -> list[RtLampCheckResult]:
    """Assess multiple primer sets from a primer_sets.json payload.

    Args:
        sets_data: The ``"primer_sets"`` list from a ``primer_sets.json`` file.
        target_na: RNA or DNA target type; applied uniformly to all sets.
        params: RT-LAMP temperature parameters.  Defaults to standard
            one-step RT-LAMP values.
        top_n: If given, assess only the first ``top_n`` sets (front of list
            = highest composite score).

    Returns:
        :class:`RtLampCheckResult` for each set, in input order.
    """
    subset = sets_data if top_n is None else sets_data[:top_n]
    return [
        check_primer_set_for_rt_lamp(s, target_na=target_na, params=params)
        for s in subset
    ]


def write_rt_check_csv(results: list[RtLampCheckResult], path: Path) -> None:
    """Write RT-LAMP check results to a CSV file.

    Args:
        results: RT-LAMP check results to serialise.
        path: Output file path.  Parent directories are created as needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "set_id",
            "target_na",
            "n_primers",
            "n_in_rt_range",
            "n_core_below_rt_min",
            "is_rt_optimized",
            "warnings",
            "recommended_action",
        ])
        for r in results:
            writer.writerow([
                r.set_id,
                str(r.target_na),
                r.n_primers,
                r.n_in_rt_range,
                r.n_core_below_rt_min,
                r.is_rt_optimized,
                "; ".join(r.warnings),
                r.recommended_action,
            ])
