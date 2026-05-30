"""YAML config loading and validation.

We deliberately keep validation simple: raise ``ValueError`` with a precise
message at the first bad field. The CLI catches and prints these so a
non-bioinformatician gets a useful error instead of a stack trace.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml

from lamp_forge.types import PipelineConfig


class ConfigError(ValueError):
    """Raised when the YAML config is malformed or contains impossible values."""


def _require(d: dict[str, Any], key: str, ctx: str) -> Any:
    """Pull ``key`` from ``d`` or raise a ``ConfigError`` naming the section."""
    if key not in d:
        raise ConfigError(f"Missing required field '{key}' in [{ctx}] section of config")
    return d[key]


def _check_range(name: str, value: float, low: float, high: float) -> None:
    """Raise if ``value`` falls outside [low, high]."""
    if not (low <= value <= high):
        raise ConfigError(f"{name}={value} out of allowed range [{low}, {high}]")


def load_config(path: str | Path) -> PipelineConfig:
    """Load and validate a LAMP-Forge YAML config.

    Args:
        path: Path to the YAML file.

    Returns:
        A fully populated :class:`PipelineConfig`.

    Raises:
        ConfigError: If the file is missing, malformed, or contains values
            that would produce an impossible design (e.g. ``tm_min > tm_max``).
    """
    path = Path(path)
    if not path.exists():
        raise ConfigError(f"Config file not found: {path}")

    with path.open() as fh:
        raw = yaml.safe_load(fh)

    if not isinstance(raw, dict):
        raise ConfigError(f"Config root must be a mapping, got {type(raw).__name__}")

    target = _require(raw, "target", "root")
    off_targets = _require(raw, "off_targets", "root")
    conservation = _require(raw, "conservation", "root")
    primer = _require(raw, "primer", "root")
    output = _require(raw, "output", "root")

    taxon_id = target.get("taxon_id")
    accessions = target.get("accessions", []) or []
    if not taxon_id and not accessions:
        raise ConfigError(
            "[target] must specify either 'taxon_id' (preferred) or 'accessions' list"
        )

    email = target.get("email") or os.environ.get("NCBI_EMAIL", "")
    if not email:
        raise ConfigError(
            "[target.email] is required (or set NCBI_EMAIL env var). "
            "NCBI requires an email for Entrez requests."
        )

    tm_min = float(_require(primer, "tm_min", "primer"))
    tm_max = float(_require(primer, "tm_max", "primer"))
    if tm_min >= tm_max:
        raise ConfigError(f"primer.tm_min ({tm_min}) must be < primer.tm_max ({tm_max})")
    _check_range("primer.tm_min", tm_min, 40.0, 80.0)
    _check_range("primer.tm_max", tm_max, 40.0, 80.0)

    gc_min = float(_require(primer, "gc_min", "primer"))
    gc_max = float(_require(primer, "gc_max", "primer"))
    if gc_min >= gc_max:
        raise ConfigError(f"primer.gc_min ({gc_min}) must be < primer.gc_max ({gc_max})")

    f2b2 = primer.get("amplicon_size", {})
    f2_b2_min = int(f2b2.get("f2_b2_min", 120))
    f2_b2_max = int(f2b2.get("f2_b2_max", 160))
    if f2_b2_min >= f2_b2_max:
        raise ConfigError("primer.amplicon_size.f2_b2_min must be < f2_b2_max")
    if f2_b2_min < 80:
        raise ConfigError(
            f"f2_b2_min={f2_b2_min} is below the LAMP geometric floor (~80bp); "
            "primers will overlap"
        )

    min_region_length = int(_require(conservation, "min_region_length", "conservation"))
    if min_region_length < f2_b2_max + 40:
        raise ConfigError(
            f"conservation.min_region_length={min_region_length} is too short for "
            f"f2_b2_max={f2_b2_max}; need at least f2_b2_max + 40bp to flank F3/B3"
        )

    return PipelineConfig(
        target_name=str(_require(target, "name", "target")),
        taxon_id=int(taxon_id) if taxon_id else None,
        accessions=list(accessions),
        gene=target.get("gene"),
        max_sequences=int(target.get("max_sequences", 20)),
        email=email,
        off_target_dir=Path(_require(off_targets, "fasta_dir", "off_targets")),
        min_identity_threshold=float(off_targets.get("min_identity_threshold", 0.80)),
        min_coverage_threshold=float(off_targets.get("min_coverage_threshold", 0.80)),
        window_size=int(conservation.get("window_size", 30)),
        entropy_threshold=float(conservation.get("entropy_threshold", 0.20)),
        min_region_length=min_region_length,
        tm_min=tm_min,
        tm_max=tm_max,
        tm_match_tolerance=float(primer.get("tm_match_tolerance", 2.0)),
        gc_min=gc_min,
        gc_max=gc_max,
        hairpin_dg_threshold=float(primer.get("hairpin_dg_threshold", -2.0)),
        dimer_dg_threshold=float(primer.get("dimer_dg_threshold", -5.0)),
        f2_b2_min=f2_b2_min,
        f2_b2_max=f2_b2_max,
        output_dir=Path(_require(output, "dir", "output")),
        top_n=int(output.get("top_n", 10)),
        generate_html=bool(output.get("generate_html", True)),
        generate_csv=bool(output.get("generate_csv", True)),
        cache_dir=Path(os.environ.get("LAMP_FORGE_CACHE", "cache")),
    )


def check_config_warnings(
    config: PipelineConfig,
    *,
    check_dirs: bool = True,
) -> list[str]:
    """Return a list of soft-warning strings for a validated config.

    These are not hard errors (``load_config`` already enforces those) but
    common parameter choices that often indicate a misconfiguration or that
    will produce suboptimal designs.  Returns an empty list when the config
    looks reasonable.

    Args:
        config: A :class:`~lamp_forge.types.PipelineConfig` already returned
            by :func:`load_config`.
        check_dirs: When *True* (default), warn if the off-target FASTA
            directory does not exist on the filesystem.

    Returns:
        A (possibly empty) list of human-readable warning strings, each
        beginning with a short tag such as ``[off_targets]`` or ``[primer]``.
    """
    warnings: list[str] = []

    if check_dirs and not config.off_target_dir.exists():
        warnings.append(
            f"[off_targets] fasta_dir '{config.off_target_dir}' does not exist "
            f"-- the pipeline will fail when it tries to build the BLAST database"
        )

    if config.entropy_threshold < 0.08:
        warnings.append(
            f"[conservation] entropy_threshold {config.entropy_threshold:.3f} bits is very "
            f"tight; most targets produce no conserved regions below 0.08 bits"
        )
    elif config.entropy_threshold > 0.50:
        warnings.append(
            f"[conservation] entropy_threshold {config.entropy_threshold:.3f} bits is very "
            f"loose; may include poorly-conserved positions -- consider <= 0.40 bits for "
            f"broad assay inclusivity"
        )

    if config.max_sequences < 5:
        warnings.append(
            f"[target] max_sequences={config.max_sequences} is very low; conservation "
            f"analysis needs >= 5 sequences to be meaningful -- increase to at least 10"
        )

    if config.tm_min < 58.0:
        warnings.append(
            f"[primer] tm_min={config.tm_min:.1f} degC is below the practical LAMP floor "
            f"(~58 degC); Bst polymerase strand-displacement activity drops sharply below "
            f"this temperature"
        )
    if config.tm_max > 72.0:
        warnings.append(
            f"[primer] tm_max={config.tm_max:.1f} degC is above 72 degC; primers this hot "
            f"are difficult to design in the GC-range specified and may mis-prime on "
            f"secondary structures"
        )

    if config.top_n > 50:
        warnings.append(
            f"[output] top_n={config.top_n} is very large; consider <= 20 to keep the "
            f"HTML report navigable"
        )

    if config.min_identity_threshold < 0.70:
        warnings.append(
            f"[off_targets] min_identity_threshold={config.min_identity_threshold:.2f} is "
            f"below 70%; at this level, near-random nucleotide matches will be flagged as "
            f"cross-reactive -- consider >= 0.75"
        )

    return warnings
