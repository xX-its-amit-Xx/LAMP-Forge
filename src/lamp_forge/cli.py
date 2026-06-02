"""Click-based CLI: ``lamp-forge run --config config.yaml`` and helpers.

The CLI is deliberately thin — every subcommand maps to one module function so
that wet-lab users see the same code path whether they invoke from the shell,
from Snakemake, or from the notebook walkthrough.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click
import yaml

from lamp_forge import __version__
from lamp_forge.config import ConfigError, load_config


def _configure_logging(verbose: bool) -> None:
    """Set up root logger with INFO/DEBUG based on ``--verbose``."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="lamp-forge")
@click.option("-v", "--verbose", is_flag=True, help="Enable DEBUG logging.")
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """LAMP-Forge: reproducible LAMP primer design for microbial diagnostics."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    _configure_logging(verbose)


@cli.command()
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to the YAML pipeline config.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Override config.output.dir.",
)
def run(config_path: Path, output_dir: Path | None) -> None:
    """Run the full pipeline end-to-end (fetch → align → design → screen → report)."""
    from lamp_forge.pipeline import run_pipeline

    try:
        config = load_config(config_path)
    except ConfigError as e:
        click.secho(f"Config error: {e}", fg="red", err=True)
        sys.exit(2)

    if output_dir is not None:
        config.output_dir = output_dir

    try:
        sets = run_pipeline(config)
    except Exception as e:
        click.secho(f"Pipeline failed: {e}", fg="red", err=True)
        if click.get_current_context().obj.get("verbose"):
            raise
        sys.exit(1)

    click.secho(
        f"Done. {len(sets)} primer sets written to {config.output_dir}",
        fg="green",
    )
    html_path = config.output_dir / "lamp_forge_report.html"
    if html_path.exists():
        click.echo(f"Open the HTML report: {html_path}")


@cli.command()
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
)
@click.option(
    "--out",
    "out_path",
    type=click.Path(path_type=Path),
    default=Path("results/sequences.fasta"),
    help="Where to write the fetched FASTA.",
)
def fetch(config_path: Path, out_path: Path) -> None:
    """Fetch sequences from NCBI according to the config (no downstream steps)."""
    from lamp_forge import fetch as fetch_mod

    config = load_config(config_path)
    records = fetch_mod.fetch_for_config(config)
    fetch_mod.write_fasta(records, out_path)
    click.echo(f"Fetched {len(records)} sequence(s) → {out_path}")


@cli.command()
@click.option(
    "--input",
    "input_fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Unaligned FASTA file.",
)
@click.option(
    "--output",
    "output_fasta",
    type=click.Path(path_type=Path),
    required=True,
    help="Aligned FASTA output path.",
)
@click.option("--strategy", default="--auto", show_default=True, help="MAFFT strategy flag.")
@click.option("--threads", default=1, type=int, show_default=True)
def align(input_fasta: Path, output_fasta: Path, strategy: str, threads: int) -> None:
    """Run MAFFT alignment on a FASTA file."""
    from lamp_forge import align as align_mod

    align_mod.run_mafft(input_fasta, output_fasta, strategy=strategy, threads=threads)
    click.echo(f"Alignment written to {output_fasta}")


@cli.command()
@click.option(
    "--alignment",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Aligned FASTA file.",
)
@click.option("--window-size", default=30, type=int, show_default=True)
@click.option("--entropy-threshold", default=0.20, type=float, show_default=True)
@click.option("--min-region-length", default=200, type=int, show_default=True)
@click.option(
    "--out-tsv",
    type=click.Path(path_type=Path),
    default=Path("conservation.tsv"),
    show_default=True,
)
def conserve(
    alignment: Path,
    window_size: int,
    entropy_threshold: float,
    min_region_length: int,
    out_tsv: Path,
) -> None:
    """Compute per-position conservation entropy and detect conserved regions."""
    from lamp_forge import align as align_mod
    from lamp_forge import conserve as conserve_mod
    from lamp_forge.report import write_conservation_tsv

    msa = align_mod.load_alignment(alignment)
    track = conserve_mod.compute_track(msa, window_size=window_size)
    regions = conserve_mod.find_conserved_regions(
        track,
        entropy_threshold=entropy_threshold,
        min_region_length=min_region_length,
    )
    write_conservation_tsv(track, out_tsv)
    click.echo(f"Wrote {len(track.raw_entropy)} positions → {out_tsv}")
    click.echo(f"Detected {len(regions)} conserved region(s) ≥{min_region_length}bp")
    for r in regions:
        click.echo(
            f"  {r.region_id}: {r.start}-{r.end} ({r.length}bp, "
            f"mean entropy {r.mean_entropy:.3f} bits)"
        )


@cli.command()
@click.option(
    "--set",
    "set_specs",
    multiple=True,
    required=True,
    metavar="LABEL=PATH",
    help=(
        "Target label and its primer_sets.json, e.g. "
        "--set ASFV=results/asfv/primer_sets.json. Repeat once per target."
    ),
)
@click.option(
    "--out",
    "out_dir",
    type=click.Path(path_type=Path),
    default=Path("results/panel"),
    show_default=True,
    help="Output directory for panel.json / panel_primers.csv / panel_report.html.",
)
@click.option(
    "--top-per-target",
    default=5,
    type=int,
    show_default=True,
    help="Candidate sets considered per target (caps combinatorial search).",
)
@click.option(
    "--dimer-dg-threshold",
    default=-5.0,
    type=float,
    show_default=True,
    help="Inter-assay heterodimer ΔG (kcal/mol) below which a pair is flagged.",
)
@click.option("--no-html", is_flag=True, default=False, help="Skip the HTML report.")
def panel(
    set_specs: tuple[str, ...],
    out_dir: Path,
    top_per_target: int,
    dimer_dg_threshold: float,
    no_html: bool,
) -> None:
    """Check multiplex compatibility across independently-designed targets.

    Each ``--set LABEL=PATH`` points at a primer_sets.json produced by a prior
    ``run``. Selects one set per target that minimises inter-assay
    cross-dimerisation and writes a panel report.
    """
    from lamp_forge.panel import run_panel

    specs: dict[str, Path] = {}
    for raw in set_specs:
        label, sep, raw_path = raw.partition("=")
        label = label.strip()
        if not sep or not label or not raw_path.strip():
            click.secho(f"--set must be LABEL=PATH, got: {raw!r}", fg="red", err=True)
            sys.exit(2)
        if label in specs:
            click.secho(f"Duplicate target label: {label!r}", fg="red", err=True)
            sys.exit(2)
        path = Path(raw_path.strip())
        if not path.exists():
            click.secho(f"primer_sets.json not found for {label!r}: {path}", fg="red", err=True)
            sys.exit(2)
        specs[label] = path

    if len(specs) < 2:
        click.secho("A multiplex panel needs at least 2 --set targets.", fg="red", err=True)
        sys.exit(2)

    try:
        result = run_panel(
            specs,
            out_dir,
            top_per_target=top_per_target,
            dimer_dg_threshold=dimer_dg_threshold,
            generate_html=not no_html,
        )
    except Exception as e:
        click.secho(f"Panel analysis failed: {e}", fg="red", err=True)
        if click.get_current_context().obj.get("verbose"):
            raise
        sys.exit(1)

    if result.is_clean:
        click.secho(
            f"Compatible panel: {len(result.selection)} targets, "
            f"worst inter-assay dG {result.worst_dg:.2f} kcal/mol ({result.search_mode} search).",
            fg="green",
        )
    else:
        click.secho(
            f"Panel flags {len(result.flagged)} inter-assay heterodimer(s); "
            f"worst dG {result.worst_dg:.2f} kcal/mol. See report before ordering.",
            fg="yellow",
        )
    click.echo(f"Outputs written to {out_dir}")
    html_path = out_dir / "panel_report.html"
    if html_path.exists():
        click.echo(f"Open the panel report: {html_path}")


@cli.command()
@click.option(
    "--sample-volume",
    "sample_volume_ul",
    type=float,
    required=True,
    help="Sample volume input to extraction (uL).",
)
@click.option(
    "--efficiency",
    "extraction_efficiency",
    type=float,
    default=0.50,
    show_default=True,
    help="Extraction efficiency as a fraction 0-1 (e.g. 0.50 for 50%).",
)
@click.option(
    "--eluate-volume",
    "eluate_volume_ul",
    type=float,
    required=True,
    help="Volume of extraction eluate (uL).",
)
@click.option(
    "--reaction-input",
    "reaction_input_ul",
    type=float,
    required=True,
    help="Eluate added to each LAMP reaction (uL).",
)
@click.option(
    "--probability",
    "probabilities",
    type=float,
    multiple=True,
    default=(0.90, 0.95, 0.99, 0.999),
    show_default=True,
    help="Detection probability threshold (repeat for multiple, e.g. --probability 0.95 --probability 0.99).",
)
@click.option(
    "--copies-per-cell",
    "copies_per_cell",
    type=int,
    default=1,
    show_default=True,
    help=(
        "Gene copy number per organism cell (default 1 for single-copy genes). "
        "Set to 4-10 for the 16S rRNA sample-adequacy control channel, "
        "or to the plasmid copy number for resistance-gene targets "
        "(e.g. blaKPC on a high-copy plasmid). "
        "When > 1, an additional LOD (cells/mL) column is shown."
    ),
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the LOD table to a CSV file.",
)
def lod(
    sample_volume_ul: float,
    extraction_efficiency: float,
    eluate_volume_ul: float,
    reaction_input_ul: float,
    probabilities: tuple[float, ...],
    copies_per_cell: int,
    out_csv: Path | None,
) -> None:
    r"""Estimate LAMP assay limit of detection (LOD) across the extraction chain.

    Computes the LOD in copies/reaction and back-calculates to copies/mL in
    the original sample using Poisson single-molecule statistics.

    Use --copies-per-cell > 1 for multi-copy gene targets to also show the
    LOD expressed in cells/mL.  The 16S rRNA control channel typically has
    4-10 copies per cell; --copies-per-cell 7 is a conservative midpoint.

    \b
    Example (200 uL blood, 50% extraction into 50 uL eluate, 5 uL to reaction):
        lamp-forge lod --sample-volume 200 --eluate-volume 50 --reaction-input 5

    \b
    Example (16S rRNA control, 1 mL produced water, 7 copies per cell):
        lamp-forge lod --sample-volume 1000 --eluate-volume 100 \
          --reaction-input 5 --copies-per-cell 7
    """
    from lamp_forge.lod import ExtractionParams, lod_table, write_lod_csv

    try:
        params = ExtractionParams(
            sample_volume_ul=sample_volume_ul,
            extraction_efficiency=extraction_efficiency,
            eluate_volume_ul=eluate_volume_ul,
            reaction_input_ul=reaction_input_ul,
        )
    except ValueError as exc:
        click.secho(f"Parameter error: {exc}", fg="red", err=True)
        import sys

        sys.exit(2)

    if copies_per_cell < 1:
        click.secho("--copies-per-cell must be >= 1.", fg="red", err=True)
        sys.exit(2)

    estimates = lod_table(params, tuple(probabilities), copies_per_cell=copies_per_cell)

    click.echo(
        f"Extraction chain: {sample_volume_ul:.0f} uL sample, "
        f"{extraction_efficiency * 100:.0f}% efficiency, "
        f"{eluate_volume_ul:.0f} uL eluate, "
        f"{reaction_input_ul:.1f} uL to reaction"
    )
    click.echo(
        f"Effective sample volume per reaction: {params.copies_per_rxn_per_copy_per_ul:.2f} uL"
    )
    if copies_per_cell > 1:
        click.echo(f"Gene copies per cell: {copies_per_cell} (LOD also shown in cells/mL)")
    click.echo("")

    if copies_per_cell > 1:
        header = (
            f"{'P(detect)':>12}  {'lambda (copies/rxn)':>20}  "
            f"{'LOD (copies/mL)':>18}  {'LOD (cells/mL)':>16}"
        )
    else:
        header = f"{'P(detect)':>12}  {'lambda (copies/rxn)':>20}  {'LOD (copies/mL)':>18}"
    click.echo(header)
    click.echo("-" * len(header))
    for e in estimates:
        if copies_per_cell > 1:
            click.echo(
                f"{e.detection_probability:>12.3f}  "
                f"{e.lod_copies_per_reaction:>20.3f}  "
                f"{e.lod_copies_per_ml:>18.1f}  "
                f"{e.lod_cells_per_ml:>16.1f}"
            )
        else:
            click.echo(
                f"{e.detection_probability:>12.3f}  "
                f"{e.lod_copies_per_reaction:>20.3f}  "
                f"{e.lod_copies_per_ml:>18.1f}"
            )

    if out_csv is not None:
        write_lod_csv(estimates, out_csv)
        click.echo(f"\nLOD table written to {out_csv}")


@cli.command(name="ttp")
@click.option(
    "--preset",
    type=click.Choice(["dna-lamp", "rt-lamp", "fast-lamp"], case_sensitive=False),
    default="dna-lamp",
    show_default=True,
    help=(
        "Chemistry preset: 'dna-lamp' (Bst 2.0 WarmStart, 65 degC; default for bacterial "
        "and DNA-virus assays), 'rt-lamp' (RTx + Bst 2.0, 63-65 degC; for RNA-target "
        "assays such as PRRSV, FMDV, avian influenza), 'fast-lamp' (high enzyme loading, "
        "early fluorescence). Overridden by --ttp-one-copy / --slope if both are given."
    ),
)
@click.option(
    "--ttp-one-copy",
    "ttp_one_copy",
    type=float,
    default=None,
    help=(
        "TTP (minutes) at 1 copy per reaction. Overrides the preset value. "
        "Typical range: 45-65 min depending on enzyme and temperature."
    ),
)
@click.option(
    "--slope",
    "slope",
    type=float,
    default=None,
    help=(
        "TTP reduction (minutes) per 10-fold increase in copy count. "
        "Overrides the preset value. Typical range: 4-8 min/decade."
    ),
)
@click.option(
    "--device-window",
    "device_window",
    type=float,
    default=60.0,
    show_default=True,
    help=(
        "Device run window in minutes. Results outside this window are flagged "
        "as 'FAIL'. Default 60 min (BioVind BioID 30-60 min run)."
    ),
)
@click.option(
    "--copies-min",
    "copies_min",
    type=float,
    default=1.0,
    show_default=True,
    help="Minimum copies per reaction for the table (default 1).",
)
@click.option(
    "--copies-max",
    "copies_max",
    type=float,
    default=1e6,
    show_default=True,
    help="Maximum copies per reaction for the table (default 1e6).",
)
@click.option(
    "--n-points",
    "n_points",
    type=int,
    default=13,
    show_default=True,
    help=(
        "Number of log-spaced evaluation points across the copy-count range. "
        "Default 13 gives a half-decade step over 1-10^6."
    ),
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the TTP table to a CSV file.",
)
def ttp(
    preset: str,
    ttp_one_copy: float | None,
    slope: float | None,
    device_window: float,
    copies_min: float,
    copies_max: float,
    n_points: int,
    out_csv: Path | None,
) -> None:
    r"""Estimate time-to-positive (TTP) across a copy-count range.

    Predicts when a LAMP or RT-LAMP reaction will turn positive for a given
    initial copy count using the empirically validated linear log10 model::

        TTP(N) = ttp_one_copy_min - slope * log10(N)

    Results are flagged PASS/FAIL against the device window (default 60 min
    for BioVind BioID). Use this before ordering primers to confirm the assay
    will read out within the device run time for the expected sample load.

    \b
    Examples:
        # DNA-LAMP SRB assay, default 60-min window:
        lamp-forge ttp --preset dna-lamp

        # One-step RT-LAMP (PRRSV), 45-min window:
        lamp-forge ttp --preset rt-lamp --device-window 45

        # Custom parameters, export to CSV:
        lamp-forge ttp --ttp-one-copy 58 --slope 5.5 --out-csv results/ttp.csv
    """
    from lamp_forge.ttp import TtpParams, TtpPreset, ttp_table, write_ttp_csv

    preset_enum = TtpPreset(preset.lower())
    base = TtpParams.from_preset(preset_enum, device_window_min=device_window)

    ttp_one_val = ttp_one_copy if ttp_one_copy is not None else base.ttp_one_copy_min
    slope_val = slope if slope is not None else base.slope_min_per_decade

    try:
        params = TtpParams(
            ttp_one_copy_min=ttp_one_val,
            slope_min_per_decade=slope_val,
            device_window_min=device_window,
            min_ttp_min=base.min_ttp_min,
        )
    except ValueError as exc:
        click.secho(f"Parameter error: {exc}", fg="red", err=True)
        sys.exit(2)

    try:
        table = ttp_table(params, copies_range=(copies_min, copies_max), n_points=n_points)
    except ValueError as exc:
        click.secho(f"Range error: {exc}", fg="red", err=True)
        sys.exit(2)

    click.echo(
        f"TTP model: preset={preset.lower()}, "
        f"TTP@1cp={params.ttp_one_copy_min:.1f} min, "
        f"slope={params.slope_min_per_decade:.1f} min/decade, "
        f"min_TTP={params.min_ttp_min:.1f} min"
    )
    click.echo(f"Device window: {params.device_window_min:.0f} min")
    click.echo("")

    header = f"{'Copies/rxn':>14}  {'TTP (min)':>10}  {'vs window':>10}"
    click.echo(header)
    click.echo("-" * len(header))
    for pt in table.points:
        status = "PASS" if pt.in_device_window else "FAIL"
        click.echo(f"{pt.copies_per_reaction:>14.1f}  {pt.ttp_minutes:>10.1f}  {status:>10}")

    if table.min_detectable_copies is not None:
        click.echo(
            f"\nMin detectable: {table.min_detectable_copies:.1f} copies/rxn "
            f"=> TTP within {params.device_window_min:.0f} min"
        )
    else:
        click.echo(
            f"\nWARNING: no copy count in range [{copies_min:.0f}, "
            f"{copies_max:.0f}] yields TTP within {params.device_window_min:.0f} min."
        )

    if out_csv is not None:
        write_ttp_csv(table, out_csv)
        click.echo(f"\nTTP table written to {out_csv}")


@cli.command(name="export")
@click.option(
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="primer_sets.json produced by 'lamp-forge run'.",
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["idt", "twist", "eurofins"], case_sensitive=False),
    default="idt",
    show_default=True,
    help="Vendor order-sheet format.",
)
@click.option(
    "--out",
    "out_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Output CSV path. Defaults to <input_stem>_<format>_order.csv.",
)
@click.option(
    "--top-n",
    default=1,
    type=int,
    show_default=True,
    help="Number of top-ranked primer sets to export (1 = best set only).",
)
@click.option(
    "--target-label",
    default=None,
    help=(
        "Short target label prepended to primer names (e.g. 'dsrB_SRB'). "
        "Makes order-sheet rows self-documenting when sharing with a vendor or lab."
    ),
)
@click.option(
    "--scale",
    default=None,
    help=(
        "Override synthesis scale for all primers "
        "(e.g. '25nm', '100nm', '250nm'). "
        "Default: role-specific (25nm for F3/B3/LF/LB, 100nm for FIP/BIP)."
    ),
)
@click.option(
    "--purification",
    default=None,
    help=(
        "Override purification for all primers "
        "(e.g. 'STD', 'HPLC', 'PAGE'). "
        "Default: role-specific (STD for outer/loop primers, HPLC for FIP/BIP)."
    ),
)
def export(
    input_path: Path,
    fmt: str,
    out_path: Path | None,
    top_n: int,
    target_label: str | None,
    scale: str | None,
    purification: str | None,
) -> None:
    r"""Export designed primers to an IDT, Twist, or Eurofins vendor order CSV.

    Reads the primer_sets.json produced by 'lamp-forge run' and writes a
    spreadsheet ready for direct upload to the IDT SDS, Twist oligo-order,
    or Eurofins Genomics portal.

    Scale and purification default to role-specific values:
    FIP/BIP get 100nm/HPLC (chimeric oligos benefit from full HPLC
    purification); F3/B3/LF/LB get 25nm/STD.  For Eurofins output,
    IDT scale strings are automatically converted (25nm -> 0.04umol,
    100nm -> 0.2umol, STD -> Desalt).

    \b
    Example (best set in IDT format):
        lamp-forge export \\
          --input results/primer_sets.json \\
          --format idt \\
          --target-label dsrB_SRB \\
          --out orders/dsrB_SRB_idt.csv

    \b
    Example (top 3 sets in Twist format):
        lamp-forge export \\
          --input results/primer_sets.json \\
          --format twist \\
          --top-n 3 \\
          --out orders/asfv_twist.csv

    \b
    Example (best set in Eurofins format):
        lamp-forge export \\
          --input results/primer_sets.json \\
          --format eurofins \\
          --target-label dsrB_SRB \\
          --out orders/dsrB_SRB_eurofins.csv
    """
    import json

    from lamp_forge.vendor_export import VendorFormat, rows_from_json_data, write_vendor_csv

    with input_path.open(encoding="utf-8") as fh:
        data: dict[str, object] = json.load(fh)

    sets_data = data.get("primer_sets", [])
    if not isinstance(sets_data, list) or not sets_data:
        click.secho("No primer sets found in the input file.", fg="red", err=True)
        sys.exit(1)

    vendor_fmt = VendorFormat(fmt.lower())
    if out_path is None:
        out_path = input_path.parent / f"{input_path.stem}_{fmt.lower()}_order.csv"

    rows = rows_from_json_data(
        sets_data,
        target_label=target_label,
        top_n=top_n,
        scale_override=scale,
        purification_override=purification,
    )
    if not rows:
        click.secho("No primers extracted — check the input file.", fg="red", err=True)
        sys.exit(1)

    write_vendor_csv(rows, out_path, vendor_fmt)

    n_sets = min(top_n, len(sets_data))
    click.secho(
        f"Wrote {len(rows)} primers ({n_sets} set(s)) -> {out_path}",
        fg="green",
    )
    click.echo(f"Format: {vendor_fmt.value.upper()} bulk order")


@cli.command(name="rt-check")
@click.option(
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="primer_sets.json produced by 'lamp-forge run'.",
)
@click.option(
    "--na-type",
    "na_type",
    type=click.Choice(["rna", "dna"], case_sensitive=False),
    default="rna",
    show_default=True,
    help=(
        "Target nucleic acid type. "
        "Use 'rna' for RNA virus targets (one-step RT-LAMP); "
        "'dna' confirms RT is not required."
    ),
)
@click.option(
    "--rt-min-tm",
    "rt_min_tm",
    type=float,
    default=63.0,
    show_default=True,
    help=(
        "Minimum Tm (degC) for core primers (F3/B3/FIP/BIP) in one-step RT-LAMP. "
        "Primers below this value may reduce reverse-transcriptase co-activity."
    ),
)
@click.option(
    "--top-n",
    "top_n",
    default=None,
    type=int,
    help="Assess only the top-N ranked sets (default: all).",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write full results table to a CSV file.",
)
def rt_check(
    input_path: Path,
    na_type: str,
    rt_min_tm: float,
    top_n: int | None,
    out_csv: Path | None,
) -> None:
    r"""Check whether primer sets are optimised for one-step RT-LAMP.

    RNA-target LAMP assays (PRRSV, avian influenza, FMDV, SARS-CoV-2, etc.)
    require a reverse transcription step before amplification.  One-step
    RT-LAMP combines the RT and LAMP reactions at 63-65 degC; primers designed
    for the broad 60-65 degC window may have Tms that are suboptimal when the
    RT enzyme must remain co-active with Bst polymerase.

    Reads a primer_sets.json produced by 'lamp-forge run' and reports, for
    each set, how many primers fall within the RT-LAMP optimal Tm range and
    whether the set is ready to order for one-step RT-LAMP.

    Sets marked NOT OPTIMIZED should be re-designed with a tighter Tm window
    (primer.tm_min >= 63.0) before ordering for RNA-virus applications.

    \b
    Example (PRRSV ORF7 assay targeting RNA):
        lamp-forge rt-check \\
          --input results/prrsv_orf7/primer_sets.json \\
          --na-type rna \\
          --out-csv results/prrsv_orf7/rt_check.csv

    \b
    Example (ASFV DNA assay, confirm RT not needed):
        lamp-forge rt-check \\
          --input results/asfv_p72/primer_sets.json \\
          --na-type dna
    """
    import json

    from lamp_forge.rt_lamp import (
        RtLampParams,
        TargetNucleicAcid,
        check_primer_sets_for_rt_lamp,
        write_rt_check_csv,
    )

    with input_path.open(encoding="utf-8") as fh:
        data: dict[str, object] = json.load(fh)

    sets_data = data.get("primer_sets", [])
    if not isinstance(sets_data, list) or not sets_data:
        click.secho("No primer sets found in the input file.", fg="red", err=True)
        sys.exit(1)

    target_na = TargetNucleicAcid(na_type.lower())
    try:
        params = RtLampParams(rt_min_tm=rt_min_tm)
    except ValueError as exc:
        click.secho(f"Parameter error: {exc}", fg="red", err=True)
        sys.exit(2)

    results = check_primer_sets_for_rt_lamp(
        sets_data,
        target_na=target_na,
        params=params,
        top_n=top_n,
    )

    click.echo(f"RT-LAMP compatibility check -- target: {na_type.upper()}")
    if target_na is TargetNucleicAcid.RNA:
        click.echo(
            f"Parameters: core primers {params.rt_min_tm:.1f}-{params.rt_max_tm:.1f} degC, "
            f"loop primers >= {params.loop_min_tm:.1f} degC"
        )
    click.echo("")

    col_id = max(len(r.set_id) for r in results) if results else 10
    col_id = max(col_id, 8)
    header = f"{'Set ID':<{col_id}}  {'N':>2}  {'In-range':>8}  {'Core-low':>8}  Status"
    click.echo(header)
    click.echo("-" * len(header))

    for r in results:
        status = "OK" if r.is_rt_optimized else "NOT OPTIMIZED"
        click.echo(
            f"{r.set_id:<{col_id}}  {r.n_primers:>2}  {r.n_in_rt_range:>8}  "
            f"{r.n_core_below_rt_min:>8}  {status}"
        )
        for w in r.warnings:
            click.echo(f"  {'':>{col_id}}  {w}")

    n_ok = sum(1 for r in results if r.is_rt_optimized)
    click.echo("")
    click.echo(f"Summary: {n_ok} of {len(results)} set(s) RT-LAMP optimized.")
    if target_na is TargetNucleicAcid.RNA and n_ok < len(results):
        click.secho(
            "Tip: tighten primer.tm_min to >= "
            f"{params.rt_min_tm:.1f} in your config and re-run for suboptimal sets.",
            fg="yellow",
        )

    if out_csv is not None:
        write_rt_check_csv(results, out_csv)
        click.echo(f"Results written to {out_csv}")


@cli.command(name="pool")
@click.option(
    "--panel",
    "panel_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="panel.json produced by 'lamp-forge panel'.",
)
@click.option(
    "--stock-conc",
    "stock_conc_um",
    type=float,
    default=100.0,
    show_default=True,
    help=(
        "Concentration of each synthesised primer stock tube (uM). "
        "IDT/Twist standard is 100 uM; request 200 uM for 3+ target panels."
    ),
)
@click.option(
    "--total-volume",
    "total_pool_volume_ul",
    type=float,
    default=500.0,
    show_default=True,
    help="Target total volume of the pooled 10x working mix (uL).",
)
@click.option(
    "--out",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Output CSV path. Defaults to <panel_dir>/pool_sheet.csv.",
)
def pool(
    panel_path: Path,
    stock_conc_um: float,
    total_pool_volume_ul: float,
    out_csv: Path | None,
) -> None:
    r"""Calculate primer volumes for a multiplex pooled working stock.

    Reads the panel.json produced by 'lamp-forge panel' and emits a
    pipetting sheet: for each primer in the selected panel, how many uL
    to transfer from its synthesis stock tube into the shared working mix.

    Standard LAMP 10x working-stock concentrations applied per target:

    \b
        FIP / BIP  : 16 uM  ->  1.6 uM final in the 1x reaction
        LF  / LB   :  4 uM  ->  0.4 uM final in the 1x reaction
        F3  / B3   :  2 uM  ->  0.2 uM final in the 1x reaction

    The pool is used at 1/10 volume per reaction (e.g. 2.5 uL in a 25 uL
    LAMP reaction).  Each target contributes the same per-reaction primer
    copy number as it would in a single-target LAMP.

    \b
    Example (oilfield souring 3-target panel, 200 uM stocks):
        lamp-forge pool \\
          --panel results/souring_panel/panel.json \\
          --stock-conc 200 \\
          --total-volume 500 \\
          --out results/souring_panel/pool_sheet.csv

    Note: the minimum required stock concentration equals the sum of all
    primer working concentrations across all targets (44 uM x N targets for
    complete 6-primer sets). Use --stock-conc to match your vendor resuspension.
    """
    import json

    from lamp_forge.pool import PoolingParams, build_pool_plan, write_pool_csv

    with panel_path.open(encoding="utf-8") as fh:
        panel_data: dict[str, object] = json.load(fh)

    selection_raw = panel_data.get("selection", [])
    if not isinstance(selection_raw, list) or not selection_raw:
        click.secho("No selection found in panel.json.", fg="red", err=True)
        sys.exit(1)

    try:
        params = PoolingParams(
            stock_conc_um=stock_conc_um,
            total_pool_volume_ul=total_pool_volume_ul,
        )
    except ValueError as exc:
        click.secho(f"Parameter error: {exc}", fg="red", err=True)
        sys.exit(2)

    try:
        plan = build_pool_plan(selection_raw, params)
    except ValueError as exc:
        click.secho(f"Pool calculation error: {exc}", fg="red", err=True)
        sys.exit(1)

    if out_csv is None:
        out_csv = panel_path.parent / "pool_sheet.csv"

    write_pool_csv(plan, out_csv)

    click.secho(
        f"Pool: {plan.n_targets} target(s), {plan.n_primers} primers, "
        f"{plan.total_pool_volume_ul:.0f} uL total "
        f"({plan.water_volume_ul:.1f} uL water).",
        fg="green",
    )
    click.echo(f"Pooling sheet: {out_csv}")
    click.echo("")

    col_t = max((len(e.target_label) for e in plan.entries), default=10)
    col_t = max(col_t, 10)
    col_n = max((len(e.primer_name) for e in plan.entries), default=12)
    col_n = max(col_n, 12)
    header = (
        f"{'Target':<{col_t}}  {'Role':<5}  {'Name':<{col_n}}  "
        f"{'Stock conc':>10}  {'Pool conc':>10}  {'Vol (uL)':>9}"
    )
    click.echo(header)
    click.echo("-" * len(header))
    for e in plan.entries:
        click.echo(
            f"{e.target_label:<{col_t}}  {e.primer_role:<5}  "
            f"{e.primer_name:<{col_n}}  "
            f"{e.stock_conc_um:>9.0f}u  {e.target_conc_um:>9.1f}u  "
            f"{e.vol_stock_ul:>9.3f}"
        )
    click.echo("-" * len(header))
    click.echo(f"{'Nuclease-free water':<{col_t + 5 + col_n + 25}}  {plan.water_volume_ul:>9.3f}")
    click.echo(f"{'TOTAL':<{col_t + 5 + col_n + 25}}  {plan.total_pool_volume_ul:>9.3f}")


@cli.command(name="validate")
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to the YAML pipeline config to validate.",
)
@click.option(
    "--no-check-dirs",
    "skip_dir_check",
    is_flag=True,
    default=False,
    help="Skip filesystem checks (off-target directory existence). Useful in CI.",
)
def validate(config_path: Path, skip_dir_check: bool) -> None:
    r"""Validate a LAMP-Forge YAML config without running the pipeline.

    Loads and validates the config (the same checks run by ``lamp-forge run``),
    prints a structured summary, and reports soft warnings for parameter choices
    that are valid but likely suboptimal.

    Exits 0 if the config is valid (warnings do not change the exit code).
    Exits 2 on a hard config error.

    \b
    Example:
        lamp-forge validate --config config/srb_dsrB.yaml
    """
    from lamp_forge.config import ConfigError, check_config_warnings, load_config

    try:
        config = load_config(config_path)
    except ConfigError as e:
        click.secho(f"Config error: {e}", fg="red", err=True)
        sys.exit(2)

    # Source description
    if config.taxon_id and config.accessions:
        source = f"taxon_id {config.taxon_id} + {len(config.accessions)} accession(s)"
    elif config.taxon_id:
        source = f"taxon_id {config.taxon_id}"
    else:
        source = f"{len(config.accessions)} explicit accession(s)"
    gene_label = config.gene if config.gene else "(whole genome / accession)"

    click.secho(f"Config OK -- {config.target_name}", fg="green")
    click.echo(f"  Target:        {config.target_name}")
    click.echo(
        f"  Source:        {source}, gene={gene_label}, max {config.max_sequences} sequences"
    )
    click.echo(
        f"  Off-target:    {config.off_target_dir} "
        f"(identity >={config.min_identity_threshold:.0%}, "
        f"coverage >={config.min_coverage_threshold:.0%})"
    )
    click.echo(
        f"  Conservation:  window={config.window_size}, "
        f"entropy<={config.entropy_threshold} bits, "
        f"min region {config.min_region_length}bp"
    )
    click.echo(
        f"  Primers:       Tm {config.tm_min:.1f}-{config.tm_max:.1f} degC "
        f"(tolerance {config.tm_match_tolerance:.1f}), "
        f"GC {config.gc_min:.0f}-{config.gc_max:.0f}%, "
        f"amplicon {config.f2_b2_min}-{config.f2_b2_max}bp"
    )
    click.echo(
        f"  Output:        {config.output_dir}, "
        f"top {config.top_n} sets, "
        f"HTML={'yes' if config.generate_html else 'no'}, "
        f"CSV={'yes' if config.generate_csv else 'no'}"
    )

    soft_warnings = check_config_warnings(config, check_dirs=not skip_dir_check)
    if soft_warnings:
        click.echo("")
        click.secho(f"Warnings ({len(soft_warnings)}):", fg="yellow")
        for w in soft_warnings:
            click.secho(f"  {w}", fg="yellow")
    else:
        click.echo("  No warnings.")


@cli.command(name="scaffold")
@click.option(
    "--target-name",
    "target_name",
    required=True,
    help=(
        "Short identifier for the target "
        "(e.g. 'prrsv_orf7', 'srb_dsrB'). Used in target.name and output paths."
    ),
)
@click.option(
    "--vertical",
    type=click.Choice(["oil-gas", "farm", "poc", "generic"], case_sensitive=False),
    default="generic",
    show_default=True,
    help=(
        "Deployment vertical: 'oil-gas' (MIC/souring), 'farm' (animal biosecurity), "
        "'poc' (human point-of-care), 'generic' (neutral defaults)."
    ),
)
@click.option(
    "--na-type",
    "na_type",
    type=click.Choice(["dna", "rna"], case_sensitive=False),
    default="dna",
    show_default=True,
    help=(
        "Target nucleic-acid type. 'rna' raises the primer Tm floor to 63 degC "
        "for one-step RT-LAMP co-activity."
    ),
)
@click.option(
    "--taxon-id",
    "taxon_id",
    type=int,
    default=None,
    help="NCBI taxonomy ID for the target organism (e.g. 28344 for PRRSV).",
)
@click.option(
    "--gene",
    default=None,
    help="Target gene name as annotated in NCBI records (e.g. 'ORF7', 'dsrB').",
)
@click.option(
    "--email",
    default="your-email@institution.edu",
    show_default=True,
    help="NCBI courtesy email written into the config.",
)
@click.option(
    "--max-sequences",
    "max_sequences",
    type=int,
    default=None,
    help="Override the vertical default maximum sequence count.",
)
@click.option(
    "--out",
    "out_path",
    type=click.Path(path_type=Path),
    default=None,
    help=(
        "Output YAML path. Omit to print to stdout "
        "(useful for piping: lamp-forge scaffold ... | less)."
    ),
)
def scaffold(
    target_name: str,
    vertical: str,
    na_type: str,
    taxon_id: int | None,
    gene: str | None,
    email: str,
    max_sequences: int | None,
    out_path: Path | None,
) -> None:
    r"""Generate a starter YAML config with vertical-appropriate defaults.

    Prints a ready-to-edit LAMP-Forge config pre-populated with
    parameter choices tuned to the selected deployment vertical and
    target nucleic-acid type.  Edit the placeholder email and taxon_id,
    drop off-target FASTAs into input/off_targets/, then validate and run.

    \b
    Example -- farm RNA virus (PRRSV):
        lamp-forge scaffold \\
          --target-name prrsv_orf7 \\
          --vertical farm \\
          --na-type rna \\
          --taxon-id 28344 \\
          --gene ORF7 \\
          --out config/prrsv_orf7.yaml

    \b
    Example -- oilfield MIC (SRB dsrB):
        lamp-forge scaffold \\
          --target-name srb_dsrB \\
          --vertical oil-gas \\
          --taxon-id 872 \\
          --gene dsrB \\
          --out config/srb_dsrB_new.yaml
    """
    from lamp_forge.scaffold import NaType as ScNaType
    from lamp_forge.scaffold import ScaffoldParams, scaffold_yaml, write_scaffold
    from lamp_forge.scaffold import Vertical as ScVertical

    params = ScaffoldParams(
        target_name=target_name,
        vertical=ScVertical(vertical.lower()),
        na_type=ScNaType(na_type.lower()),
        taxon_id=taxon_id,
        gene=gene,
        email=email,
        max_sequences=max_sequences,
    )

    if out_path is None:
        click.echo(scaffold_yaml(params), nl=False)
    else:
        write_scaffold(params, out_path)
        click.secho(f"Scaffold written to {out_path}", fg="green")
        click.echo("Next: edit target.email and target.taxon_id, then run:")
        click.echo(f"  lamp-forge validate --config {out_path}")


@cli.command(name="panel-export")
@click.option(
    "--panel",
    "panel_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="panel.json produced by 'lamp-forge panel'.",
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["idt", "twist", "eurofins"], case_sensitive=False),
    default="idt",
    show_default=True,
    help="Vendor order-sheet format.",
)
@click.option(
    "--out",
    "out_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Output CSV path. Defaults to <panel_dir>/panel_<format>_order.csv.",
)
@click.option(
    "--scale",
    default=None,
    help=(
        "Override synthesis scale for all primers "
        "(e.g. '25nm', '100nm', '250nm'). "
        "Default: role-specific (25nm for F3/B3/LF/LB, 100nm for FIP/BIP)."
    ),
)
@click.option(
    "--purification",
    default=None,
    help=(
        "Override purification for all primers "
        "(e.g. 'STD', 'HPLC', 'PAGE'). "
        "Default: role-specific (STD for outer/loop primers, HPLC for FIP/BIP)."
    ),
)
def panel_export(
    panel_path: Path,
    fmt: str,
    out_path: Path | None,
    scale: str | None,
    purification: str | None,
) -> None:
    r"""Export all primers from a multiplex panel to a single vendor order CSV.

    Reads the panel.json produced by 'lamp-forge panel' and exports the
    selected primer set for every target into one IDT, Twist, or Eurofins
    order sheet -- removing the need to run 'lamp-forge export' separately
    for each target.

    Each primer name is prefixed with its target label so rows are
    self-documenting in the vendor portal (e.g. SRB_region_01_set_001_FIP).

    Note: Tm values are not stored in panel.json; the Tm column in Twist and
    Eurofins output will show 0.0.  For per-primer Tm, use 'lamp-forge export
    --input primer_sets.json' for each target individually.

    \b
    Example (5-channel oilfield MIC panel, IDT format):
        lamp-forge panel-export \
          --panel results/mic_5plex_panel/panel.json \
          --format idt \
          --out orders/mic_5plex_idt_order.csv

    \b
    Example (farm biosecurity 4-plex, Twist format):
        lamp-forge panel-export \
          --panel results/farm_biosecurity_4plex/panel.json \
          --format twist \
          --out orders/farm_biosecurity_twist.csv

    \b
    Example (oilfield souring 3-plex, Eurofins format):
        lamp-forge panel-export \
          --panel results/souring_panel/panel.json \
          --format eurofins \
          --out orders/souring_eurofins_order.csv
    """
    import json

    from lamp_forge.vendor_export import VendorFormat, rows_from_panel_json, write_vendor_csv

    with panel_path.open(encoding="utf-8") as fh:
        panel_data: dict[str, object] = json.load(fh)

    selection_raw = panel_data.get("selection", [])
    if not isinstance(selection_raw, list) or not selection_raw:
        click.secho("No selection found in panel.json.", fg="red", err=True)
        sys.exit(1)

    vendor_fmt = VendorFormat(fmt.lower())
    if out_path is None:
        out_path = panel_path.parent / f"panel_{fmt.lower()}_order.csv"

    rows = rows_from_panel_json(
        selection_raw,
        scale_override=scale,
        purification_override=purification,
    )
    if not rows:
        click.secho("No primers extracted from panel.json.", fg="red", err=True)
        sys.exit(1)

    write_vendor_csv(rows, out_path, vendor_fmt)

    n_targets = len(selection_raw)
    click.secho(
        f"Wrote {len(rows)} primers ({n_targets} target(s)) -> {out_path}",
        fg="green",
    )
    click.echo(f"Format: {vendor_fmt.value.upper()} bulk order")


@cli.command(name="mic-risk")
@click.option("--srb", "srb", is_flag=True, default=False, help="SRB (dsrB) positive.")
@click.option("--mcr", "mcr", is_flag=True, default=False, help="Methanogens (mcrA) positive.")
@click.option("--irb", "irb", is_flag=True, default=False, help="IRB (omcA) positive.")
@click.option(
    "--apb", "apb", is_flag=True, default=False, help="APB/homoacetogens (fthfs) positive."
)
@click.option("--nrb", "nrb", is_flag=True, default=False, help="NRB (narG) positive.")
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read guild flags from a JSON file instead of individual flags. "
        "Expected keys: srb, mcr, irb, apb, nrb (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def mic_risk(
    srb: bool,
    mcr: bool,
    irb: bool,
    apb: bool,
    nrb: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a five-guild oilfield LAMP panel as a MIC risk assessment.

    Takes the positivity flags for the five functional guilds monitored by a
    BioVind-style oilfield corrosion panel (SRB, methanogens, IRB, APB, NRB)
    and returns a structured risk assessment with risk level, numeric score,
    interpretation, and recommended action.

    Guild markers:
    \b
        --srb   Sulfate-reducing bacteria       dsrB   (Recipe 6)
        --mcr   Methanogens                     mcrA   (Recipe 10)
        --irb   Iron-reducing bacteria          omcA   (Recipe 15)
        --apb   Acid-producing bacteria         fthfs  (Recipe 16)
        --nrb   Nitrate-reducing bacteria       narG   (Recipe 17)

    \b
    Example -- SRB + IRB co-detected (FeS scale risk):
        lamp-forge mic-risk --srb --irb

    \b
    Example -- full five-guild panel from a JSON file:
        lamp-forge mic-risk --input-json results/guild_flags.json \
          --out-json results/mic_assessment.json

    \b
    Example -- all four corrosion guilds active:
        lamp-forge mic-risk --srb --mcr --irb --apb
    """
    import json as json_mod

    from lamp_forge.mic_risk import (
        GuildFlags,
        MICRiskLevel,
        assess_mic_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = GuildFlags(srb=srb, mcr=mcr, irb=irb, apb=apb, nrb=nrb)

    assessment = assess_mic_risk(flags)

    # --- Guild table ----------------------------------------------------------
    click.echo("Oilfield MIC guild panel results:")
    guild_rows = [
        ("SRB", "dsrB", assessment.flags.srb),
        ("MCR", "mcrA", assessment.flags.mcr),
        ("IRB", "omcA", assessment.flags.irb),
        ("APB", "fthfs", assessment.flags.apb),
        ("NRB", "narG", assessment.flags.nrb),
    ]
    for label, gene, positive in guild_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<5} ({gene:<5})  [{symbol}]")

    click.echo("")

    # --- Risk level (colour-coded) -------------------------------------------
    level_color = {
        MICRiskLevel.CRITICAL: "red",
        MICRiskLevel.HIGH: "red",
        MICRiskLevel.MODERATE: "yellow",
        MICRiskLevel.LOW: "green",
        MICRiskLevel.MINIMAL: "green",
    }
    color = level_color[assessment.risk_level]
    click.secho(
        f"Risk level : {assessment.risk_level.value}  (score {assessment.risk_score}/100)",
        fg=color,
        bold=(assessment.risk_level in (MICRiskLevel.CRITICAL, MICRiskLevel.HIGH)),
    )

    if assessment.nrb_suppression_active:
        click.secho("  Note: NRB+ / SRB- -- nitrate injection appears effective.", fg="green")
    if assessment.treatment_underdosed:
        click.secho(
            "  Warning: NRB and SRB both detected -- treatment may be under-dosed.",
            fg="yellow",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if assessment.corrosion_drivers:
        click.echo("")
        click.echo("Active corrosion pathways:")
        for driver in assessment.corrosion_drivers:
            click.echo(f"  - {driver}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="farm-risk")
@click.option("--asfv", "asfv", is_flag=True, default=False, help="ASFV (B646L / p72) positive.")
@click.option("--fmdv", "fmdv", is_flag=True, default=False, help="FMDV (3Dpol) positive.")
@click.option(
    "--aiv", "aiv", is_flag=True, default=False, help="Avian influenza A (M gene) positive."
)
@click.option(
    "--ndv", "ndv", is_flag=True, default=False, help="Newcastle disease virus (M gene) positive."
)
@click.option("--prrsv", "prrsv", is_flag=True, default=False, help="PRRSV (ORF7) positive.")
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: asfv, fmdv, aiv, ndv, prrsv (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def farm_risk(
    asfv: bool,
    fmdv: bool,
    aiv: bool,
    ndv: bool,
    prrsv: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a farm-biosecurity LAMP panel as a structured alert.

    Takes the positivity flags for the five pathogens monitored by a
    BioVind-style farm-biosecurity panel and returns a structured alert
    with level, score, interpretation, and recommended action.

    WOAH-listed notifiable pathogens (ASFV, FMDV, AIV, NDV) trigger an
    immediate-report flag -- confirm at a national reference laboratory before
    depopulation or trade restrictions are enforced.

    Panel targets:
    \b
        --asfv    African swine fever virus        B646L / p72  (Recipe 7)
        --fmdv    Foot-and-mouth disease virus     3Dpol        (Recipe 14)
        --aiv     Avian influenza A (pan-IAV)      M gene       (Recipe 12)
        --ndv     Newcastle disease virus          M gene       (Recipe 18)
        --prrsv   PRRS virus                       ORF7         (Recipe 11)

    \b
    Example -- ASFV detected on a pig farm:
        lamp-forge farm-risk --asfv

    \b
    Example -- full panel from a JSON flags file:
        lamp-forge farm-risk --input-json results/farm_flags.json \
          --out-json results/farm_assessment.json

    \b
    Example -- PRRSV + NDV co-detected:
        lamp-forge farm-risk --prrsv --ndv
    """
    import json as json_mod

    from lamp_forge.farm_risk import (
        FarmAlertLevel,
        FarmPanelFlags,
        assess_farm_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = FarmPanelFlags(asfv=asfv, fmdv=fmdv, aiv=aiv, ndv=ndv, prrsv=prrsv)

    assessment = assess_farm_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Farm-biosecurity LAMP panel results:")
    target_rows = [
        ("ASFV", "B646L", assessment.flags.asfv),
        ("FMDV", "3Dpol", assessment.flags.fmdv),
        ("AIV", "M gene", assessment.flags.aiv),
        ("NDV", "M gene", assessment.flags.ndv),
        ("PRRSV", "ORF7", assessment.flags.prrsv),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<6} ({gene:<6})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        FarmAlertLevel.CRITICAL: "red",
        FarmAlertLevel.HIGH: "red",
        FarmAlertLevel.MODERATE: "yellow",
        FarmAlertLevel.LOW: "yellow",
        FarmAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level in (FarmAlertLevel.CRITICAL, FarmAlertLevel.HIGH)),
    )

    if assessment.immediate_report_required:
        click.secho(
            "  IMMEDIATE REPORT REQUIRED: notify national veterinary authority.",
            fg="red",
            bold=True,
        )
        click.secho(
            f"  Notifiable detection(s): {', '.join(assessment.notifiable_targets)}",
            fg="red",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="poc-risk")
@click.option("--mtb", "mtb", is_flag=True, default=False, help="MTB (rpoB) positive.")
@click.option(
    "--cdiff", "cdiff", is_flag=True, default=False, help="C. difficile toxin B (tcdB) positive."
)
@click.option("--sars2", "sars2", is_flag=True, default=False, help="SARS-CoV-2 (N gene) positive.")
@click.option(
    "--flu-a", "flu_a", is_flag=True, default=False, help="Influenza A (M gene) positive."
)
@click.option(
    "--gas", "gas", is_flag=True, default=False, help="Group A Streptococcus (speB) positive."
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: mtb, cdiff, sars2, flu_a, gas (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def poc_risk(
    mtb: bool,
    cdiff: bool,
    sars2: bool,
    flu_a: bool,
    gas: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a human POC LAMP panel as a structured clinical alert.

    Takes the positivity flags for the five pathogens monitored by a
    BioVind-style human point-of-care LAMP panel and returns a structured
    alert with level, score, interpretation, clinical action flags, and
    recommended action.

    MTB (tuberculosis) is mandatorily notifiable in virtually all
    jurisdictions -- a positive triggers IMMEDIATE public-health reporting
    and airborne isolation.  Confirm at a national reference laboratory
    before initiating therapy.

    Panel targets:
    \b
        --mtb     Mycobacterium tuberculosis     rpoB    (Recipe 1)
        --cdiff   Clostridioides difficile tcdB  tcdB    (Recipe 22)
        --sars2   SARS-CoV-2                     N gene  (Recipe 2)
        --flu-a   Influenza A (pan-IAV)          M gene  (Recipe 12)
        --gas     Group A Streptococcus          speB    (Recipe 21)

    \b
    Example -- C. difficile detected in a hospital patient:
        lamp-forge poc-risk --cdiff

    \b
    Example -- SARS-CoV-2 + GAS pharyngitis co-detection:
        lamp-forge poc-risk --sars2 --gas

    \b
    Example -- full panel from a JSON flags file:
        lamp-forge poc-risk --input-json results/poc_flags.json \
          --out-json results/poc_assessment.json
    """
    import json as json_mod

    from lamp_forge.poc_risk import (
        PocAlertLevel,
        PocPanelFlags,
        assess_poc_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = PocPanelFlags(mtb=mtb, cdiff=cdiff, sars2=sars2, flu_a=flu_a, gas=gas)

    assessment = assess_poc_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Human POC LAMP panel results:")
    target_rows = [
        ("MTB", "rpoB", assessment.flags.mtb),
        ("CDIFF", "tcdB", assessment.flags.cdiff),
        ("SARS2", "N gene", assessment.flags.sars2),
        ("FLU_A", "M gene", assessment.flags.flu_a),
        ("GAS", "speB", assessment.flags.gas),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<6} ({gene:<6})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        PocAlertLevel.CRITICAL: "red",
        PocAlertLevel.HIGH: "red",
        PocAlertLevel.MODERATE: "yellow",
        PocAlertLevel.LOW: "yellow",
        PocAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level in (PocAlertLevel.CRITICAL, PocAlertLevel.HIGH)),
    )

    if assessment.immediate_report_required:
        click.secho(
            "  IMMEDIATE REPORT REQUIRED: notify national public health authority.",
            fg="red",
            bold=True,
        )
        click.secho(
            f"  Notifiable detection(s): {', '.join(assessment.notifiable_targets)}",
            fg="red",
        )

    if assessment.contact_precautions_required:
        click.secho(
            "  Contact precautions required"
            + (" + airborne isolation (MTB)" if assessment.flags.mtb else "")
            + ".",
            fg="yellow",
        )

    if assessment.antiviral_indicated:
        click.secho(
            "  Antiviral treatment indicated -- check treatment window.",
            fg="yellow",
        )

    if assessment.antibiotic_indicated:
        click.secho(
            "  Antibiotic therapy indicated (GAS pharyngitis).",
            fg="yellow",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="bov-risk")
@click.option("--brsv", "brsv", is_flag=True, default=False, help="BRSV (N gene) positive.")
@click.option("--bcov", "bcov", is_flag=True, default=False, help="BCoV (N gene) positive.")
@click.option(
    "--bvdv", "bvdv", is_flag=True, default=False, help="BVDV (5-UTR, types 1+2) positive."
)
@click.option(
    "--ibr", "ibr", is_flag=True, default=False, help="IBR / BoHV-1 (gB / UL27) positive."
)
@click.option(
    "--mhae",
    "mhae",
    is_flag=True,
    default=False,
    help="Mannheimia haemolytica (lktA leukotoxin A) positive.",
)
@click.option(
    "--mbovis",
    "mbovis",
    is_flag=True,
    default=False,
    help="Mycoplasma bovis (uvrC DNA repair gene) positive.",
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: brsv, bcov, bvdv, ibr, mhae, mbovis (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def bov_risk(
    brsv: bool,
    bcov: bool,
    bvdv: bool,
    ibr: bool,
    mhae: bool,
    mbovis: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a bovine respiratory LAMP panel as a structured BRD alert.

    Takes the positivity flags for the five pathogens monitored by a
    BioVind-style bovine respiratory panel and returns a structured alert with
    level, score, interpretation, and recommended action.

    Bovine respiratory disease (BRD / "shipping fever") is the leading cause
    of morbidity and mortality in beef feedlots worldwide.  This command
    surfaces the highest-mortality co-infection pattern (viral trigger +
    M. haemolytica bacterial co-infection) and flags IBR (BoHV-1) detections
    that may be subject to mandatory national control programme reporting.

    Panel targets:
    \b
        --brsv    Bovine respiratory syncytial virus    N gene
        --bcov    Bovine coronavirus                    N gene
        --bvdv    Bovine viral diarrhea virus (1+2)     5-UTR
        --ibr     Inf. bovine rhinotracheitis (BoHV-1)  gB / UL27
        --mhae    Mannheimia haemolytica                lktA (leukotoxin A)
        --mbovis  Mycoplasma bovis (6th channel)        uvrC (DNA repair)

    \b
    Example -- BRSV + M. haemolytica co-infection (highest-mortality BRD pattern):
        lamp-forge bov-risk --brsv --mhae

    \b
    Example -- BVDV detected alone (screen for PI animals):
        lamp-forge bov-risk --bvdv

    \b
    Example -- full panel result from a JSON flags file:
        lamp-forge bov-risk --input-json results/bov_flags.json \
          --out-json results/bov_assessment.json
    """
    import json as json_mod

    from lamp_forge.bov_risk import (
        BovAlertLevel,
        BovPanelFlags,
        assess_bov_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = BovPanelFlags(brsv=brsv, bcov=bcov, bvdv=bvdv, ibr=ibr, mhae=mhae, mbovis=mbovis)

    assessment = assess_bov_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Bovine respiratory LAMP panel results:")
    target_rows = [
        ("BRSV", "N gene", assessment.flags.brsv),
        ("BCOV", "N gene", assessment.flags.bcov),
        ("BVDV", "5-UTR", assessment.flags.bvdv),
        ("IBR", "gB", assessment.flags.ibr),
        ("MHAE", "lktA", assessment.flags.mhae),
        ("MBOVIS", "uvrC", assessment.flags.mbovis),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<6} ({gene:<6})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        BovAlertLevel.CRITICAL: "red",
        BovAlertLevel.HIGH: "red",
        BovAlertLevel.MODERATE: "yellow",
        BovAlertLevel.LOW: "yellow",
        BovAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level in (BovAlertLevel.CRITICAL, BovAlertLevel.HIGH)),
    )

    if assessment.bacterial_coinfection:
        click.secho(
            "  Viral-bacterial co-infection: virus + M. haemolytica -- highest-mortality BRD "
            "pattern. Initiate antimicrobial therapy immediately.",
            fg="red",
            bold=True,
        )

    if assessment.mbovis_pan_resistant:
        click.secho(
            "  M. bovis detected: DO NOT use beta-lactam or aminoglycoside protocols. "
            "Submit for susceptibility testing (MIC panel) before selecting antimicrobials.",
            fg="red",
            bold=True,
        )

    if assessment.ibr_mandatory_notification:
        click.secho(
            "  IBR (BoHV-1) detected: mandatory notification may be required in your "
            "jurisdiction (EU, UK, Scandinavia). Check national competent authority rules.",
            fg="yellow",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="version")
def version() -> None:
    """Print version and exit."""
    click.echo(f"lamp-forge {__version__}")


@cli.command(name="sampling-plan")
@click.option(
    "--prevalence",
    "prevalences",
    type=float,
    multiple=True,
    required=True,
    metavar="FLOAT",
    help=(
        "Minimum design prevalence to detect, as a fraction 0-1 "
        "(e.g. 0.05 for 5%). Repeat to generate a table across multiple "
        "prevalences (e.g. --prevalence 0.05 --prevalence 0.10)."
    ),
)
@click.option(
    "--confidence",
    "confidence",
    type=float,
    default=0.95,
    show_default=True,
    help=(
        "Required detection confidence as a fraction 0-1 "
        "(e.g. 0.95 for 95% probability of finding at least one "
        "positive sample if the true prevalence >= design prevalence)."
    ),
)
@click.option(
    "--sensitivity",
    "sensitivity",
    type=float,
    default=0.95,
    show_default=True,
    help=(
        "LAMP assay analytical sensitivity as a fraction 0-1 "
        "(probability of a positive result given target is truly present "
        "above the assay LOD). Default 0.95."
    ),
)
@click.option(
    "--population-size",
    "population_size",
    type=int,
    default=None,
    help=(
        "Total number of sampling units in the population (e.g. number of "
        "monitored wells, animals in a herd). Omit for effectively infinite "
        "populations. Supplying this value applies the hypergeometric "
        "finite-population correction, which reduces required n."
    ),
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the sampling plan table to a CSV file.",
)
def sampling_plan_cmd(
    prevalences: tuple[float, ...],
    confidence: float,
    sensitivity: float,
    population_size: int | None,
    out_csv: Path | None,
) -> None:
    r"""Calculate surveillance sample sizes for a LAMP panel deployment.

    Answers: "How many wells / animals must I test per site visit to be
    confident of detecting the target pathogen if it is present at or
    above a minimum prevalence?"

    Uses the standard freedom-from-disease binomial formula (infinite
    population) with an optional hypergeometric finite-population
    correction when --population-size is supplied.

    \b
    Formula (infinite population):
        n = ceil[ log(1-confidence) / log(1 - prevalence * sensitivity) ]

    \b
    Example -- oilfield SRB monitoring, 5%/10%/20% design prevalence:
        lamp-forge sampling-plan \\
          --prevalence 0.05 \\
          --prevalence 0.10 \\
          --prevalence 0.20 \\
          --confidence 0.95 \\
          --sensitivity 0.95

    \b
    Example -- farm ASFV screening, finite herd of 200 pigs:
        lamp-forge sampling-plan \\
          --prevalence 0.05 \\
          --prevalence 0.10 \\
          --population-size 200 \\
          --out-csv results/asfv_sampling_plan.csv

    \b
    Example -- POC clinic surveillance, detect if >= 1% of patients carry MTB:
        lamp-forge sampling-plan \\
          --prevalence 0.01 \\
          --confidence 0.95 \\
          --sensitivity 0.95
    """
    from lamp_forge.sampling import SamplingParams, sample_size, write_sampling_csv

    try:
        estimates = []
        for prev in prevalences:
            params = SamplingParams(
                target_prevalence=prev,
                confidence=confidence,
                test_sensitivity=sensitivity,
                population_size=population_size,
            )
            estimates.append(sample_size(params))
    except ValueError as exc:
        click.secho(f"Parameter error: {exc}", fg="red", err=True)
        sys.exit(2)

    pop_label = f"{population_size} units" if population_size is not None else "infinite"
    click.echo(
        f"Surveillance sampling plan: "
        f"confidence={confidence * 100:.0f}%, "
        f"sensitivity={sensitivity * 100:.0f}%, "
        f"population={pop_label}"
    )
    click.echo("")

    col_prev = 18
    col_n = 10
    col_ach = 24
    header = (
        f"{'Prevalence (%)':>{col_prev}}  {'n_samples':>{col_n}}  {'Achieved conf (%)':>{col_ach}}"
    )
    if population_size is not None:
        header += "  FPC"
    click.echo(header)
    click.echo("-" * len(header))

    for e in estimates:
        fpc_flag = "  yes" if e.finite_population_corrected else "  no"
        row = (
            f"{e.target_prevalence * 100:>{col_prev}.2f}  "
            f"{e.n_samples:>{col_n}}  "
            f"{e.achieved_confidence * 100:>{col_ach}.2f}"
        )
        if population_size is not None:
            row += fpc_flag
        click.echo(row)

    if out_csv is not None:
        write_sampling_csv(estimates, out_csv)
        click.echo(f"\nSampling plan written to {out_csv}")


@cli.command(name="preorder")
@click.option(
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="primer_sets.json produced by 'lamp-forge run'.",
)
@click.option(
    "--na-type",
    "na_type",
    type=click.Choice(["rna", "dna"], case_sensitive=False),
    default="dna",
    show_default=True,
    help=(
        "Target nucleic-acid type. 'rna' enables the RT-LAMP Tm check "
        "(use for PRRSV, FMDV, avian influenza, NDV, SARS-CoV-2, etc.)."
    ),
)
@click.option(
    "--sample-volume",
    "sample_volume_ul",
    type=float,
    required=True,
    help="Sample volume input to extraction (uL).",
)
@click.option(
    "--efficiency",
    "extraction_efficiency",
    type=float,
    default=0.50,
    show_default=True,
    help="Extraction efficiency as a fraction 0-1 (e.g. 0.50 for 50%).",
)
@click.option(
    "--eluate-volume",
    "eluate_volume_ul",
    type=float,
    required=True,
    help="Volume of extraction eluate (uL).",
)
@click.option(
    "--reaction-input",
    "reaction_input_ul",
    type=float,
    required=True,
    help="Eluate added to each LAMP reaction (uL).",
)
@click.option(
    "--preset",
    type=click.Choice(["dna-lamp", "rt-lamp", "fast-lamp"], case_sensitive=False),
    default=None,
    help=("TTP chemistry preset. Defaults to 'rt-lamp' when --na-type rna, otherwise 'dna-lamp'."),
)
@click.option(
    "--device-window",
    "device_window",
    type=float,
    default=60.0,
    show_default=True,
    help="Device run window (minutes). Default 60 min (BioVind BioID).",
)
@click.option(
    "--rt-min-tm",
    "rt_min_tm",
    type=float,
    default=63.0,
    show_default=True,
    help=(
        "Minimum Tm (degC) for core primers in one-step RT-LAMP. Only applies when --na-type rna."
    ),
)
@click.option(
    "--top-n",
    "top_n",
    default=None,
    type=int,
    help="Assess only the top-N ranked sets (default: all).",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write result summary to a key-value CSV file.",
)
def preorder(
    input_path: Path,
    na_type: str,
    sample_volume_ul: float,
    extraction_efficiency: float,
    eluate_volume_ul: float,
    reaction_input_ul: float,
    preset: str | None,
    device_window: float,
    rt_min_tm: float,
    top_n: int | None,
    out_csv: Path | None,
) -> None:
    r"""Combined pre-order readiness check: LOD + TTP + RT-LAMP in one pass.

    Chains three checks into a single GO / WARNING / NO-GO verdict:

    \b
      1. LOD_95 -- Poisson model: how many copies/mL are detectable at 95%?
      2. TTP at LOD_95 -- will the reaction turn positive in time on the device?
      3. RT-LAMP Tm check (RNA targets only) -- are core primers >= 63 degC?

    Use this before ordering primers to confirm the designed assay will work
    end-to-end on your specific device / extraction chain.

    \b
    Example (PRRSV ORF7, 1 mL oral fluid, 50% RNA extraction, 60-min device):
        lamp-forge preorder \\
          --input results/prrsv_orf7/primer_sets.json \\
          --na-type rna \\
          --sample-volume 1000 \\
          --efficiency 0.50 \\
          --eluate-volume 100 \\
          --reaction-input 5 \\
          --device-window 60

    \b
    Example (SRB dsrB, 1 mL produced water, 40% DNA extraction):
        lamp-forge preorder \\
          --input results/srb_dsrB/primer_sets.json \\
          --na-type dna \\
          --sample-volume 1000 \\
          --efficiency 0.40 \\
          --eluate-volume 100 \\
          --reaction-input 5
    """
    import json

    from lamp_forge.lod import ExtractionParams
    from lamp_forge.preorder import PreorderStatus, run_preorder, write_preorder_csv
    from lamp_forge.rt_lamp import RtLampParams, TargetNucleicAcid
    from lamp_forge.ttp import TtpParams, TtpPreset

    with input_path.open(encoding="utf-8") as fh:
        data: dict[str, object] = json.load(fh)

    sets_data = data.get("primer_sets", [])
    if not isinstance(sets_data, list) or not sets_data:
        click.secho("No primer sets found in the input file.", fg="red", err=True)
        sys.exit(1)

    try:
        extraction = ExtractionParams(
            sample_volume_ul=sample_volume_ul,
            extraction_efficiency=extraction_efficiency,
            eluate_volume_ul=eluate_volume_ul,
            reaction_input_ul=reaction_input_ul,
        )
    except ValueError as exc:
        click.secho(f"Extraction parameter error: {exc}", fg="red", err=True)
        sys.exit(2)

    target_na = TargetNucleicAcid(na_type.lower())

    resolved_preset = preset or ("rt-lamp" if target_na is TargetNucleicAcid.RNA else "dna-lamp")
    preset_enum = TtpPreset(resolved_preset.lower())
    ttp_params = TtpParams.from_preset(preset_enum, device_window_min=device_window)

    try:
        rt_params = RtLampParams(rt_min_tm=rt_min_tm)
    except ValueError as exc:
        click.secho(f"RT-LAMP parameter error: {exc}", fg="red", err=True)
        sys.exit(2)

    result = run_preorder(
        sets_data,
        extraction,
        ttp_params,
        target_na=target_na,
        rt_params=rt_params,
        top_n=top_n,
    )

    # --- header ---------------------------------------------------------------
    click.echo(
        f"Pre-order readiness check -- {na_type.upper()} target, "
        f"{result.n_sets_assessed} set(s) assessed"
    )
    click.echo(
        f"Extraction: {sample_volume_ul:.0f} uL sample, "
        f"{extraction_efficiency * 100:.0f}% efficiency, "
        f"{eluate_volume_ul:.0f} uL eluate, "
        f"{reaction_input_ul:.1f} uL to reaction"
    )
    click.echo(
        f"TTP model: {resolved_preset} "
        f"(TTP@1cp={ttp_params.ttp_one_copy_min:.1f} min, "
        f"slope={ttp_params.slope_min_per_decade:.1f} min/decade)"
    )
    click.echo(f"Device window: {device_window:.0f} min")
    click.echo("")

    # --- LOD ------------------------------------------------------------------
    eff_ul = extraction.copies_per_rxn_per_copy_per_ul
    click.echo(
        f"LOD_95 :  {result.lod_95_copies_per_rxn:.3f} copies/reaction  =>  "
        f"{result.lod_95_copies_per_ml:.1f} copies/mL  "
        f"(effective sample {eff_ul:.2f} uL/reaction)"
    )

    # --- TTP ------------------------------------------------------------------
    ttp_status = "PASS" if result.ttp_in_window else "FAIL"
    click.echo(f"TTP at LOD_95 :  {result.ttp_at_lod_min:.1f} min  [{ttp_status}]")

    # --- RT-LAMP --------------------------------------------------------------
    if target_na is TargetNucleicAcid.RNA:
        rt_ok = result.n_sets_rt_ok
        rt_total = result.n_sets_assessed
        rt_label = "OK" if rt_ok == rt_total else "WARNING"
        click.echo(f"RT-LAMP check :  {rt_ok}/{rt_total} sets optimized  [{rt_label}]")
    else:
        click.echo("RT-LAMP check :  N/A (DNA target)")

    click.echo("")

    # --- Verdict --------------------------------------------------------------
    if result.status is PreorderStatus.GO:
        click.secho(f"Overall status: {result.status}", fg="green")
    elif result.status is PreorderStatus.WARNING:
        click.secho(f"Overall status: {result.status}", fg="yellow")
    else:
        click.secho(f"Overall status: {result.status}", fg="red")

    for reason in result.reasons:
        click.echo(f"  {reason}")

    if out_csv is not None:
        write_preorder_csv(result, out_csv)
        click.echo(f"\nResult written to {out_csv}")


@cli.command(name="souring-trend")
@click.option(
    "--mic-result",
    "mic_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge mic-risk --out-json' for one monitoring "
        "interval.  Repeat once per sample in chronological order (oldest first). "
        "Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "srb, mcr, irb, apb, nrb.  Rows must be in chronological order "
        "(oldest first).  Mutually exclusive with --mic-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological guild-flag timeline to a CSV file.",
)
def souring_trend(
    mic_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse souring trajectory across multiple consecutive LAMP panel results.

    Takes a chronological sequence of oilfield MIC LAMP panel results and
    produces a trend assessment: IMPROVING, STABLE_LOW, STABLE_HIGH,
    DETERIORATING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - SRB newly detected (first positive after consecutive negatives)
      - SRB clearing (positive in early samples, negative in most recent)
      - Treatment under-dose (SRB + NRB co-detected in multiple intervals)

    Input (pick one):

    \b
      --mic-result FILE   JSON from 'lamp-forge mic-risk --out-json', repeated
                          once per monitoring interval, oldest first.
      --csv FILE          Monitoring-spreadsheet CSV (header required):
                            sample_id, date (opt), srb, mcr, irb, apb, nrb

    \b
    Example -- three-month trend from mic-risk JSON outputs:
        lamp-forge souring-trend \\
          --mic-result results/W1/2026-01.json \\
          --mic-result results/W1/2026-02.json \\
          --mic-result results/W1/2026-03.json \\
          --out-json results/W1/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge souring-trend \\
          --csv monitoring/W1_2026.csv \\
          --out-json results/W1/trend_2026.json
    """
    from lamp_forge.souring_trend import (
        analyse_souring_trend,
        records_from_csv,
        records_from_mic_json,
        write_trend_csv,
        write_trend_json,
    )

    # Validate mutual exclusion
    if csv_path is not None and mic_result_paths:
        click.secho(
            "--csv and --mic-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not mic_result_paths:
        click.secho(
            "Provide input via --mic-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_mic_json(list(mic_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_souring_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Display
    # -----------------------------------------------------------------------
    direction_colour = {
        "IMPROVING": "green",
        "STABLE_LOW": "green",
        "STABLE_HIGH": "yellow",
        "DETERIORATING": "red",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Souring trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    # Timeline table
    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = (
        f"{'Sample ID':<{col_id}}  {'Date':<12}  "
        f"{'SRB':>4} {'MCR':>4} {'IRB':>4} {'APB':>4} {'NRB':>4}  Risk"
    )
    click.echo(header)
    click.echo("-" * len(header))
    from lamp_forge.mic_risk import assess_mic_risk

    for r in trend.records:
        lvl = assess_mic_risk(r.flags).risk_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.srb else '[-]':>4} "
            f"{'[+]' if r.flags.mcr else '[-]':>4} "
            f"{'[+]' if r.flags.irb else '[-]':>4} "
            f"{'[+]' if r.flags.apb else '[-]':>4} "
            f"{'[+]' if r.flags.nrb else '[-]':>4}  {lvl}"
        )
    click.echo("")

    # Summary
    click.secho(f"Trend direction:  {trend.direction.value}", fg=col)
    click.echo(f"SRB positives:    {trend.srb_detection_count} of {trend.n_samples} interval(s)")
    if trend.nrb_detection_count > 0:
        click.echo(
            f"NRB positives:    {trend.nrb_detection_count} of {trend.n_samples} interval(s)"
        )
    if trend.treatment_underdosed_count > 0:
        click.secho(
            f"Under-dose flag:  SRB+NRB co-detected in {trend.treatment_underdosed_count} "
            "interval(s) -- nitrate treatment insufficient.",
            fg="yellow",
        )
    if trend.srb_newly_detected:
        click.secho("Event:            SRB newly detected (onset event).", fg="red")
    if trend.srb_clearing:
        click.secho("Event:            SRB clearing (treatment effective).", fg="green")
    click.echo(f"Worst risk level: {trend.worst_risk_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="farm-trend")
@click.option(
    "--farm-result",
    "farm_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge farm-risk --out-json' for one monitoring "
        "interval. Repeat once per sample in chronological order (oldest first). "
        "Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "asfv, fmdv, aiv, ndv, prrsv. Rows must be in chronological order "
        "(oldest first). Mutually exclusive with --farm-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def farm_trend(
    farm_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse biosecurity trajectory across consecutive farm LAMP panel results.

    Takes a chronological sequence of farm-biosecurity LAMP panel results and
    produces a trend assessment: EMERGING, STABLE_CLEAR, STABLE_ENDEMIC,
    RESOLVING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - WOAH-notifiable pathogen (ASFV/FMDV/AIV/NDV) newly detected
      - WOAH-notifiable pathogen cleared
      - PRRSV endemic herd status (consistently positive)

    Input (pick one):

    \b
      --farm-result FILE   JSON from 'lamp-forge farm-risk --out-json', repeated
                           once per monitoring interval, oldest first.
      --csv FILE           Monitoring-spreadsheet CSV (header required):
                             sample_id, date (opt), asfv, fmdv, aiv, ndv, prrsv

    \b
    Example -- three-month outbreak track on a poultry farm (AIV detected):
        lamp-forge farm-trend \\
          --farm-result results/farm-A/2026-01.json \\
          --farm-result results/farm-A/2026-02.json \\
          --farm-result results/farm-A/2026-03.json \\
          --out-json results/farm-A/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge farm-trend \\
          --csv monitoring/farm_A_2026.csv \\
          --out-json results/farm-A/trend_2026.json
    """
    from lamp_forge.farm_risk import assess_farm_risk as _assess_farm_risk
    from lamp_forge.farm_trend import (
        analyse_farm_trend,
        records_from_csv,
        records_from_farm_json,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and farm_result_paths:
        click.secho(
            "--csv and --farm-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not farm_result_paths:
        click.secho(
            "Provide input via --farm-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_farm_json(list(farm_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_farm_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Farm biosecurity trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = (
        f"{'Sample ID':<{col_id}}  {'Date':<12}  "
        f"{'ASFV':>5} {'FMDV':>5} {'AIV':>5} {'NDV':>5} {'PRRSV':>5}  Alert"
    )
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess_farm_risk(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.asfv else '[-]':>5} "
            f"{'[+]' if r.flags.fmdv else '[-]':>5} "
            f"{'[+]' if r.flags.aiv else '[-]':>5} "
            f"{'[+]' if r.flags.ndv else '[-]':>5} "
            f"{'[+]' if r.flags.prrsv else '[-]':>5}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.notifiable_newly_detected:
        click.secho(
            f"  NEWLY DETECTED (notifiable): {', '.join(trend.notifiable_newly_detected)}",
            fg="red",
            bold=True,
        )
    if trend.notifiable_cleared:
        click.secho(
            f"  Cleared (notifiable): {', '.join(trend.notifiable_cleared)}",
            fg="green",
        )
    if trend.prrsv_endemic:
        click.secho("  PRRSV endemic: detected in every interval.", fg="yellow")
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="poc-trend")
@click.option(
    "--poc-result",
    "poc_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge poc-risk --out-json' for one monitoring "
        "interval. Repeat once per sample in chronological order (oldest first). "
        "Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "mtb, cdiff, sars2, flu_a, gas. Rows must be in chronological order "
        "(oldest first). Mutually exclusive with --poc-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def poc_trend(
    poc_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse outbreak trajectory across consecutive human POC LAMP panel results.

    Takes a chronological sequence of human point-of-care LAMP panel results
    and produces a clinic surveillance trend assessment: EMERGING,
    STABLE_CLEAR, STABLE_ENDEMIC, RESOLVING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - MTB (TB) newly detected -- mandatory public-health notification
      - MTB clearing -- treatment response signal
      - Respiratory wave active (SARS-CoV-2 and/or Influenza A in recent samples)

    Input (pick one):

    \b
      --poc-result FILE   JSON from 'lamp-forge poc-risk --out-json', repeated
                          once per monitoring interval, oldest first.
      --csv FILE          Monitoring-spreadsheet CSV (header required):
                            sample_id, date (opt), mtb, cdiff, sars2, flu_a, gas

    \b
    Example -- three-month respiratory wave track at a rural clinic:
        lamp-forge poc-trend \\
          --poc-result results/clinic-A/2026-01.json \\
          --poc-result results/clinic-A/2026-02.json \\
          --poc-result results/clinic-A/2026-03.json \\
          --out-json results/clinic-A/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge poc-trend \\
          --csv monitoring/clinic_A_2026.csv \\
          --out-json results/clinic-A/trend_2026.json
    """
    from lamp_forge.poc_risk import assess_poc_risk as _assess_poc_risk
    from lamp_forge.poc_trend import (
        analyse_poc_trend,
        records_from_csv,
        records_from_poc_json,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and poc_result_paths:
        click.secho(
            "--csv and --poc-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not poc_result_paths:
        click.secho(
            "Provide input via --poc-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_poc_json(list(poc_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_poc_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"POC clinic surveillance trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = (
        f"{'Sample ID':<{col_id}}  {'Date':<12}  "
        f"{'MTB':>5} {'CDIFF':>5} {'SARS2':>5} {'FLU_A':>5} {'GAS':>5}  Alert"
    )
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess_poc_risk(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.mtb else '[-]':>5} "
            f"{'[+]' if r.flags.cdiff else '[-]':>5} "
            f"{'[+]' if r.flags.sars2 else '[-]':>5} "
            f"{'[+]' if r.flags.flu_a else '[-]':>5} "
            f"{'[+]' if r.flags.gas else '[-]':>5}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.tb_newly_detected:
        click.secho(
            "  TB NEWLY DETECTED: mandatory public-health notification required.",
            fg="red",
            bold=True,
        )
    if trend.tb_clearing:
        click.secho(
            "  TB clearing: MTB negative in most-recent sample (was positive earlier).",
            fg="green",
        )
    if trend.respiratory_wave_active:
        click.secho(
            "  Respiratory wave active: SARS-CoV-2 / Influenza A dominant in recent samples.",
            fg="yellow",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="bov-trend")
@click.option(
    "--bov-result",
    "bov_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge bov-risk --out-json' for one monitoring "
        "interval. Repeat once per sample in chronological order (oldest first). "
        "Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "brsv, bcov, bvdv, ibr, mhae. Rows must be in chronological order "
        "(oldest first). Mutually exclusive with --bov-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def bov_trend(
    bov_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse BRD trajectory across consecutive bovine respiratory LAMP panel results.

    Takes a chronological sequence of bovine respiratory LAMP panel results and
    produces a trend assessment: EMERGING, STABLE_CLEAR, STABLE_ENDEMIC,
    RESOLVING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - IBR (BoHV-1) newly detected -- mandatory notification may be required
        in EU/UK/Scandinavian eradication programmes
      - IBR cleared -- relevant for mandatory control programme compliance
      - BVDV endemic (positive every interval) -- persistently infected (PI)
        animal shedding suspected; individual antigen-ELISA screening indicated
      - Bacterial co-infection count -- M. haemolytica + viral co-detection
        intervals (highest-mortality BRD pattern)

    Input (pick one):

    \b
      --bov-result FILE   JSON from 'lamp-forge bov-risk --out-json', repeated
                          once per monitoring interval, oldest first.
      --csv FILE          Monitoring-spreadsheet CSV (header required):
                            sample_id, date (opt), brsv, bcov, bvdv, ibr, mhae

    \b
    Example -- three-month BRD track on a beef feedlot:
        lamp-forge bov-trend \\
          --bov-result results/herd-A/2026-01.json \\
          --bov-result results/herd-A/2026-02.json \\
          --bov-result results/herd-A/2026-03.json \\
          --out-json results/herd-A/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge bov-trend \\
          --csv monitoring/herd_A_2026.csv \\
          --out-json results/herd-A/trend_2026.json
    """
    from lamp_forge.bov_risk import assess_bov_risk as _assess_bov_risk
    from lamp_forge.bov_trend import (
        analyse_bov_trend,
        records_from_bov_json,
        records_from_csv,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and bov_result_paths:
        click.secho(
            "--csv and --bov-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not bov_result_paths:
        click.secho(
            "Provide input via --bov-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_bov_json(list(bov_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_bov_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Bovine respiratory BRD trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = (
        f"{'Sample ID':<{col_id}}  {'Date':<12}  "
        f"{'BRSV':>5} {'BCOV':>5} {'BVDV':>5} {'IBR':>5} {'MHAE':>5} {'MBOVIS':>6}  Alert"
    )
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess_bov_risk(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.brsv else '[-]':>5} "
            f"{'[+]' if r.flags.bcov else '[-]':>5} "
            f"{'[+]' if r.flags.bvdv else '[-]':>5} "
            f"{'[+]' if r.flags.ibr else '[-]':>5} "
            f"{'[+]' if r.flags.mhae else '[-]':>5} "
            f"{'[+]' if r.flags.mbovis else '[-]':>6}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.ibr_newly_detected:
        click.secho(
            "  IBR (BoHV-1) NEWLY DETECTED: check mandatory notification obligations "
            "(EU, UK, Scandinavia).",
            fg="red",
            bold=True,
        )
    if trend.ibr_cleared:
        click.secho(
            "  IBR cleared: confirm with serology before lifting eradication programme "
            "obligations.",
            fg="green",
        )
    if trend.bvdv_endemic:
        click.secho(
            "  BVDV endemic: detected in every interval -- PI animal screening (antigen "
            "ELISA) is indicated.",
            fg="yellow",
        )
    if trend.bacterial_coinfection_count > 0:
        click.secho(
            f"  Bacterial co-infection: M. haemolytica + viral co-detection in "
            f"{trend.bacterial_coinfection_count} interval(s) "
            f"(highest-mortality BRD pattern).",
            fg="yellow",
        )
    if trend.mbovis_detected_count > 0:
        click.secho(
            f"  M. bovis detected in {trend.mbovis_detected_count} interval(s): "
            "beta-lactam/aminoglycoside protocols are INEFFECTIVE; submit for "
            "susceptibility testing.",
            fg="yellow",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="brucella-risk")
@click.option(
    "--is711/--no-is711",
    "is711_positive",
    default=False,
    help="IS711 LAMP assay result (positive / negative).",
)
@click.option(
    "--context",
    "context",
    type=click.Choice(
        ["cattle", "small_ruminant", "swine", "human", "environmental"],
        case_sensitive=False,
    ),
    default="cattle",
    show_default=True,
    help=(
        "Sample context for species-level interpretation. "
        "cattle=B.abortus, small_ruminant=B.melitensis, swine=B.suis, "
        "human=human brucellosis POC, environmental=fomite/water/soil."
    ),
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read flags from a JSON file instead of CLI options. "
        "Expected keys: is711_positive (bool), context (str). "
        "Missing keys use defaults."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def brucella_risk(
    is711_positive: bool,
    context: str,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a Brucella IS711 LAMP result as a structured alert.

    IS711 is a pan-Brucella insertion sequence (5-40 copies per genome) that
    enables genus-level detection of all nine Brucella species with exceptional
    analytical sensitivity.  The alert level depends on the sample context:
    livestock samples trigger WOAH-notifiable HIGH alerts; a human sample
    triggers a CRITICAL public-health alert; an environmental positive triggers
    a LOW alert requiring follow-up animal testing.

    \b
    Alert levels:
        CRITICAL  -- IS711+ in a human sample (notifiable zoonosis; treat immediately)
        HIGH      -- IS711+ in cattle, small ruminant, or swine (WOAH-listed; quarantine)
        LOW       -- IS711+ in an environmental sample (follow-up animal testing required)
        NEGATIVE  -- IS711 negative (Brucella not detected)

    \b
    Example -- positive cattle blood sample:
        lamp-forge brucella-risk --is711 --context cattle

    \b
    Example -- positive human POC sample:
        lamp-forge brucella-risk --is711 --context human

    \b
    Example -- read flags from JSON and write outputs:
        lamp-forge brucella-risk --input-json results/brucella_flags.json \
          --out-json results/brucella_assessment.json \
          --out-csv results/brucella_assessment.csv
    """
    import json as json_mod

    from lamp_forge.brucella_risk import (
        BrucellaAlertLevel,
        BrucellaContext,
        BrucellaFlags,
        assess_brucella_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = BrucellaFlags(
            is711_positive=is711_positive,
            context=BrucellaContext(context.lower()),
        )

    assessment = assess_brucella_risk(flags)

    # --- Result table ---------------------------------------------------------
    click.echo("Brucella IS711 LAMP result:")
    symbol = "+" if assessment.flags.is711_positive else "-"
    click.echo(f"  IS711  (pan-Brucella)  [{symbol}]")
    click.echo(f"  Context: {assessment.flags.context.value}")
    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color: dict[BrucellaAlertLevel, str] = {
        BrucellaAlertLevel.CRITICAL: "red",
        BrucellaAlertLevel.HIGH: "red",
        BrucellaAlertLevel.LOW: "yellow",
        BrucellaAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level in (BrucellaAlertLevel.CRITICAL, BrucellaAlertLevel.HIGH)),
    )

    if assessment.immediate_report_required:
        click.secho(
            "  IMMEDIATE REPORT REQUIRED: notify the national veterinary or public "
            "health authority today.",
            fg="red",
            bold=True,
        )

    if assessment.notifiable:
        click.secho(
            f"  NOTIFIABLE: zoonotic_risk={assessment.zoonotic_risk}",
            fg="yellow",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="brucella-trend")
@click.option(
    "--brucella-result",
    "brucella_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge brucella-risk --out-json' for one "
        "monitoring interval. Repeat once per sample in chronological order (oldest "
        "first). Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, is711_positive, "
        "date (optional), context (optional; default cattle). "
        "Rows must be in chronological order (oldest first). "
        "Mutually exclusive with --brucella-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological IS711 timeline to a CSV file.",
)
def brucella_trend(
    brucella_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse Brucella IS711 surveillance trajectory across consecutive LAMP results.

    Takes a chronological sequence of Brucella IS711 LAMP results and produces a
    trend assessment: NEWLY_DETECTED, PERSISTENT, RESOLVING, STABLE_CLEAR, or
    INSUFFICIENT_DATA.

    Key events detected:

    \b
      - IS711 newly detected -- WOAH-notifiable event; herd quarantine required
      - IS711 cleared -- confirm with serology before lifting quarantine
      - Persistent infection likely -- >= 75% of >= 3 samples positive

    Input (pick one):

    \b
      --brucella-result FILE   JSON from 'lamp-forge brucella-risk --out-json',
                               repeated once per monitoring interval, oldest first.
      --csv FILE               Monitoring-spreadsheet CSV (header required):
                                 sample_id, is711_positive, date (opt), context (opt)

    \b
    Example -- quarterly herd surveillance, cattle:
        lamp-forge brucella-trend \
          --brucella-result results/herd-A/2026-01.json \
          --brucella-result results/herd-A/2026-02.json \
          --brucella-result results/herd-A/2026-03.json \
          --out-json results/herd-A/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge brucella-trend \
          --csv monitoring/herd_A_2026.csv \
          --out-json results/herd-A/trend_2026.json
    """
    from lamp_forge.brucella_risk import assess_brucella_risk as _assess
    from lamp_forge.brucella_trend import (
        analyse_brucella_trend,
        records_from_brucella_json,
        records_from_csv,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and brucella_result_paths:
        click.secho(
            "--csv and --brucella-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not brucella_result_paths:
        click.secho(
            "Provide input via --brucella-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_brucella_json(list(brucella_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_brucella_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "NEWLY_DETECTED": "red",
        "PERSISTENT": "red",
        "RESOLVING": "green",
        "STABLE_CLEAR": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Brucella IS711 surveillance trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = f"{'Sample ID':<{col_id}}  {'Date':<12}  {'IS711':>5}  {'Context':<14}  Alert"
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.is711_positive else '[-]':>5}  "
            f"{r.flags.context.value:<14}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.newly_detected:
        click.secho(
            "  IS711 NEWLY DETECTED: notify the national veterinary or public health "
            "authority immediately and initiate herd quarantine.",
            fg="red",
            bold=True,
        )
    if trend.cleared:
        click.secho(
            "  IS711 cleared: confirm with serology before lifting quarantine.",
            fg="green",
        )
    if trend.persistent_infection_likely:
        click.secho(
            "  Persistent infection likely: conduct whole-herd serology to identify "
            "and remove the source animal.",
            fg="yellow",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="csfv-risk")
@click.option(
    "--ns5b/--no-ns5b",
    "ns5b_positive",
    default=False,
    help="CSFV NS5B RT-LAMP assay result (positive / negative).",
)
@click.option(
    "--context",
    "context",
    type=click.Choice(["on_farm", "border_post", "market"], case_sensitive=False),
    default="on_farm",
    show_default=True,
    help=(
        "Deployment context for the recommended action. "
        "on_farm=herd surveillance or post-morbidity investigation, "
        "border_post=border inspection post or port of entry, "
        "market=livestock market or auction site."
    ),
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read flags from a JSON file instead of CLI options. "
        "Expected keys: ns5b_positive (bool), context (str). "
        "Missing keys use defaults."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def csfv_risk(
    ns5b_positive: bool,
    context: str,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a CSFV NS5B RT-LAMP result as a structured alert.

    Classical Swine Fever Virus (CSFV; hog cholera) is a WOAH-listed Tier 1
    notifiable disease causing near-100% case fatality in naive pig herds.
    It is clinically indistinguishable from African Swine Fever Virus (ASFV) --
    a BioVind portable device running CSFV and ASFV assays in the same
    cartridge resolves this differential at the farm gate in under 60 minutes.

    NS5B (RNA-directed RNA polymerase) is the most conserved CSFV gene across
    all three genotype groups (1/2/3) and requires one-step RT-LAMP.

    Alert levels:
    \b
        CRITICAL  -- NS5B positive (any context): WOAH-notifiable, near-100% CFR;
                     immediate herd quarantine and national-authority notification.
        NEGATIVE  -- NS5B negative: CSFV not detected.

    \b
    Example -- positive on-farm detection (differential with ASFV):
        lamp-forge csfv-risk --ns5b --context on_farm

    \b
    Example -- positive at a border post:
        lamp-forge csfv-risk --ns5b --context border_post

    \b
    Example -- read flags from JSON and write outputs:
        lamp-forge csfv-risk --input-json results/csfv_flags.json \
          --out-json results/csfv_assessment.json \
          --out-csv results/csfv_assessment.csv
    """
    import json as json_mod

    from lamp_forge.csfv_risk import (
        CsfvAlertLevel,
        CsfvContext,
        CsfvFlags,
        assess_csfv_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = CsfvFlags(
            ns5b_positive=ns5b_positive,
            context=CsfvContext(context.lower()),
        )

    assessment = assess_csfv_risk(flags)

    # --- Result table ---------------------------------------------------------
    click.echo("CSFV NS5B RT-LAMP result:")
    symbol = "+" if assessment.flags.ns5b_positive else "-"
    click.echo(f"  NS5B  (pan-CSFV RNA polymerase)  [{symbol}]")
    click.echo(f"  Context: {assessment.flags.context.value}")
    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color: dict[CsfvAlertLevel, str] = {
        CsfvAlertLevel.CRITICAL: "red",
        CsfvAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level is CsfvAlertLevel.CRITICAL),
    )

    if assessment.immediate_report_required:
        click.secho(
            "  IMMEDIATE REPORT REQUIRED: notify the national veterinary authority today.",
            fg="red",
            bold=True,
        )

    if assessment.woah_notifiable:
        click.secho(
            "  WOAH-NOTIFIABLE: CSFV is a Tier 1 listed disease -- official "
            "confirmation at a reference laboratory is required before depopulation.",
            fg="yellow",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="csfv-trend")
@click.option(
    "--csfv-result",
    "csfv_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge csfv-risk --out-json' for one "
        "monitoring interval. Repeat once per sample in chronological order (oldest "
        "first). Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, ns5b_positive, "
        "date (optional), context (optional; default on_farm). "
        "Rows must be in chronological order (oldest first). "
        "Mutually exclusive with --csfv-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological NS5B timeline to a CSV file.",
)
def csfv_trend(
    csfv_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse CSFV NS5B RT-LAMP surveillance trajectory across consecutive results.

    Takes a chronological sequence of CSFV NS5B RT-LAMP results and produces a
    trend assessment: NEWLY_DETECTED, ONGOING_OUTBREAK, CONTROLLED, STABLE_CLEAR,
    or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - NS5B newly detected   -- WOAH-notifiable; emergency herd quarantine required
      - NS5B controlled       -- apparent containment; confirm before lifting quarantine
      - Sustained outbreak    -- >= 75% of >= 3 samples positive; active circulation

    Input (pick one):

    \b
      --csfv-result FILE   JSON from 'lamp-forge csfv-risk --out-json',
                           repeated once per monitoring interval, oldest first.
      --csv FILE           Monitoring-spreadsheet CSV (header required):
                             sample_id, ns5b_positive, date (opt), context (opt)

    \b
    Example -- monthly farm surveillance:
        lamp-forge csfv-trend \
          --csfv-result results/farm-A/2026-01.json \
          --csfv-result results/farm-A/2026-02.json \
          --csfv-result results/farm-A/2026-03.json \
          --out-json results/farm-A/csfv_trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge csfv-trend \
          --csv monitoring/farm_A_2026.csv \
          --out-json results/farm-A/csfv_trend_2026.json
    """
    from lamp_forge.csfv_risk import assess_csfv_risk as _assess
    from lamp_forge.csfv_trend import (
        analyse_csfv_trend,
        records_from_csfv_json,
        records_from_csv,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and csfv_result_paths:
        click.secho(
            "--csv and --csfv-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not csfv_result_paths:
        click.secho(
            "Provide input via --csfv-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_csfv_json(list(csfv_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_csfv_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "NEWLY_DETECTED": "red",
        "ONGOING_OUTBREAK": "red",
        "CONTROLLED": "green",
        "STABLE_CLEAR": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"CSFV NS5B RT-LAMP surveillance trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = f"{'Sample ID':<{col_id}}  {'Date':<12}  {'NS5B':>5}  {'Context':<12}  Alert"
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.ns5b_positive else '[-]':>5}  "
            f"{r.flags.context.value:<12}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.newly_detected:
        click.secho(
            "  NS5B NEWLY DETECTED: initiate emergency herd quarantine and notify "
            "the national veterinary authority immediately (WOAH-listed disease).",
            fg="red",
            bold=True,
        )
    if trend.controlled:
        click.secho(
            "  NS5B controlled: apparent containment -- maintain quarantine and "
            "confirm with at least three consecutive negatives before review.",
            fg="green",
        )
    if trend.sustained_outbreak:
        click.secho(
            "  Sustained outbreak: >= 75% of samples positive -- active CSFV "
            "circulation; review containment and escalate to national authority.",
            fg="yellow",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="swine-enteric-risk")
@click.option(
    "--pedv",
    "pedv",
    is_flag=True,
    default=False,
    help="PEDV (N gene) positive.",
)
@click.option(
    "--pdcov",
    "pdcov",
    is_flag=True,
    default=False,
    help="Porcine deltacoronavirus (N gene) positive.",
)
@click.option(
    "--rota-a",
    "rota_a",
    is_flag=True,
    default=False,
    help="Porcine rotavirus A (VP6) positive.",
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: pedv, pdcov, rota_a (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def swine_enteric_risk(
    pedv: bool,
    pdcov: bool,
    rota_a: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a swine neonatal enteric LAMP panel as a structured alert.

    Takes the positivity flags for the three pathogens monitored by a
    BioVind-style swine neonatal enteric panel and returns a structured alert
    with level, score, interpretation, and recommended action.

    PEDV (Porcine Epidemic Diarrhea Virus) causes near-100% mortality in
    neonatal piglets and triggers immediate quarantine and barn lockdown.

    Panel targets:

    \b
        --pedv     Porcine epidemic diarrhea virus    N gene  (~100% CFR neonates)
        --pdcov    Porcine deltacoronavirus           N gene  (40-80% CFR neonates)
        --rota-a   Porcine rotavirus A                VP6     (10-40% CFR neonates)

    \b
    Example -- PEDV detected in a farrowing barn:
        lamp-forge swine-enteric-risk --pedv

    \b
    Example -- PDCoV and Rotavirus A co-detected:
        lamp-forge swine-enteric-risk --pdcov --rota-a

    \b
    Example -- full panel from a JSON flags file:
        lamp-forge swine-enteric-risk --input-json results/enteric_flags.json \
          --out-json results/enteric_assessment.json
    """
    import json as json_mod

    from lamp_forge.swine_enteric_risk import (
        SwineEntericAlertLevel,
        SwineEntericFlags,
        assess_swine_enteric_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = SwineEntericFlags(pedv=pedv, pdcov=pdcov, rota_a=rota_a)

    assessment = assess_swine_enteric_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Swine neonatal enteric LAMP panel results:")
    target_rows = [
        ("PEDV", "N gene", assessment.flags.pedv),
        ("PDCoV", "N gene", assessment.flags.pdcov),
        ("RotaA", "VP6", assessment.flags.rota_a),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<6} ({gene:<6})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        SwineEntericAlertLevel.CRITICAL: "red",
        SwineEntericAlertLevel.HIGH: "red",
        SwineEntericAlertLevel.MODERATE: "yellow",
        SwineEntericAlertLevel.LOW: "yellow",
        SwineEntericAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(
            assessment.alert_level in (SwineEntericAlertLevel.CRITICAL, SwineEntericAlertLevel.HIGH)
        ),
    )

    if assessment.immediate_quarantine_required:
        click.secho(
            "  IMMEDIATE QUARANTINE REQUIRED: isolate affected pens; stop all pig movement.",
            fg="red",
            bold=True,
        )

    if assessment.biosecurity_lockdown:
        click.secho(
            "  BIOSECURITY LOCKDOWN: enhanced PPE, boot dips, one-way traffic; "
            "decontaminate before leaving the farrowing unit.",
            fg="red",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="swine-enteric-trend")
@click.option(
    "--enteric-result",
    "enteric_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge swine-enteric-risk --out-json' for one "
        "monitoring interval. Repeat once per sample in chronological order (oldest "
        "first). Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "pedv, pdcov, rota_a. Rows must be in chronological order "
        "(oldest first). Mutually exclusive with --enteric-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def swine_enteric_trend(
    enteric_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse enteric outbreak trajectory across consecutive swine neonatal LAMP panel results.

    Takes a chronological sequence of swine neonatal enteric LAMP panel results and
    produces a trend assessment: EMERGING, STABLE_CLEAR, STABLE_ENDEMIC,
    RESOLVING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - PEDV newly detected -- near-100% neonatal mortality; barn lockdown required
      - PEDV cleared -- confirm with RT-PCR before relaxing controls
      - Overall pathogen burden trend across the monitoring period

    Input (pick one):

    \b
      --enteric-result FILE   JSON from 'lamp-forge swine-enteric-risk --out-json',
                              repeated once per monitoring interval, oldest first.
      --csv FILE              Monitoring-spreadsheet CSV (header required):
                                sample_id, date (opt), pedv, pdcov, rota_a

    \b
    Example -- three-month PEDV outbreak track in a farrowing barn:
        lamp-forge swine-enteric-trend \\
          --enteric-result results/barn-A/2026-01.json \\
          --enteric-result results/barn-A/2026-02.json \\
          --enteric-result results/barn-A/2026-03.json \\
          --out-json results/barn-A/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge swine-enteric-trend \\
          --csv monitoring/barn_A_2026.csv \\
          --out-json results/barn-A/trend_2026.json
    """
    from lamp_forge.swine_enteric_risk import assess_swine_enteric_risk as _assess
    from lamp_forge.swine_enteric_trend import (
        analyse_swine_enteric_trend,
        records_from_csv,
        records_from_enteric_json,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and enteric_result_paths:
        click.secho(
            "--csv and --enteric-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not enteric_result_paths:
        click.secho(
            "Provide input via --enteric-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_enteric_json(list(enteric_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_swine_enteric_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Swine neonatal enteric trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = f"{'Sample ID':<{col_id}}  {'Date':<12}  {'PEDV':>5} {'PDCoV':>5} {'RotaA':>5}  Alert"
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.pedv else '[-]':>5} "
            f"{'[+]' if r.flags.pdcov else '[-]':>5} "
            f"{'[+]' if r.flags.rota_a else '[-]':>5}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.pedv_newly_detected:
        click.secho(
            "  PEDV NEWLY DETECTED: near-100% neonatal mortality -- initiate barn "
            "lockdown and quarantine all affected pens immediately.",
            fg="red",
            bold=True,
        )
    if trend.pedv_cleared:
        click.secho(
            "  PEDV cleared: confirm with RT-PCR before relaxing biosecurity controls.",
            fg="green",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="poultry-risk")
@click.option(
    "--aiv",
    "aiv",
    is_flag=True,
    default=False,
    help="Avian influenza A (M gene, pan-IAV) positive.",
)
@click.option(
    "--ndv",
    "ndv",
    is_flag=True,
    default=False,
    help="Newcastle disease virus (M gene, pan-NDV) positive.",
)
@click.option(
    "--ibdv",
    "ibdv",
    is_flag=True,
    default=False,
    help="Infectious bursal disease virus / Gumboro (VP2 gene) positive.",
)
@click.option(
    "--ibv",
    "ibv",
    is_flag=True,
    default=False,
    help="Infectious bronchitis virus (N gene, pan-IBV) positive.",
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: aiv, ndv, ibdv, ibv (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def poultry_risk(
    aiv: bool,
    ndv: bool,
    ibdv: bool,
    ibv: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a poultry biosecurity LAMP panel as a structured alert.

    Takes the positivity flags for four pathogens monitored by a BioVind-style
    portable poultry biosecurity panel and returns a structured alert with level,
    score, interpretation, and recommended action.

    AIV (Avian Influenza A) triggers a CRITICAL alert with immediate reporting
    and zoonotic-risk flags.  NDV triggers a HIGH alert.  IBDV and IBV trigger
    MODERATE and LOW alerts respectively.

    Panel targets:

    \b
        --aiv    Avian influenza A        M gene    WOAH Tier 1 notifiable; zoonotic
        --ndv    Newcastle disease virus  M gene    WOAH-listed notifiable
        --ibdv   Infect. bursal disease   VP2       Gumboro; immunosuppressive
        --ibv    Infect. bronchitis       N gene    most prevalent avian coronavirus

    \b
    Example -- AIV detected (presumptive HPAI):
        lamp-forge poultry-risk --aiv

    \b
    Example -- NDV and IBDV co-detected:
        lamp-forge poultry-risk --ndv --ibdv

    \b
    Example -- full panel from a JSON flags file:
        lamp-forge poultry-risk --input-json results/poultry_flags.json \
          --out-json results/poultry_assessment.json
    """
    import json as json_mod

    from lamp_forge.poultry_risk import (
        PoultryAlertLevel,
        PoultryFlags,
        assess_poultry_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = PoultryFlags(aiv=aiv, ndv=ndv, ibdv=ibdv, ibv=ibv)

    assessment = assess_poultry_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Poultry biosecurity LAMP panel results:")
    target_rows = [
        ("AIV", "M gene", assessment.flags.aiv),
        ("NDV", "M gene", assessment.flags.ndv),
        ("IBDV", "VP2  ", assessment.flags.ibdv),
        ("IBV", "N gene", assessment.flags.ibv),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<5} ({gene})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        PoultryAlertLevel.CRITICAL: "red",
        PoultryAlertLevel.HIGH: "red",
        PoultryAlertLevel.MODERATE: "yellow",
        PoultryAlertLevel.LOW: "yellow",
        PoultryAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level in (PoultryAlertLevel.CRITICAL, PoultryAlertLevel.HIGH)),
    )

    if assessment.immediate_report_required:
        click.secho(
            "  REPORT REQUIRED: notify the national veterinary authority "
            "before any birds are moved or culled.",
            fg="red",
            bold=True,
        )

    if assessment.zoonotic_risk:
        click.secho(
            "  ZOONOTIC RISK: all personnel in the affected house must wear "
            "N95 respirators, eye protection, and gloves.",
            fg="red",
            bold=True,
        )

    if assessment.depopulation_risk:
        click.secho(
            "  DEPOPULATION RISK: consult the national authority before "
            "any birds are moved; emergency depopulation may be required.",
            fg="red",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="poultry-trend")
@click.option(
    "--poultry-result",
    "poultry_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge poultry-risk --out-json' for one "
        "monitoring interval. Repeat once per sample in chronological order (oldest "
        "first). Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "aiv, ndv, ibdv, ibv. Rows must be in chronological order "
        "(oldest first). Mutually exclusive with --poultry-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def poultry_trend(
    poultry_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse biosecurity trajectory across consecutive poultry LAMP panel results.

    Takes a chronological sequence of poultry biosecurity LAMP panel results and
    produces a trend assessment: EMERGING, STABLE_CLEAR, STABLE_ENDEMIC,
    RESOLVING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - WOAH-notifiable pathogen (AIV/NDV) newly detected
      - WOAH-notifiable pathogen cleared
      - IBV endemic flock status (consistently positive)

    Input (pick one):

    \b
      --poultry-result FILE   JSON from 'lamp-forge poultry-risk --out-json',
                              repeated once per monitoring interval, oldest first.
      --csv FILE              Monitoring-spreadsheet CSV (header required):
                                sample_id, date (opt), aiv, ndv, ibdv, ibv

    \b
    Example -- quarterly AIV surveillance on a broiler farm:
        lamp-forge poultry-trend \\
          --poultry-result results/flock-A/2026-01.json \\
          --poultry-result results/flock-A/2026-02.json \\
          --poultry-result results/flock-A/2026-03.json \\
          --out-json results/flock-A/trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge poultry-trend \\
          --csv monitoring/flock_A_2026.csv \\
          --out-json results/flock-A/trend_2026.json
    """
    from lamp_forge.poultry_risk import assess_poultry_risk as _assess
    from lamp_forge.poultry_trend import (
        analyse_poultry_trend,
        records_from_csv,
        records_from_poultry_json,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and poultry_result_paths:
        click.secho(
            "--csv and --poultry-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not poultry_result_paths:
        click.secho(
            "Provide input via --poultry-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_poultry_json(list(poultry_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_poultry_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Poultry biosecurity trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = (
        f"{'Sample ID':<{col_id}}  {'Date':<12}  "
        f"{'AIV':>5} {'NDV':>5} {'IBDV':>5} {'IBV':>5}  Alert"
    )
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.aiv else '[-]':>5} "
            f"{'[+]' if r.flags.ndv else '[-]':>5} "
            f"{'[+]' if r.flags.ibdv else '[-]':>5} "
            f"{'[+]' if r.flags.ibv else '[-]':>5}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.notifiable_newly_detected:
        targets = ", ".join(trend.notifiable_newly_detected)
        click.secho(
            f"  {targets} NEWLY DETECTED: notify the national veterinary authority "
            "immediately; quarantine the affected house before any birds are moved.",
            fg="red",
            bold=True,
        )
    if trend.notifiable_cleared:
        targets_c = ", ".join(trend.notifiable_cleared)
        click.secho(
            f"  {targets_c} cleared: confirm with a WOAH Reference Laboratory "
            "before lifting movement restrictions.",
            fg="green",
        )
    if trend.ibv_endemic:
        click.echo(
            "  IBV endemic flock: consistently positive -- review IBV vaccination "
            "programme and serotype match with the flock health veterinarian."
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


@cli.command(name="cartridge")
@click.argument(
    "manifest_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--out-dir",
    "out_dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Output directory for cartridge_manifest.json and every per-target artifact.",
)
@click.option(
    "--no-html",
    "no_html",
    is_flag=True,
    default=False,
    help="Skip the HTML reports (panel_report.html and build_report.html).",
)
@click.option(
    "--vendor",
    "vendor_override",
    type=click.Choice(["idt", "twist"], case_sensitive=False),
    default=None,
    help="Override the manifest's vendor field for the vendor_order.csv export.",
)
def cartridge(
    manifest_path: Path,
    out_dir: Path,
    no_html: bool,
    vendor_override: str | None,
) -> None:
    """One-shot BioID cartridge build: panel + pool + LoD + RT + preorder + vendor.

    Loads a cartridge manifest YAML, composes every LAMP-Forge primitive
    needed for a multi-target SKU, and writes a signed cartridge_manifest.json
    plus the eight standard build artifacts.

    Exits 0 on GO/WARNING and non-zero on NO_GO -- wire this into a CI
    pre-order gate so dirty panels and failed RT checks block vendor order
    submission.

    \b
    Example -- 18-channel respiratory cartridge:
        lamp-forge cartridge config/cartridge_respiratory_18ch.yaml \\
          --out-dir results/resp_18ch_v1
    """
    from lamp_forge.cartridge import (
        load_cartridge_manifest,
        run_cartridge,
    )

    try:
        manifest = load_cartridge_manifest(manifest_path)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        click.secho(f"Manifest error: {exc}", fg="red", err=True)
        sys.exit(2)

    if vendor_override is not None:
        from dataclasses import replace

        manifest = replace(manifest, vendor=vendor_override.lower())  # type: ignore[arg-type]

    try:
        result = run_cartridge(manifest, out_dir, write_html=not no_html)
    except Exception as exc:
        click.secho(f"Cartridge build failed: {exc}", fg="red", err=True)
        if click.get_current_context().obj.get("verbose"):
            raise
        sys.exit(1)

    click.echo(
        f"Cartridge {manifest.cartridge_id}: "
        f"{len(manifest.targets)}/{manifest.max_channels} channel(s), "
        f"chemistry={manifest.chemistry}"
    )
    color = {"GO": "green", "WARNING": "yellow", "NO_GO": "red"}[result.build_status]
    click.secho(f"Build status: {result.build_status}", fg=color)
    for r in result.reasons:
        click.echo(f"  {r}")
    click.echo(f"Manifest hash: {result.manifest_hash}")
    click.echo(f"Outputs written to {out_dir}")

    if result.build_status == "NO_GO":
        sys.exit(1)


@cli.command(name="bov-enteric-risk")
@click.option(
    "--etec-k99",
    "etec_k99",
    is_flag=True,
    default=False,
    help="ETEC K99 / F5 fimbria (faeG gene) positive.",
)
@click.option(
    "--salmonella",
    "salmonella",
    is_flag=True,
    default=False,
    help="Salmonella sp. (invA gene) positive.",
)
@click.option(
    "--crypto",
    "crypto",
    is_flag=True,
    default=False,
    help="Cryptosporidium parvum (18S / hsp70 gene) positive.",
)
@click.option(
    "--bcov",
    "bcov",
    is_flag=True,
    default=False,
    help="Bovine coronavirus enteric pathotype (N gene) positive.",
)
@click.option(
    "--rota-a",
    "rota_a",
    is_flag=True,
    default=False,
    help="Bovine rotavirus A (VP6 gene) positive.",
)
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: etec_k99, salmonella, crypto, bcov, rota_a (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def bov_enteric_risk(
    etec_k99: bool,
    salmonella: bool,
    crypto: bool,
    bcov: bool,
    rota_a: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a bovine neonatal calf diarrhea LAMP panel as a structured alert.

    Takes the positivity flags for the five pathogens monitored by a
    BioVind-style bovine neonatal enteric panel and returns a structured alert
    with level, score, interpretation, and recommended action.

    The critical clinical question: does this calf need intravenous antibiotics
    in the next two hours (ETEC K99 / Salmonella) or supportive care only
    (BCoV / Rotavirus A / Cryptosporidium)?  The LAMP panel resolves this
    at the pen-side in under 60 minutes.

    Panel targets:

    \b
        --etec-k99    ETEC K99 / F5 fimbria      faeG   (~100% CFR neonates <4d)
        --salmonella  Salmonella sp.             invA   (40-70% CFR; zoonotic)
        --crypto      Cryptosporidium parvum     18S    (10-50% CFR; zoonotic)
        --bcov        Bovine coronavirus         N gene (5-20% CFR; supportive)
        --rota-a      Bovine rotavirus A         VP6    (5-15% CFR; supportive)

    \b
    Example -- ETEC K99 detected in a neonatal calf (<4 days old):
        lamp-forge bov-enteric-risk --etec-k99

    \b
    Example -- BCoV and Rotavirus A co-detected (supportive care):
        lamp-forge bov-enteric-risk --bcov --rota-a

    \b
    Example -- Salmonella detected (bacteremia risk; zoonotic):
        lamp-forge bov-enteric-risk --salmonella

    \b
    Example -- full panel from a JSON flags file:
        lamp-forge bov-enteric-risk --input-json results/enteric_flags.json \
          --out-json results/calf_assessment.json
    """
    import json as json_mod

    from lamp_forge.bov_enteric_risk import (
        BovEntericAlertLevel,
        BovEntericFlags,
        assess_bov_enteric_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = BovEntericFlags(
            etec_k99=etec_k99,
            salmonella=salmonella,
            crypto=crypto,
            bcov=bcov,
            rota_a=rota_a,
        )

    assessment = assess_bov_enteric_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Bovine neonatal calf diarrhea LAMP panel results:")
    target_rows = [
        ("ETEC_K99", "faeG", assessment.flags.etec_k99),
        ("Salmonella", "invA", assessment.flags.salmonella),
        ("Crypto", "18S", assessment.flags.crypto),
        ("BCoV", "N gene", assessment.flags.bcov),
        ("RotaA", "VP6", assessment.flags.rota_a),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<10} ({gene:<6})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        BovEntericAlertLevel.CRITICAL: "red",
        BovEntericAlertLevel.HIGH: "red",
        BovEntericAlertLevel.MODERATE: "yellow",
        BovEntericAlertLevel.LOW: "yellow",
        BovEntericAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level in (BovEntericAlertLevel.CRITICAL, BovEntericAlertLevel.HIGH)),
    )

    if assessment.antibiotic_indicated:
        click.secho(
            "  ANTIBIOTIC INDICATED: bacterial pathogen(s) detected. "
            "Initiate IV antibiotic therapy -- do NOT wait for culture results.",
            fg="red",
            bold=True,
        )

    if assessment.zoonotic_risk:
        click.secho(
            "  ZOONOTIC RISK: Salmonella and/or Cryptosporidium detected. "
            "Use gloves and strict handwashing when handling affected calves.",
            fg="yellow",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="bov-enteric-trend")
@click.option(
    "--enteric-result",
    "enteric_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge bov-enteric-risk --out-json' for one monitoring "
        "interval. Repeat once per sample in chronological order (oldest first). "
        "Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "etec_k99, salmonella, crypto, bcov, rota_a. Rows must be in chronological "
        "order (oldest first). Mutually exclusive with --enteric-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def bov_enteric_trend(
    enteric_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse neonatal calf diarrhea trajectory across consecutive enteric LAMP panel results.

    Takes a chronological sequence of bovine neonatal enteric LAMP panel results and
    produces a trend assessment: EMERGING, STABLE_CLEAR, STABLE_ENDEMIC,
    RESOLVING, or INSUFFICIENT_DATA.

    Key events detected:

    \b
      - ETEC K99 newly detected -- IV antibiotic therapy required within 2 h;
        near-certain mortality in calves <4 days old without treatment
      - ETEC K99 cleared -- current treatment is working; maintain hygiene
      - Salmonella persistent (positive in every interval) -- herd-level
        contamination or chronic shedder suspected; culture and investigation
      - Zoonotic risk count (Salmonella / Cryptosporidium intervals)

    Input (pick one):

    \b
      --enteric-result FILE   JSON from 'lamp-forge bov-enteric-risk --out-json',
                              repeated once per monitoring interval, oldest first.
      --csv FILE              Monitoring-spreadsheet CSV (header required):
                                sample_id, date (opt), etec_k99, salmonella,
                                crypto, bcov, rota_a

    \b
    Example -- three-week neonatal scours track on a dairy herd:
        lamp-forge bov-enteric-trend \
          --enteric-result results/calf-grp-1/week-01.json \
          --enteric-result results/calf-grp-1/week-02.json \
          --enteric-result results/calf-grp-1/week-03.json \
          --out-json results/calf-grp-1/trend_wk1-3.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge bov-enteric-trend \
          --csv monitoring/calf_grp1_2026.csv \
          --out-json results/calf-grp-1/trend_2026.json
    """
    from lamp_forge.bov_enteric_risk import assess_bov_enteric_risk as _assess_enteric
    from lamp_forge.bov_enteric_trend import (
        analyse_bov_enteric_trend,
        records_from_csv,
        records_from_enteric_json,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and enteric_result_paths:
        click.secho(
            "--csv and --enteric-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not enteric_result_paths:
        click.secho(
            "Provide input via --enteric-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_enteric_json(list(enteric_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_bov_enteric_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Bovine neonatal enteric trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = (
        f"{'Sample ID':<{col_id}}  {'Date':<12}  "
        f"{'ETEC':>5} {'SALM':>5} {'CRPT':>5} {'BCOV':>5} {'ROTA':>5}  Alert"
    )
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess_enteric(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.etec_k99 else '[-]':>5} "
            f"{'[+]' if r.flags.salmonella else '[-]':>5} "
            f"{'[+]' if r.flags.crypto else '[-]':>5} "
            f"{'[+]' if r.flags.bcov else '[-]':>5} "
            f"{'[+]' if r.flags.rota_a else '[-]':>5}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.etec_newly_detected:
        click.secho(
            "  ETEC K99 NEWLY DETECTED: IV antibiotic therapy required within 2 h; "
            "isolate affected calves immediately.",
            fg="red",
            bold=True,
        )
    if trend.etec_cleared:
        click.secho(
            "  ETEC K99 cleared: current treatment is working. Maintain hygiene and "
            "pen-sanitation protocols for at least three consecutive negative intervals.",
            fg="green",
        )
    if trend.salmonella_persistent:
        click.secho(
            "  SALMONELLA PERSISTENT: detected in every interval -- herd-level "
            "contamination or chronic shedder suspected; submit for culture and "
            "herd investigation.",
            fg="red",
        )
    if trend.zoonotic_risk_count > 0:
        click.secho(
            f"  Zoonotic risk in {trend.zoonotic_risk_count} interval(s) "
            "(Salmonella and/or Cryptosporidium): enforce gloves and handwashing.",
            fg="yellow",
        )
    if trend.antibiotic_required_count > 0:
        click.secho(
            f"  Antibiotics indicated in {trend.antibiotic_required_count} interval(s) "
            "(ETEC K99 and/or Salmonella detected).",
            fg="yellow",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")

    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"CSV timeline written to {out_csv}")


@cli.command(name="poc-respiratory-risk")
@click.option(
    "--hrsv", "hrsv", is_flag=True, default=False, help="hRSV (N gene, pan-RSV-A/B) positive."
)
@click.option(
    "--iav", "iav", is_flag=True, default=False, help="Influenza A (M gene, pan-IAV) positive."
)
@click.option("--sars2", "sars2", is_flag=True, default=False, help="SARS-CoV-2 (N gene) positive.")
@click.option(
    "--input-json",
    "input_json",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Read pathogen flags from a JSON file instead of individual flags. "
        "Expected keys: hrsv, iav, sars2 (bool). "
        "Missing keys default to false."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write a flat key-value CSV summary to a file.",
)
def poc_respiratory_risk(
    hrsv: bool,
    iav: bool,
    sars2: bool,
    input_json: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Interpret a three-target respiratory RT-LAMP panel as a clinical alert.

    Resolves the 'which antiviral?' differential for acute lower respiratory
    illness by interpreting a BioVind-style portable panel of three RT-LAMP
    assays: human RSV (hRSV), Influenza A (pan-IAV), and SARS-CoV-2.

    Detection of IAV or SARS-CoV-2 triggers antiviral-indicated guidance;
    hRSV detection triggers RSV-specific guidance (nirsevimab for eligible
    infants; supportive care for adults). Co-infection of any two or three
    targets yields a HIGH alert with combined antiviral recommendations.

    This command covers the respiratory sub-panel described in Recipe 34 of
    the cookbook (human RSV via RT-LAMP, POC vertical).  For the full five-
    target POC panel (MTB, CDiff, SARS2, Flu A, GAS) use 'lamp-forge poc-risk'.

    Panel targets:
    \b
        --hrsv    Human RSV (N gene, pan-RSV-A/B)    RT-LAMP (Recipe 34)
        --iav     Influenza A (M gene, pan-IAV)       RT-LAMP (Recipe 12)
        --sars2   SARS-CoV-2 (N gene)                 RT-LAMP (Recipe 2)

    \b
    Example -- hRSV detected in an infant presenting with bronchiolitis:
        lamp-forge poc-respiratory-risk --hrsv

    \b
    Example -- IAV detected, 30 hours into illness (oseltamivir window open):
        lamp-forge poc-respiratory-risk --iav

    \b
    Example -- SARS-CoV-2 + hRSV co-infection in an elderly patient:
        lamp-forge poc-respiratory-risk --sars2 --hrsv

    \b
    Example -- full panel from a JSON flags file:
        lamp-forge poc-respiratory-risk \
          --input-json results/resp_flags.json \
          --out-json results/resp_assessment.json
    """
    import json as json_mod

    from lamp_forge.poc_respiratory_risk import (
        PocRespAlertLevel,
        PocRespFlags,
        assess_poc_resp_risk,
        flags_from_dict,
        write_assessment_csv,
        write_assessment_json,
    )

    if input_json is not None:
        with input_json.open(encoding="utf-8") as fh:
            raw: dict[str, object] = json_mod.load(fh)
        flags = flags_from_dict(raw)
    else:
        flags = PocRespFlags(hrsv=hrsv, iav=iav, sars2=sars2)

    assessment = assess_poc_resp_risk(flags)

    # --- Panel table ----------------------------------------------------------
    click.echo("Respiratory differential RT-LAMP panel results:")
    target_rows = [
        ("hRSV", "N gene", assessment.flags.hrsv),
        ("IAV", "M gene", assessment.flags.iav),
        ("SARS2", "N gene", assessment.flags.sars2),
    ]
    for label, gene, positive in target_rows:
        symbol = "+" if positive else "-"
        click.echo(f"  {label:<6} ({gene:<6})  [{symbol}]")

    click.echo("")

    # --- Alert level (colour-coded) ------------------------------------------
    level_color = {
        PocRespAlertLevel.HIGH: "red",
        PocRespAlertLevel.MODERATE: "yellow",
        PocRespAlertLevel.LOW: "yellow",
        PocRespAlertLevel.NEGATIVE: "green",
    }
    color = level_color[assessment.alert_level]
    click.secho(
        f"Alert level : {assessment.alert_level.value}  (score {assessment.alert_score}/100)",
        fg=color,
        bold=(assessment.alert_level is PocRespAlertLevel.HIGH),
    )

    if assessment.antiviral_indicated:
        click.secho(
            "  Antiviral treatment indicated -- check treatment window before prescribing.",
            fg="yellow",
        )

    if assessment.rsv_indicated:
        click.secho(
            "  hRSV detected: assess nirsevimab eligibility for paediatric patients.",
            fg="yellow",
        )

    if assessment.dual_infection:
        click.secho(
            "  Dual/triple respiratory viral infection -- combined clinical management required.",
            fg="red",
        )

    click.echo("")
    click.echo("Interpretation:")
    click.echo(f"  {assessment.interpretation}")
    click.echo("")
    click.echo("Recommended action:")
    click.echo(f"  {assessment.recommended_action}")

    if out_json is not None:
        write_assessment_json(assessment, out_json)
        click.echo(f"\nAssessment written to {out_json}")

    if out_csv is not None:
        write_assessment_csv(assessment, out_csv)
        click.echo(f"CSV summary written to {out_csv}")


@cli.command(name="poc-respiratory-trend")
@click.option(
    "--resp-result",
    "resp_result_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    metavar="PATH",
    help=(
        "JSON file produced by 'lamp-forge poc-respiratory-risk --out-json' for one "
        "monitoring interval. Repeat once per sample in chronological order (oldest "
        "first). Mutually exclusive with --csv."
    ),
)
@click.option(
    "--csv",
    "csv_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Monitoring-spreadsheet CSV with columns: sample_id, date (optional), "
        "hrsv, iav, sars2. Rows must be in chronological order (oldest first). "
        "Mutually exclusive with --resp-result."
    ),
)
@click.option(
    "--out-json",
    "out_json",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the full trend assessment to a JSON file.",
)
@click.option(
    "--out-csv",
    "out_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the chronological pathogen-flag timeline to a CSV file.",
)
def poc_respiratory_trend(
    resp_result_paths: tuple[Path, ...],
    csv_path: Path | None,
    out_json: Path | None,
    out_csv: Path | None,
) -> None:
    r"""Analyse respiratory wave trajectory across consecutive RT-LAMP panel results.

    Takes a chronological sequence of three-target respiratory RT-LAMP panel
    results (hRSV, Influenza A, SARS-CoV-2) and produces a clinic surveillance
    trend: EMERGING, RESOLVING, STABLE_CLEAR, STABLE_ENDEMIC, or
    INSUFFICIENT_DATA.

    Key events detected:

    \b
      - IAV newly detected -- influenza season onset at this site
      - Dual-wave active (IAV + SARS-CoV-2 co-dominant in recent samples)

    Input (pick one):

    \b
      --resp-result FILE  JSON from 'lamp-forge poc-respiratory-risk --out-json',
                          repeated once per monitoring interval, oldest first.
      --csv FILE          Monitoring-spreadsheet CSV (header required):
                            sample_id, date (opt), hrsv, iav, sars2

    \b
    Example -- three-month respiratory wave track at a POC clinic:
        lamp-forge poc-respiratory-trend \\
          --resp-result results/clinic-A/2026-01.json \\
          --resp-result results/clinic-A/2026-02.json \\
          --resp-result results/clinic-A/2026-03.json \\
          --out-json results/clinic-A/resp_trend_Q1.json

    \b
    Example -- from a monitoring-spreadsheet CSV:
        lamp-forge poc-respiratory-trend \\
          --csv monitoring/clinic_A_resp_2026.csv \\
          --out-json results/clinic-A/resp_trend_2026.json
    """
    from lamp_forge.poc_respiratory_risk import assess_poc_resp_risk as _assess_resp_risk
    from lamp_forge.poc_respiratory_trend import (
        analyse_poc_resp_trend,
        records_from_csv,
        records_from_resp_json,
        write_trend_csv,
        write_trend_json,
    )

    if csv_path is not None and resp_result_paths:
        click.secho(
            "--csv and --resp-result are mutually exclusive. Use one or the other.",
            fg="red",
            err=True,
        )
        sys.exit(2)
    if csv_path is None and not resp_result_paths:
        click.secho(
            "Provide input via --resp-result (one or more JSON files) or --csv.",
            fg="red",
            err=True,
        )
        sys.exit(2)

    try:
        if csv_path is not None:
            records = records_from_csv(csv_path)
        else:
            records = records_from_resp_json(list(resp_result_paths))
    except (ValueError, OSError) as exc:
        click.secho(f"Input error: {exc}", fg="red", err=True)
        sys.exit(1)

    try:
        trend = analyse_poc_resp_trend(records)
    except ValueError as exc:
        click.secho(f"Analysis error: {exc}", fg="red", err=True)
        sys.exit(1)

    direction_colour = {
        "EMERGING": "red",
        "STABLE_CLEAR": "green",
        "STABLE_ENDEMIC": "yellow",
        "RESOLVING": "green",
        "INSUFFICIENT_DATA": "yellow",
    }
    col = direction_colour.get(trend.direction.value, "white")
    click.echo(f"Respiratory RT-LAMP panel surveillance trajectory: {trend.n_samples} sample(s)")
    click.echo("")

    col_id = max((len(r.sample_id) for r in trend.records), default=10)
    col_id = max(col_id, 9)
    header = f"{'Sample ID':<{col_id}}  {'Date':<12}  {'hRSV':>6} {'IAV':>6} {'SARS2':>6}  Alert"
    click.echo(header)
    click.echo("-" * len(header))

    for r in trend.records:
        lvl = _assess_resp_risk(r.flags).alert_level.value
        date_str = (r.date or "")[:12]
        click.echo(
            f"{r.sample_id:<{col_id}}  {date_str:<12}  "
            f"{'[+]' if r.flags.hrsv else '[-]':>6} "
            f"{'[+]' if r.flags.iav else '[-]':>6} "
            f"{'[+]' if r.flags.sars2 else '[-]':>6}  {lvl}"
        )
    click.echo("")

    click.secho(f"Trend direction: {trend.direction.value}", fg=col)

    if trend.iav_newly_detected:
        click.secho(
            "  IAV NEWLY DETECTED: influenza season onset -- check antiviral stock.",
            fg="yellow",
            bold=True,
        )
    if trend.dual_wave_active:
        click.secho(
            "  Dual-wave active: IAV + SARS-CoV-2 co-circulating in recent samples.",
            fg="yellow",
        )
    click.echo(f"  Worst alert level: {trend.worst_alert_level.value}")
    click.echo("")
    click.echo(f"Interpretation: {trend.interpretation}")
    click.echo("")
    click.echo(f"Recommended action: {trend.recommended_action}")

    if out_json is not None:
        write_trend_json(trend, out_json)
        click.echo(f"\nTrend assessment written to {out_json}")
    if out_csv is not None:
        write_trend_csv(trend, out_csv)
        click.echo(f"Timeline CSV written to {out_csv}")


if __name__ == "__main__":
    cli()
