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
    type=click.Choice(["idt", "twist"], case_sensitive=False),
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
    r"""Export designed primers to an IDT or Twist vendor order CSV.

    Reads the primer_sets.json produced by 'lamp-forge run' and writes a
    spreadsheet ready for direct upload to the IDT SDS or Twist oligo-order
    portal.

    Scale and purification default to role-specific values:
    FIP/BIP get 100nm/HPLC (chimeric oligos benefit from full HPLC
    purification); F3/B3/LF/LB get 25nm/STD.

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
    type=click.Choice(["idt", "twist"], case_sensitive=False),
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
    selected primer set for every target into one IDT or Twist order sheet --
    removing the need to run 'lamp-forge export' separately for each target.

    Each primer name is prefixed with its target label so rows are
    self-documenting in the vendor portal (e.g. SRB_region_01_set_001_FIP).

    Note: Tm values are not stored in panel.json; the Tm column in Twist
    output will show 0.0.  For per-primer Tm, use 'lamp-forge export
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


@cli.command(name="version")
def version() -> None:
    """Print version and exit."""
    click.echo(f"lamp-forge {__version__}")


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


if __name__ == "__main__":
    cli()
