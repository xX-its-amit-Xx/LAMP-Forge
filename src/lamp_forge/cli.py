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
    out_csv: Path | None,
) -> None:
    r"""Estimate LAMP assay limit of detection (LOD) across the extraction chain.

    Computes the LOD in copies/reaction and back-calculates to copies/mL in
    the original sample using Poisson single-molecule statistics.

    Example (200 uL blood, 50% extraction into 50 uL eluate, 5 uL to reaction):

    \b
        lamp-forge lod --sample-volume 200 --eluate-volume 50 --reaction-input 5
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

    estimates = lod_table(params, tuple(probabilities))

    click.echo(
        f"Extraction chain: {sample_volume_ul:.0f} uL sample, "
        f"{extraction_efficiency * 100:.0f}% efficiency, "
        f"{eluate_volume_ul:.0f} uL eluate, "
        f"{reaction_input_ul:.1f} uL to reaction"
    )
    click.echo(
        f"Effective sample volume per reaction: {params.copies_per_rxn_per_copy_per_ul:.2f} uL"
    )
    click.echo("")
    header = f"{'P(detect)':>12}  {'lambda (copies/rxn)':>20}  {'LOD (copies/mL)':>18}"
    click.echo(header)
    click.echo("-" * len(header))
    for e in estimates:
        click.echo(
            f"{e.detection_probability:>12.3f}  "
            f"{e.lod_copies_per_reaction:>20.3f}  "
            f"{e.lod_copies_per_ml:>18.1f}"
        )

    if out_csv is not None:
        write_lod_csv(estimates, out_csv)
        click.echo(f"\nLOD table written to {out_csv}")


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


@cli.command(name="version")
def version() -> None:
    """Print version and exit."""
    click.echo(f"lamp-forge {__version__}")


if __name__ == "__main__":
    cli()
