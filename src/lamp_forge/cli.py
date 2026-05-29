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


@cli.command(name="version")
def version() -> None:
    """Print version and exit."""
    click.echo(f"lamp-forge {__version__}")


if __name__ == "__main__":
    cli()
