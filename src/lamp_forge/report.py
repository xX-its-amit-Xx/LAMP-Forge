"""HTML + JSON + CSV report generation.

The HTML report is the primary deliverable for wet-lab users. It uses inline
matplotlib SVGs (no JS dependency) so the file is self-contained and can be
emailed around. JSON and CSV are for downstream automation.

Three plot types, all matplotlib → inline base64 PNG:

1. Conservation track: smoothed entropy across the alignment, with conserved
   regions highlighted.
2. Specificity heatmap: rows = primer sets, cols = off-target genomes, cells
   shaded by max identity × coverage.
3. Primer ladder: for the top N sets, show F3/F2/F1c/B1c/B2/B3 positions on
   the consensus.
"""

from __future__ import annotations

import base64
import csv
import io
import json
import logging
from dataclasses import asdict, is_dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from jinja2 import BaseLoader, Environment
from matplotlib.figure import Figure

from lamp_forge.conserve import ConservationTrack
from lamp_forge.types import ConservedRegion, LampPrimerSet, PipelineConfig

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------


def _fig_to_base64_png(fig: Figure) -> str:
    """Serialise a matplotlib figure to a base64 data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=130)
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def plot_conservation(
    track: ConservationTrack,
    regions: list[ConservedRegion],
    entropy_threshold: float,
) -> str:
    """Build a conservation-entropy track with conserved regions shaded.

    Args:
        track: Output of :func:`lamp_forge.conserve.compute_track`.
        regions: List of detected conserved regions (will be shaded green).
        entropy_threshold: Drawn as a horizontal reference line.

    Returns:
        Base64 PNG data URI suitable for embedding in HTML.
    """
    fig, ax = plt.subplots(figsize=(12, 3.5))
    x = np.arange(len(track.smoothed_entropy))
    ax.plot(
        x, track.smoothed_entropy, color="#0369a1", linewidth=0.8, label="smoothed entropy (bits)"
    )
    ax.axhline(
        entropy_threshold,
        color="#dc2626",
        linestyle="--",
        linewidth=0.8,
        label=f"threshold = {entropy_threshold:.2f}",
    )
    for region in regions:
        ax.axvspan(region.start, region.end, color="#86efac", alpha=0.35)
    ax.set_xlabel("Alignment position (bp)")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_title("Per-position conservation")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)

    if regions:
        patch = mpatches.Patch(color="#86efac", alpha=0.35, label="conserved region")
        handles, labels = ax.get_legend_handles_labels()
        handles.append(patch)
        labels.append("conserved region")
        ax.legend(handles, labels, loc="upper right", fontsize=8)

    return _fig_to_base64_png(fig)


def plot_specificity_heatmap(primer_sets: list[LampPrimerSet], top_n: int = 20) -> str:
    """Heatmap of max (identity × coverage) per (set, off-target).

    If no off-target hits exist (clean panel), returns an empty placeholder
    image rather than failing.
    """
    top_sets = primer_sets[:top_n]
    off_targets = sorted({h.off_target_name for s in top_sets for h in s.cross_reactivity_hits})
    if not off_targets or not top_sets:
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(
            0.5,
            0.5,
            "No off-target hits detected at any threshold.",
            ha="center",
            va="center",
            fontsize=11,
            color="#15803d",
        )
        ax.axis("off")
        return _fig_to_base64_png(fig)

    matrix = np.zeros((len(top_sets), len(off_targets)), dtype=np.float64)
    for i, s in enumerate(top_sets):
        for h in s.cross_reactivity_hits:
            j = off_targets.index(h.off_target_name)
            risk = h.percent_identity * h.coverage
            matrix[i, j] = max(matrix[i, j], risk)

    fig, ax = plt.subplots(figsize=(max(6, len(off_targets) * 0.6), max(4, len(top_sets) * 0.3)))
    im = ax.imshow(matrix, aspect="auto", cmap="Reds", vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(off_targets)))
    ax.set_xticklabels(off_targets, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(top_sets)))
    ax.set_yticklabels([s.set_id for s in top_sets], fontsize=7)
    ax.set_xlabel("Off-target genome")
    ax.set_ylabel("Primer set")
    ax.set_title(f"Cross-reactivity risk (top {len(top_sets)} sets)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("max(identity × coverage)")
    return _fig_to_base64_png(fig)


def plot_primer_ladder(primer_set: LampPrimerSet, region: ConservedRegion) -> str:
    """Draw F3/F2/F1c/B1c/B2/B3/LF/LB locations along the conserved region."""
    fig, ax = plt.subplots(figsize=(11, 2.5))
    ax.axhline(0, color="#475569", linewidth=1)
    ax.set_xlim(region.start, region.end)
    ax.set_ylim(-1.5, 1.5)

    role_color = {
        "F3": "#1e40af",
        "B3": "#1e40af",
        "FIP": "#7c3aed",
        "BIP": "#7c3aed",
        "LF": "#059669",
        "LB": "#059669",
    }
    for p in primer_set.primers:
        y = 0.5 if p.strand == "+" else -0.5
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (p.start, y - 0.18),
                p.length,
                0.36,
                boxstyle="round,pad=0.02",
                facecolor=role_color.get(p.role, "#64748b"),
                edgecolor="black",
                linewidth=0.5,
            )
        )
        ax.text(
            (p.start + p.end) / 2,
            y,
            p.role,
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            weight="bold",
        )
    ax.set_xlabel("Consensus position (bp)")
    ax.set_yticks([])
    ax.set_title(
        f"{primer_set.set_id} — mean Tm {primer_set.mean_tm:.1f}°C, "
        f"spread {primer_set.tm_spread:.1f}°C, F2-B2 {primer_set.f2_b2_distance} bp"
    )
    return _fig_to_base64_png(fig)


# ----------------------------------------------------------------------------
# Serialisation
# ----------------------------------------------------------------------------


def _to_json_safe(obj: Any) -> Any:
    """Recursively make an object JSON-serialisable (dataclasses, Paths, etc)."""
    if is_dataclass(obj) and not isinstance(obj, type):
        return {k: _to_json_safe(v) for k, v in asdict(obj).items()}
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (list, tuple)):
        return [_to_json_safe(v) for v in obj]
    if isinstance(obj, dict):
        return {k: _to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def write_json(primer_sets: list[LampPrimerSet], path: Path) -> None:
    """Persist the full primer-set list to JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": "0.1.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "primer_sets": [_to_json_safe(s) for s in primer_sets],
    }
    with path.open("w") as fh:
        json.dump(payload, fh, indent=2)
    logger.info("Wrote %d primer sets to %s", len(primer_sets), path)


def write_csv(primer_sets: list[LampPrimerSet], path: Path) -> None:
    """Order-ready CSV: one row per primer, with metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "set_id",
                "region_id",
                "primer_role",
                "primer_name",
                "sequence",
                "length",
                "tm_celsius",
                "gc_percent",
                "strand",
                "start",
                "end",
                "hairpin_dg_kcal_mol",
                "homodimer_dg_kcal_mol",
                "composite_score",
            ]
        )
        for s in primer_sets:
            for p in s.primers:
                writer.writerow(
                    [
                        s.set_id,
                        s.region_id,
                        p.role,
                        f"{s.set_id}_{p.role}",
                        p.sequence,
                        p.length,
                        f"{p.tm:.2f}",
                        f"{p.gc_percent:.1f}",
                        p.strand,
                        p.start,
                        p.end,
                        f"{p.hairpin_dg:.2f}",
                        f"{p.homodimer_dg:.2f}",
                        f"{s.composite_score:.4f}",
                    ]
                )


def write_conservation_tsv(track: ConservationTrack, path: Path) -> None:
    """Per-position entropy + consensus + coverage track."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["position", "raw_entropy", "smoothed_entropy", "consensus", "coverage"])
        for i in range(len(track.raw_entropy)):
            writer.writerow(
                [
                    i,
                    f"{track.raw_entropy[i]:.4f}",
                    f"{track.smoothed_entropy[i]:.4f}",
                    track.consensus[i],
                    f"{track.coverage[i]:.4f}",
                ]
            )


def write_manifest(
    config: PipelineConfig,
    input_manifest: dict[str, Any],
    n_regions: int,
    n_sets: int,
    path: Path,
) -> None:
    """Reproducibility manifest: config snapshot, input checksums, software versions."""
    from lamp_forge import __version__

    payload = {
        "lamp_forge_version": __version__,
        "generated_at": datetime.now(UTC).isoformat(),
        "config": _to_json_safe(asdict(config) if is_dataclass(config) else config),
        "input": input_manifest,
        "n_conserved_regions": n_regions,
        "n_primer_sets": n_sets,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(payload, fh, indent=2, default=str)


# ----------------------------------------------------------------------------
# HTML report
# ----------------------------------------------------------------------------


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>LAMP-Forge report — {{ target_name }}</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1200px; margin: 2em auto; padding: 0 1em; color: #1f2937;
         line-height: 1.5; }
  h1 { border-bottom: 3px solid #0369a1; padding-bottom: 0.3em; }
  h2 { margin-top: 2em; color: #0369a1; }
  h3 { color: #475569; }
  .meta { background: #f1f5f9; padding: 1em; border-radius: 6px; font-size: 0.9em; }
  .meta dt { font-weight: bold; display: inline-block; min-width: 180px; }
  .meta dd { display: inline; margin-left: 0; }
  .meta dd::after { content: "\A"; white-space: pre; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 0.88em; }
  th, td { border: 1px solid #cbd5e1; padding: 6px 10px; text-align: left; }
  th { background: #e0f2fe; color: #075985; }
  tr:nth-child(even) { background: #f8fafc; }
  .score { font-weight: bold; color: #15803d; }
  .score-low { color: #b45309; }
  .score-bad { color: #b91c1c; }
  .seq { font-family: ui-monospace, "SF Mono", Monaco, Consolas, monospace; font-size: 0.85em; }
  .plot { margin: 1em 0; text-align: center; }
  .plot img { max-width: 100%; border: 1px solid #e2e8f0; border-radius: 4px; }
  .warning { background: #fef3c7; border-left: 4px solid #f59e0b; padding: 0.5em 1em;
             margin: 0.5em 0; font-size: 0.9em; }
  .footer { margin-top: 3em; padding-top: 1em; border-top: 1px solid #e2e8f0;
            font-size: 0.85em; color: #64748b; }
  details { margin: 0.5em 0; }
  summary { cursor: pointer; font-weight: bold; color: #0369a1; }
</style>
</head>
<body>

<h1>LAMP-Forge primer design report</h1>

<div class="meta">
  <dt>Target:</dt><dd>{{ target_name }}</dd>
  <dt>NCBI taxon ID:</dt><dd>{{ taxon_id or "—" }}</dd>
  <dt>Gene filter:</dt><dd>{{ gene or "(whole-assembly search)" }}</dd>
  <dt>Sequences aligned:</dt><dd>{{ n_input_seqs }}</dd>
  <dt>Conserved regions found:</dt><dd>{{ n_regions }}</dd>
  <dt>Primer sets designed:</dt><dd>{{ n_sets }}</dd>
  <dt>Off-target panel:</dt><dd>{{ n_off_targets }} genome(s)</dd>
  <dt>Generated:</dt><dd>{{ generated_at }}</dd>
  <dt>LAMP-Forge version:</dt><dd>{{ version }}</dd>
</div>

<h2>1. Conservation analysis</h2>
<p>
  Smoothed Shannon entropy across the {{ alignment_width }} bp MSA. Shaded
  regions meet both the entropy threshold ({{ entropy_threshold }} bits) and
  the minimum length cutoff ({{ min_region_length }} bp); these are the
  candidate windows for LAMP primer placement.
</p>
<div class="plot"><img src="{{ conservation_plot }}" alt="conservation plot"></div>

<table>
  <thead>
    <tr><th>Region</th><th>Start</th><th>End</th><th>Length (bp)</th>
        <th>Mean entropy (bits)</th><th>Sets designed</th></tr>
  </thead>
  <tbody>
  {% for r in regions %}
    <tr>
      <td>{{ r.region_id }}</td>
      <td>{{ r.start }}</td>
      <td>{{ r.end }}</td>
      <td>{{ r.length }}</td>
      <td>{{ "%.3f"|format(r.mean_entropy) }}</td>
      <td>{{ sets_per_region[r.region_id] }}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>

<h2>2. Specificity screen</h2>
<p>
  Each primer was BLASTed against {{ n_off_targets }} user-supplied off-target
  genome(s). Cells in the heatmap below show <code>max(identity × coverage)</code>
  for each (primer set, off-target) pair. A set is flagged at any cell ≥
  {{ id_threshold }} × {{ cov_threshold }}.
</p>
<div class="plot"><img src="{{ specificity_plot }}" alt="specificity heatmap"></div>

<h2>3. Ranked primer sets (top {{ top_n }})</h2>
<p>
  Composite score = 0.35 × conservation + 0.45 × specificity + 0.20 × thermo balance.
</p>

<table>
  <thead>
    <tr>
      <th>Rank</th><th>Set ID</th><th>Region</th>
      <th>Composite</th><th>Conservation</th><th>Specificity</th><th>Thermo</th>
      <th>Mean Tm</th><th>Tm spread</th><th>F2-B2 (bp)</th><th>Loop primers</th>
    </tr>
  </thead>
  <tbody>
  {% for s in top_sets %}
    <tr>
      <td>{{ loop.index }}</td>
      <td>{{ s.set_id }}</td>
      <td>{{ s.region_id }}</td>
      <td class="score">{{ "%.3f"|format(s.composite_score) }}</td>
      <td>{{ "%.3f"|format(s.conservation_score) }}</td>
      <td>{{ "%.3f"|format(s.specificity_score) }}</td>
      <td>{{ "%.3f"|format(s.thermodynamic_score) }}</td>
      <td>{{ "%.1f"|format(s.mean_tm) }}°C</td>
      <td>{{ "%.1f"|format(s.tm_spread) }}°C</td>
      <td>{{ s.f2_b2_distance }}</td>
      <td>{{ "yes" if s.has_loop_primers else "no" }}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>

{% for entry in detailed %}
<details>
<summary>{{ entry.set.set_id }} — sequences and details</summary>
<div class="plot"><img src="{{ entry.ladder }}" alt="primer ladder"></div>
<table>
  <thead><tr><th>Role</th><th>Sequence (5'→3')</th><th>Length</th><th>Tm</th>
              <th>GC%</th><th>Strand</th><th>Hairpin ΔG</th><th>Dimer ΔG</th></tr></thead>
  <tbody>
  {% for p in entry.set.primers %}
    <tr>
      <td><strong>{{ p.role }}</strong></td>
      <td class="seq">{{ p.sequence }}</td>
      <td>{{ p.length }}</td>
      <td>{{ "%.1f"|format(p.tm) }}°C</td>
      <td>{{ "%.1f"|format(p.gc_percent) }}</td>
      <td>{{ p.strand }}</td>
      <td>{{ "%.2f"|format(p.hairpin_dg) }}</td>
      <td>{{ "%.2f"|format(p.homodimer_dg) }}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>
{% if entry.set.warnings %}
  {% for w in entry.set.warnings %}
    <div class="warning">{{ w }}</div>
  {% endfor %}
{% endif %}
</details>
{% endfor %}

<h2>4. Reproducibility</h2>
<p>
  Full config snapshot and input checksums are in <code>run_manifest.json</code>.
  Re-running LAMP-Forge with the same config + cached input will produce
  byte-identical primer sets. To pin the inputs further, copy the cached
  FASTAs out of <code>cache/</code> and reference them directly via
  <code>target.accessions</code> in your next config.
</p>

<div class="footer">
  Generated by LAMP-Forge v{{ version }}. <br>
  This is an in silico screen — primers must be validated in the wet lab
  before clinical or diagnostic use. See README for recommended validation
  steps and citations.
</div>

</body>
</html>
"""


def render_html(
    primer_sets: list[LampPrimerSet],
    regions: list[ConservedRegion],
    track: ConservationTrack,
    config: PipelineConfig,
    n_input_seqs: int,
    n_off_targets: int,
    path: Path,
    *,
    detail_count: int = 5,
) -> Path:
    """Render the full HTML report and write it to ``path``.

    ``detail_count`` controls how many of the top-ranked sets get an expandable
    detail card (sequences + ladder plot). Keep this small (~5) — each ladder
    plot adds ~30KB to the HTML and they're individually inspectable from the
    CSV/JSON otherwise.

    Returns the path written to.
    """
    from lamp_forge import __version__

    env = Environment(loader=BaseLoader(), autoescape=False, trim_blocks=True, lstrip_blocks=True)
    template = env.from_string(HTML_TEMPLATE)

    conservation_plot = plot_conservation(track, regions, config.entropy_threshold)
    specificity_plot = plot_specificity_heatmap(primer_sets, top_n=config.top_n)

    region_lookup = {r.region_id: r for r in regions}
    top_sets = primer_sets[: config.top_n]
    detailed = []
    for s in top_sets[:detail_count]:
        region = region_lookup.get(s.region_id)
        if region is None:
            continue
        detailed.append({"set": s, "ladder": plot_primer_ladder(s, region)})

    sets_per_region: dict[str, int] = {}
    for s in primer_sets:
        sets_per_region[s.region_id] = sets_per_region.get(s.region_id, 0) + 1
    for r in regions:
        sets_per_region.setdefault(r.region_id, 0)

    html = template.render(
        target_name=config.target_name,
        taxon_id=config.taxon_id,
        gene=config.gene,
        n_input_seqs=n_input_seqs,
        n_regions=len(regions),
        n_sets=len(primer_sets),
        n_off_targets=n_off_targets,
        generated_at=datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"),
        version=__version__,
        alignment_width=len(track.smoothed_entropy),
        entropy_threshold=config.entropy_threshold,
        min_region_length=config.min_region_length,
        regions=regions,
        sets_per_region=sets_per_region,
        id_threshold=f"{config.min_identity_threshold:.0%}",
        cov_threshold=f"{config.min_coverage_threshold:.0%}",
        top_n=min(config.top_n, len(primer_sets)),
        top_sets=top_sets,
        detailed=detailed,
        conservation_plot=conservation_plot,
        specificity_plot=specificity_plot,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    logger.info("HTML report written to %s", path)
    return path
