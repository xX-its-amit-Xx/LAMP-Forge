"""Multiplex panel compatibility analysis.

LAMP-Forge designs one target at a time (``pipeline.run_pipeline``). A *panel*
runs several LAMP assays together in a single reaction on a single sample — the
mode portable platforms use to answer several questions from one swab or one
water sample at once.

The risk a multiplex introduces is **inter-assay cross-dimerisation**. LAMP is
primer-dense: six oligos per target, so a 3-plex shares 18 oligos in one tube
and a higher-plex panel proportionally more. Every pair of oligos drawn from
*different* assays is a chance for a heterodimer that quietly suppresses
amplification. Intra-set dimers are already vetted at design time
(:mod:`lamp_forge.primer_design`); this module screens the pairs that design
never sees — those *between* independently-designed assays.

Approach:

1. Load the top-ranked primer sets for each target from its ``primer_sets.json``.
2. For every pair of candidate sets belonging to different targets, compute the
   worst (most negative) heterodimer ΔG across all of their oligo pairs, reusing
   :func:`lamp_forge.primer_design.heterodimer_dg` (the primer3 engine the rest
   of the pipeline already uses).
3. Select one set per target so that the panel's *worst* inter-assay ΔG is as
   least-negative as possible (a max-min objective): the most thermodynamically
   independent combination. Exhaustive for small panels; a coordinate-ascent
   heuristic for large ones so the workflow still scales to high-plex panels.

The ΔG engine is injectable (``dg_func``) so the selection logic can be tested
without primer3, and so future work can swap in a different thermodynamic model.
"""

from __future__ import annotations

import itertools
import json
import logging
import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

logger = logging.getLogger(__name__)

# The serialised LampPrimerSet exposes its oligos as these named fields (the
# ``primers`` list is a property and is not present in the JSON). LF/LB may be
# absent (4-primer sets) and are skipped when missing.
_PRIMER_FIELDS = ("f3", "b3", "fip", "bip", "lf", "lb")

DimerFunc = Callable[[str, str], float]
"""Signature of a heterodimer ΔG function: ``(seq_a, seq_b) -> ΔG`` in kcal/mol."""


def _pair_key(uid_a: str, uid_b: str) -> tuple[str, str]:
    """Order-independent 2-tuple key for a pair of set uids."""
    return (uid_a, uid_b) if uid_a <= uid_b else (uid_b, uid_a)


# ----------------------------------------------------------------------------
# Data model
# ----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PanelOligo:
    """One synthesised oligo belonging to one candidate set.

    ``sequence`` is the actual ordered oligo in 5'→3' orientation — for FIP/BIP
    that is the assembled chimera, because the chimera is what is physically in
    the tube and what cross-dimerises.
    """

    target: str
    set_id: str
    role: str
    sequence: str


@dataclass(frozen=True, slots=True)
class PanelCandidateSet:
    """One ranked LAMP primer set for one target, loaded from ``primer_sets.json``."""

    target: str
    set_id: str
    composite_score: float
    oligos: tuple[PanelOligo, ...]

    @property
    def uid(self) -> str:
        """Globally unique id across all targets (``target::set_id``)."""
        return f"{self.target}::{self.set_id}"


@dataclass(frozen=True, slots=True)
class CrossDimer:
    """A heterodimer between two oligos belonging to *different* targets."""

    target_a: str
    set_id_a: str
    role_a: str
    sequence_a: str
    target_b: str
    set_id_b: str
    role_b: str
    sequence_b: str
    dg: float

    @property
    def label(self) -> str:
        """Human-readable ``TARGET/SET/ROLE × TARGET/SET/ROLE`` description."""
        return (
            f"{self.target_a}/{self.set_id_a}/{self.role_a} × "
            f"{self.target_b}/{self.set_id_b}/{self.role_b}"
        )


@dataclass(slots=True)
class PanelResult:
    """The selected multiplex combination and its cross-reactivity profile."""

    selection: dict[str, PanelCandidateSet]
    worst_dg: float
    flagged: list[CrossDimer]
    dimer_dg_threshold: float
    n_combinations_evaluated: int
    search_mode: str  # "exhaustive" | "greedy"

    @property
    def is_clean(self) -> bool:
        """True when no inter-assay pair breaches the ΔG threshold."""
        return not self.flagged

    @property
    def targets(self) -> list[str]:
        """Target names in selection order."""
        return list(self.selection)


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------


def load_candidate_sets(
    path: str | Path, target: str, *, top_n: int = 5
) -> list[PanelCandidateSet]:
    """Load the top-N ranked primer sets for one target from a ``primer_sets.json``.

    Sets are sorted by ``composite_score`` (descending) before truncating to
    ``top_n``, because the per-region design output is only locally sorted.

    Args:
        path: Path to a ``primer_sets.json`` written by :func:`report.write_json`.
        target: Short target label used to tag oligos in panel output.
        top_n: Keep at most this many top-scoring sets as panel candidates.

    Returns:
        Candidate sets for this target, best first.

    Raises:
        ValueError: If the file has no usable primer sets (e.g. every set is
            empty or the file is the wrong shape).
    """
    path = Path(path)
    with path.open() as fh:
        payload = json.load(fh)

    raw_sets = payload.get("primer_sets", [])
    if not isinstance(raw_sets, list):
        raise ValueError(f"{path}: 'primer_sets' is not a list")

    raw_sorted = sorted(raw_sets, key=lambda s: float(s.get("composite_score", 0.0)), reverse=True)

    candidates: list[PanelCandidateSet] = []
    for raw in raw_sorted[:top_n]:
        set_id = str(raw.get("set_id", f"{target}_unknown"))
        oligos: list[PanelOligo] = []
        for fieldname in _PRIMER_FIELDS:
            primer = raw.get(fieldname)
            if not primer:
                continue
            seq = primer.get("sequence")
            if not seq:
                continue
            oligos.append(
                PanelOligo(
                    target=target,
                    set_id=set_id,
                    role=str(primer.get("role", fieldname.upper())),
                    sequence=str(seq).upper(),
                )
            )
        if not oligos:
            continue
        candidates.append(
            PanelCandidateSet(
                target=target,
                set_id=set_id,
                composite_score=float(raw.get("composite_score", 0.0)),
                oligos=tuple(oligos),
            )
        )

    if not candidates:
        raise ValueError(
            f"{path}: no usable primer sets for target '{target}'. "
            "Was the upstream run successful and did it produce primer sets?"
        )
    return candidates


# ----------------------------------------------------------------------------
# Cross-dimer scoring
# ----------------------------------------------------------------------------


def _default_dg_func() -> DimerFunc:
    """Return the primer3-backed heterodimer ΔG function (imported lazily).

    Kept lazy so importing this module — and unit-testing the selection logic
    with an injected ΔG function — does not require primer3.
    """
    from lamp_forge.primer_design import heterodimer_dg

    return heterodimer_dg


def _set_pair_worst(a: PanelCandidateSet, b: PanelCandidateSet, dg_func: DimerFunc) -> CrossDimer:
    """Most-negative heterodimer across all oligo pairs between two sets.

    The returned :class:`CrossDimer` records which oligo pair was worst, so the
    report can point the user at the specific offending oligos.
    """
    worst: CrossDimer | None = None
    for oa in a.oligos:
        for ob in b.oligos:
            dg = dg_func(oa.sequence, ob.sequence)
            if worst is None or dg < worst.dg:
                worst = CrossDimer(
                    target_a=oa.target,
                    set_id_a=oa.set_id,
                    role_a=oa.role,
                    sequence_a=oa.sequence,
                    target_b=ob.target,
                    set_id_b=ob.set_id,
                    role_b=ob.role,
                    sequence_b=ob.sequence,
                    dg=dg,
                )
    # a and b each always carry >= 1 oligo (load_candidate_sets guarantees it).
    assert worst is not None
    return worst


def _build_pair_cache(
    candidates_by_target: Mapping[str, list[PanelCandidateSet]], dg_func: DimerFunc
) -> dict[tuple[str, str], CrossDimer]:
    """Precompute worst cross-dimer for every inter-target candidate-set pair.

    Keyed by the sorted pair of set uids. This is the only stage that calls the
    (potentially expensive) ΔG function; selection then reads the cache.
    """
    cache: dict[tuple[str, str], CrossDimer] = {}
    targets = list(candidates_by_target)
    for i in range(len(targets)):
        for j in range(i + 1, len(targets)):
            for a in candidates_by_target[targets[i]]:
                for b in candidates_by_target[targets[j]]:
                    cache[_pair_key(a.uid, b.uid)] = _set_pair_worst(a, b, dg_func)
    return cache


def _evaluate(
    selection: Mapping[str, PanelCandidateSet],
    cache: dict[tuple[str, str], CrossDimer],
    threshold: float,
) -> tuple[float, list[CrossDimer]]:
    """Worst inter-assay ΔG and flagged pairs for one one-set-per-target choice."""
    chosen = list(selection.values())
    worst = math.inf
    flagged: list[CrossDimer] = []
    for i in range(len(chosen)):
        for j in range(i + 1, len(chosen)):
            cd = cache[_pair_key(chosen[i].uid, chosen[j].uid)]
            worst = min(worst, cd.dg)
            if cd.dg < threshold:
                flagged.append(cd)
    flagged.sort(key=lambda c: c.dg)
    return worst, flagged


# ----------------------------------------------------------------------------
# Selection
# ----------------------------------------------------------------------------


def select_panel(
    candidates_by_target: Mapping[str, list[PanelCandidateSet]],
    *,
    dimer_dg_threshold: float = -5.0,
    dg_func: DimerFunc | None = None,
    max_combinations: int = 50_000,
) -> PanelResult:
    """Choose one primer set per target to maximise inter-assay independence.

    Objective is max-min: pick the combination whose *worst* inter-assay
    heterodimer ΔG is the least negative. Exhaustive search is used when the
    number of combinations is within ``max_combinations``; otherwise a
    coordinate-ascent heuristic (start from each target's top-scoring set, then
    swap one target at a time while it improves the global worst) keeps the
    search tractable for high-plex panels.

    Args:
        candidates_by_target: Mapping of target label → its candidate sets.
        dimer_dg_threshold: ΔG (kcal/mol) below which a pair is flagged. Matches
            the homodimer threshold convention used at design time.
        dg_func: Heterodimer ΔG function; defaults to the primer3-backed one.
        max_combinations: Switch to the heuristic above this product of
            per-target candidate counts.

    Returns:
        A :class:`PanelResult`.

    Raises:
        ValueError: If fewer than two targets are supplied.
    """
    targets = list(candidates_by_target)
    if len(targets) < 2:
        raise ValueError("A multiplex panel needs at least 2 targets")
    if dg_func is None:
        dg_func = _default_dg_func()

    cache = _build_pair_cache(candidates_by_target, dg_func)
    n_comb = math.prod(len(candidates_by_target[t]) for t in targets)

    if n_comb <= max_combinations:
        return _select_exhaustive(candidates_by_target, cache, dimer_dg_threshold)
    logger.info(
        "Panel has %d combinations (> %d); using greedy coordinate-ascent search",
        n_comb,
        max_combinations,
    )
    return _select_greedy(candidates_by_target, cache, dimer_dg_threshold)


def _select_exhaustive(
    candidates_by_target: Mapping[str, list[PanelCandidateSet]],
    cache: dict[tuple[str, str], CrossDimer],
    threshold: float,
) -> PanelResult:
    """Evaluate every one-set-per-target combination and keep the max-min best."""
    targets = list(candidates_by_target)
    best_selection: dict[str, PanelCandidateSet] | None = None
    best_flagged: list[CrossDimer] = []
    best_worst = -math.inf
    evaluated = 0
    for combo in itertools.product(*(candidates_by_target[t] for t in targets)):
        selection = {c.target: c for c in combo}
        worst, flagged = _evaluate(selection, cache, threshold)
        evaluated += 1
        if worst > best_worst:
            best_worst = worst
            best_selection = selection
            best_flagged = flagged
    assert best_selection is not None
    return PanelResult(
        selection=best_selection,
        worst_dg=best_worst,
        flagged=best_flagged,
        dimer_dg_threshold=threshold,
        n_combinations_evaluated=evaluated,
        search_mode="exhaustive",
    )


def _select_greedy(
    candidates_by_target: Mapping[str, list[PanelCandidateSet]],
    cache: dict[tuple[str, str], CrossDimer],
    threshold: float,
    *,
    max_iterations: int = 1000,
) -> PanelResult:
    """Coordinate-ascent: from each target's best set, swap one target at a time."""
    targets = list(candidates_by_target)
    selection = {t: max(candidates_by_target[t], key=lambda c: c.composite_score) for t in targets}
    worst, flagged = _evaluate(selection, cache, threshold)
    evaluated = 1

    for _ in range(max_iterations):
        improved = False
        for t in targets:
            current = selection[t]
            for alt in candidates_by_target[t]:
                if alt.uid == current.uid:
                    continue
                selection[t] = alt
                cand_worst, cand_flagged = _evaluate(selection, cache, threshold)
                evaluated += 1
                if cand_worst > worst:
                    worst, flagged, current = cand_worst, cand_flagged, alt
                    improved = True
                else:
                    selection[t] = current
        if not improved:
            break

    return PanelResult(
        selection=selection,
        worst_dg=worst,
        flagged=flagged,
        dimer_dg_threshold=threshold,
        n_combinations_evaluated=evaluated,
        search_mode="greedy",
    )


# ----------------------------------------------------------------------------
# Selected-panel matrix (for report + JSON)
# ----------------------------------------------------------------------------


def selected_oligos(result: PanelResult) -> list[PanelOligo]:
    """Flatten the chosen panel's oligos in target order."""
    oligos: list[PanelOligo] = []
    for target in result.selection:
        oligos.extend(result.selection[target].oligos)
    return oligos


def cross_dimer_matrix(
    oligos: list[PanelOligo], dg_func: DimerFunc | None = None
) -> list[list[float]]:
    """Full N×N heterodimer ΔG matrix among a flat oligo list.

    Diagonal and intra-set (same target) cells are set to ``nan`` — only
    inter-assay pairs are meaningful for panel screening.
    """
    if dg_func is None:
        dg_func = _default_dg_func()
    n = len(oligos)
    matrix = [[math.nan] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if oligos[i].target == oligos[j].target:
                continue
            dg = dg_func(oligos[i].sequence, oligos[j].sequence)
            matrix[i][j] = dg
            matrix[j][i] = dg
    return matrix


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------


def write_panel_json(result: PanelResult, path: Path) -> None:
    """Persist the panel selection, flagged cross-dimers, and metadata to JSON."""
    from lamp_forge import __version__

    payload = {
        "version": __version__,
        "generated_at": datetime.now(UTC).isoformat(),
        "search_mode": result.search_mode,
        "n_combinations_evaluated": result.n_combinations_evaluated,
        "dimer_dg_threshold": result.dimer_dg_threshold,
        "worst_inter_assay_dg": result.worst_dg,
        "is_clean": result.is_clean,
        "selection": [
            {
                "target": target,
                "set_id": cset.set_id,
                "composite_score": cset.composite_score,
                "oligos": [{"role": o.role, "sequence": o.sequence} for o in cset.oligos],
            }
            for target, cset in result.selection.items()
        ],
        "flagged_cross_dimers": [
            {
                "pair": cd.label,
                "target_a": cd.target_a,
                "role_a": cd.role_a,
                "sequence_a": cd.sequence_a,
                "target_b": cd.target_b,
                "role_b": cd.role_b,
                "sequence_b": cd.sequence_b,
                "dg_kcal_mol": cd.dg,
            }
            for cd in result.flagged
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(payload, fh, indent=2)
    logger.info("Wrote panel JSON (%d targets) to %s", len(result.selection), path)


def write_panel_csv(result: PanelResult, path: Path) -> None:
    """Order-ready CSV for the selected panel: one row per oligo, target-tagged."""
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["target", "set_id", "role", "primer_name", "sequence", "length"])
        for target, cset in result.selection.items():
            for o in cset.oligos:
                writer.writerow(
                    [target, cset.set_id, o.role, f"{target}_{o.role}", o.sequence, len(o.sequence)]
                )


_PANEL_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>LAMP-Forge multiplex panel — {{ targets|join(" + ") }}</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #1f2937; line-height: 1.5; }
  h1 { border-bottom: 3px solid #0369a1; padding-bottom: 0.3em; }
  h2 { margin-top: 2em; color: #0369a1; }
  .verdict { padding: 0.8em 1em; border-radius: 6px; font-weight: bold; margin: 1em 0; }
  .clean { background: #dcfce7; border-left: 4px solid #15803d; color: #14532d; }
  .dirty { background: #fee2e2; border-left: 4px solid #b91c1c; color: #7f1d1d; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 0.88em; }
  th, td { border: 1px solid #cbd5e1; padding: 6px 10px; text-align: left; }
  th { background: #e0f2fe; color: #075985; }
  tr:nth-child(even) { background: #f8fafc; }
  .seq { font-family: ui-monospace, Consolas, monospace; font-size: 0.85em; }
  .plot { margin: 1em 0; text-align: center; }
  .plot img { max-width: 100%; border: 1px solid #e2e8f0; border-radius: 4px; }
  .footer { margin-top: 3em; padding-top: 1em; border-top: 1px solid #e2e8f0;
            font-size: 0.85em; color: #64748b; }
</style>
</head>
<body>
<h1>LAMP-Forge multiplex panel report</h1>

<div class="verdict {{ 'clean' if is_clean else 'dirty' }}">
  {% if is_clean %}
  ✓ Compatible panel — no inter-assay heterodimer breaches the ΔG threshold
  ({{ "%.1f"|format(threshold) }} kcal/mol). Worst inter-assay ΔG:
  {{ "%.2f"|format(worst_dg) }} kcal/mol.
  {% else %}
  ✗ {{ flagged|length }} inter-assay heterodimer(s) breach the threshold
  ({{ "%.1f"|format(threshold) }} kcal/mol). Worst inter-assay ΔG:
  {{ "%.2f"|format(worst_dg) }} kcal/mol. Consider splitting the panel across
  tubes or re-designing the offending assay in a different conserved window.
  {% endif %}
</div>

<p>
  Selected one set per target via <strong>{{ search_mode }}</strong> search
  ({{ n_evaluated }} combination(s) evaluated). The heatmap shows pairwise
  heterodimer ΔG among all panel oligos; intra-assay cells are blank (vetted at
  design time).
</p>

<div class="plot"><img src="{{ heatmap }}" alt="cross-dimer heatmap"></div>

<h2>Selected panel</h2>
<table>
  <thead><tr><th>Target</th><th>Set</th><th>Composite</th><th>Role</th>
             <th>Sequence (5'→3')</th><th>Length</th></tr></thead>
  <tbody>
  {% for row in oligo_rows %}
    <tr>
      <td>{{ row.target }}</td><td>{{ row.set_id }}</td>
      <td>{{ "%.3f"|format(row.score) }}</td><td><strong>{{ row.role }}</strong></td>
      <td class="seq">{{ row.sequence }}</td><td>{{ row.length }}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>

{% if flagged %}
<h2>Flagged inter-assay heterodimers</h2>
<table>
  <thead><tr><th>Pair</th><th>ΔG (kcal/mol)</th></tr></thead>
  <tbody>
  {% for cd in flagged %}
    <tr><td>{{ cd.label }}</td><td>{{ "%.2f"|format(cd.dg) }}</td></tr>
  {% endfor %}
  </tbody>
</table>
{% endif %}

<div class="footer">
  Generated by LAMP-Forge v{{ version }}. In silico cross-dimer screen — confirm
  multiplex performance empirically (equimolar pooling, no-template control,
  per-target LoD) before deploying a panel.
</div>
</body>
</html>
"""


def render_panel_html(result: PanelResult, path: Path, dg_func: DimerFunc | None = None) -> Path:
    """Render the panel HTML report (cross-dimer heatmap + selected panel table)."""
    import base64
    import io

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from jinja2 import BaseLoader, Environment

    from lamp_forge import __version__

    oligos = selected_oligos(result)
    matrix = cross_dimer_matrix(oligos, dg_func)
    labels = [f"{o.target}/{o.role}" for o in oligos]

    arr = np.array(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(max(6, len(oligos) * 0.45), max(5, len(oligos) * 0.45)))
    masked = np.ma.masked_invalid(arr)
    cmap = plt.get_cmap("RdYlGn")
    im = ax.imshow(masked, cmap=cmap, vmin=result.dimer_dg_threshold * 2, vmax=0.0, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("Inter-assay heterodimer ΔG (kcal/mol)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("ΔG (more negative = stronger dimer)")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    plt.close(fig)
    heatmap = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")

    oligo_rows = [
        {
            "target": o.target,
            "set_id": o.set_id,
            "score": result.selection[o.target].composite_score,
            "role": o.role,
            "sequence": o.sequence,
            "length": len(o.sequence),
        }
        for o in oligos
    ]

    env = Environment(loader=BaseLoader(), autoescape=True, trim_blocks=True, lstrip_blocks=True)
    html = env.from_string(_PANEL_HTML_TEMPLATE).render(
        targets=result.targets,
        is_clean=result.is_clean,
        threshold=result.dimer_dg_threshold,
        worst_dg=result.worst_dg,
        search_mode=result.search_mode,
        n_evaluated=result.n_combinations_evaluated,
        heatmap=heatmap,
        oligo_rows=oligo_rows,
        flagged=result.flagged,
        version=__version__,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    logger.info("Panel HTML report written to %s", path)
    return path


def run_panel(
    set_specs: Mapping[str, Path],
    out_dir: Path,
    *,
    top_per_target: int = 5,
    dimer_dg_threshold: float = -5.0,
    generate_html: bool = True,
    dg_func: DimerFunc | None = None,
) -> PanelResult:
    """Load per-target sets, select a compatible panel, and write all outputs.

    Mirrors :func:`lamp_forge.pipeline.run_pipeline` so the CLI stays thin.

    Args:
        set_specs: Mapping of target label → path to that target's
            ``primer_sets.json``.
        out_dir: Directory for ``panel.json`` / ``panel_primers.csv`` /
            ``panel_report.html``.
        top_per_target: Candidate sets to consider per target.
        dimer_dg_threshold: Flagging threshold (kcal/mol).
        generate_html: Whether to render the HTML report.
        dg_func: Heterodimer ΔG function; defaults to the primer3-backed one.

    Returns:
        The :class:`PanelResult`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    candidates_by_target = {
        target: load_candidate_sets(path, target, top_n=top_per_target)
        for target, path in set_specs.items()
    }
    for target, cands in candidates_by_target.items():
        logger.info("Target '%s': %d candidate set(s) loaded", target, len(cands))

    result = select_panel(
        candidates_by_target,
        dimer_dg_threshold=dimer_dg_threshold,
        dg_func=dg_func,
    )

    write_panel_json(result, out_dir / "panel.json")
    write_panel_csv(result, out_dir / "panel_primers.csv")
    if generate_html:
        render_panel_html(result, out_dir / "panel_report.html", dg_func)

    return result
