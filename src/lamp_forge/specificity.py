"""In silico specificity screening via local BLAST.

We deliberately do NOT hit NCBI BLAST web service — at the scale of one run
(dozens of primer sets × dozens of off-target genomes × six primers per set =
thousands of queries) that would be slow, rate-limited, and rude. Instead we
build a local nucleotide BLAST database from user-supplied off-target FASTAs.

The user is expected to curate the off-target panel: close relatives, expected
contaminants, host DNA (if a clinical sample). This shifts the responsibility
to where it belongs — the assay developer knows the application context
better than NCBI nt's heterogeneous noise.

Scoring:

- For each primer, run ``blastn-short`` (optimised for sequences <50 nt).
- A "flagged hit" exceeds the user-configured identity AND coverage thresholds
  (defaults ≥80% identity over ≥80% of primer length, after Wang et al. 2018's
  empirical threshold for PCR cross-amplification risk).
- Set-level specificity score = 1 - (count of flagged hits / max_possible),
  where max_possible scales with set size and panel size. A set with zero
  flagged hits scores 1.0; a set whose every primer hits every off-target
  scores 0.0.
"""

from __future__ import annotations

import csv
import io
import logging
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from lamp_forge.types import LampPrimerSet, PipelineConfig, Primer, SpecificityHit

logger = logging.getLogger(__name__)


class BlastNotFoundError(RuntimeError):
    """Raised when ``makeblastdb`` or ``blastn`` cannot be located on PATH."""


class BlastError(RuntimeError):
    """Raised when a BLAST subprocess exits non-zero."""


def _resolve_binary(name: str) -> str:
    """Locate a BLAST+ binary or raise with install instructions."""
    binary = shutil.which(name)
    if not binary:
        raise BlastNotFoundError(
            f"{name} not found on PATH. Install BLAST+ with:\n"
            "  Ubuntu/Debian:  sudo apt install ncbi-blast+\n"
            "  macOS (brew):   brew install blast\n"
            "  Or use the LAMP-Forge Docker image."
        )
    return binary


def discover_off_targets(directory: Path) -> list[Path]:
    """Find off-target FASTA files in a directory.

    Recognised extensions: ``.fa``, ``.fasta``, ``.fna``, ``.ffn``. Symlinks
    and gzip'd files are not handled here — decompress them first.
    """
    if not directory.exists():
        logger.warning("Off-target directory does not exist: %s", directory)
        return []
    suffixes = {".fa", ".fasta", ".fna", ".ffn"}
    files = sorted(p for p in directory.iterdir() if p.suffix.lower() in suffixes)
    logger.info("Found %d off-target FASTA file(s) in %s", len(files), directory)
    return files


def build_blast_db(
    fasta_files: list[Path],
    db_path: Path,
) -> Path:
    """Concatenate off-target FASTAs and build a single BLAST nucleotide DB.

    Args:
        fasta_files: Paths to user-supplied off-target FASTA files.
        db_path: Output DB path (the actual files will be ``db_path.{nhr,nin,nsq}``).

    Returns:
        The path passed in ``db_path`` (for chaining).

    Raises:
        BlastNotFoundError: If ``makeblastdb`` is not on PATH.
        BlastError: If makeblastdb exits non-zero.
        ValueError: If ``fasta_files`` is empty.
    """
    if not fasta_files:
        raise ValueError("Cannot build BLAST DB from empty fasta_files list")

    makeblastdb = _resolve_binary("makeblastdb")
    db_path.parent.mkdir(parents=True, exist_ok=True)

    # Concatenate all off-target FASTAs into one input. We rewrite headers to
    # include the source filename so cross-reactivity hits can be traced back
    # to which genome they came from.
    concatenated = db_path.with_suffix(".combined.fasta")
    with concatenated.open("w") as out:
        for f in fasta_files:
            stem = f.stem
            for rec in SeqIO.parse(f, "fasta"):
                # Sanitize the original ID, then prefix with source genome.
                # ``|`` is a BLAST-recognised separator so we keep it.
                rec.id = f"{stem}|{rec.id}"
                rec.description = ""
                SeqIO.write([rec], out, "fasta")

    cmd = [
        makeblastdb,
        "-in",
        str(concatenated),
        "-dbtype",
        "nucl",
        "-out",
        str(db_path),
        "-parse_seqids",
    ]
    logger.info("Building BLAST DB: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, check=False)
    if proc.returncode != 0:
        raise BlastError(
            f"makeblastdb failed (exit {proc.returncode}):\n"
            f"stdout: {proc.stdout.decode(errors='replace')}\n"
            f"stderr: {proc.stderr.decode(errors='replace')}"
        )

    return db_path


def _blastn_short(
    primers: list[Primer],
    db_path: Path,
    *,
    word_size: int = 7,
    evalue: float = 1.0,
) -> list[dict[str, str]]:
    """Run blastn-short against ``db_path`` for the supplied primers.

    Outputs a flat list of dicts (one per HSP). Empty query results just
    aren't included — no hits means no rows for that primer.

    Args:
        primers: Primers to query.
        db_path: Path to a BLAST nucleotide DB (built by :func:`build_blast_db`).
        word_size: BLAST seed word size. 7 is the blastn-short default and
            works well for 18-22 nt primers.
        evalue: E-value cutoff. We keep this loose (1.0) because primer-length
            queries naturally have high e-values; we filter by identity + coverage
            downstream, not by e-value.

    Returns:
        List of dicts with keys: qseqid, sseqid, pident, length, qlen,
        evalue, bitscore. All values are strings as BLAST emits them.
    """
    blastn = _resolve_binary("blastn")

    query_fasta = io.StringIO()
    # Use the primer role + start as the query id so BLAST output traces back
    # cleanly. Duplicate-sequence primers (FIP/BIP arms) would collide on
    # role-only ids.
    role_lookup: dict[str, Primer] = {}
    for i, p in enumerate(primers):
        qid = f"{p.role}_{i}"
        role_lookup[qid] = p
        query_fasta.write(f">{qid}\n{p.sequence}\n")

    cmd = [
        blastn,
        "-task",
        "blastn-short",
        "-db",
        str(db_path),
        "-word_size",
        str(word_size),
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length qlen evalue bitscore",
        "-max_target_seqs",
        "50",
    ]
    logger.debug("Running %s on %d primers", " ".join(cmd), len(primers))
    proc = subprocess.run(
        cmd,
        input=query_fasta.getvalue().encode(),
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise BlastError(
            f"blastn failed (exit {proc.returncode}):\n"
            f"stderr: {proc.stderr.decode(errors='replace')}"
        )

    rows: list[dict[str, str]] = []
    reader = csv.DictReader(
        io.StringIO(proc.stdout.decode()),
        fieldnames=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "evalue",
            "bitscore",
        ],
        delimiter="\t",
    )
    for row in reader:
        row["_primer"] = role_lookup[row["qseqid"]].sequence  # type: ignore[index]
        row["_role"] = role_lookup[row["qseqid"]].role  # type: ignore[index]
        rows.append(row)
    return rows


def _row_to_hit(row: dict[str, str], cfg: PipelineConfig) -> SpecificityHit:
    """Convert a raw BLAST tabular row into a typed :class:`SpecificityHit`."""
    pident = float(row["pident"]) / 100.0
    aln_len = int(row["length"])
    qlen = int(row["qlen"])
    coverage = aln_len / qlen if qlen else 0.0
    off_target = row["sseqid"].split("|")[0]
    return SpecificityHit(
        primer_role=row["_role"],  # type: ignore[arg-type]
        primer_sequence=row["_primer"],
        off_target_name=off_target,
        percent_identity=pident,
        alignment_length=aln_len,
        primer_length=qlen,
        coverage=coverage,
        e_value=float(row["evalue"]),
        bitscore=float(row["bitscore"]),
        is_flagged=(
            pident >= cfg.min_identity_threshold and coverage >= cfg.min_coverage_threshold
        ),
    )


def screen_set(
    primer_set: LampPrimerSet,
    db_path: Path,
    cfg: PipelineConfig,
) -> LampPrimerSet:
    """BLAST every primer in a set, fill cross_reactivity_hits + specificity_score.

    Mutates and returns ``primer_set``. Mutation is intentional — the design
    stage created a partially scored set, and this stage finishes it.
    """
    rows = _blastn_short(primer_set.primers, db_path)
    hits = [_row_to_hit(r, cfg) for r in rows]
    primer_set.cross_reactivity_hits = hits

    flagged = [h for h in hits if h.is_flagged]
    n_primers = len(primer_set.primers)
    # Worst case: every primer hits every off-target. We measure penalty
    # against the count of unique (primer, off_target) flagged pairs so
    # multiple HSPs to the same genome don't double-penalise.
    flagged_pairs = {(h.primer_role, h.off_target_name) for h in flagged}
    # Count distinct off-target genomes seen at all (flagged or not) — that's
    # the size of the panel.
    panel_size = len({h.off_target_name for h in hits}) or 1
    max_penalty = n_primers * panel_size
    score = 1.0 - (len(flagged_pairs) / max_penalty) if max_penalty else 1.0
    primer_set.specificity_score = max(0.0, score)

    if flagged_pairs:
        primer_set.warnings.append(
            f"{len(flagged_pairs)} primer→off-target pair(s) crossed the "
            f"{cfg.min_identity_threshold:.0%} identity / "
            f"{cfg.min_coverage_threshold:.0%} coverage threshold"
        )

    return primer_set


def screen_all(
    primer_sets: list[LampPrimerSet],
    db_path: Path,
    cfg: PipelineConfig,
) -> list[LampPrimerSet]:
    """Run specificity screening across all primer sets, re-sort by composite."""
    from lamp_forge.primer_design import _composite_score

    for s in primer_sets:
        screen_set(s, db_path, cfg)
        s.composite_score = _composite_score(s)

    primer_sets.sort(key=lambda s: s.composite_score, reverse=True)
    return primer_sets


def write_specificity_tsv(primer_sets: Iterable[LampPrimerSet], path: Path) -> None:
    """Write a flat TSV of every (primer, off-target) hit across all sets."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "set_id",
                "primer_role",
                "primer_sequence",
                "off_target",
                "percent_identity",
                "coverage",
                "alignment_length",
                "primer_length",
                "evalue",
                "bitscore",
                "is_flagged",
            ]
        )
        for ps in primer_sets:
            for h in ps.cross_reactivity_hits:
                writer.writerow(
                    [
                        ps.set_id,
                        h.primer_role,
                        h.primer_sequence,
                        h.off_target_name,
                        f"{h.percent_identity:.4f}",
                        f"{h.coverage:.4f}",
                        h.alignment_length,
                        h.primer_length,
                        f"{h.e_value:.2e}",
                        f"{h.bitscore:.1f}",
                        int(h.is_flagged),
                    ]
                )


def _records_have_content(records: list[SeqRecord]) -> bool:
    """Sanity check: are any of the parsed records non-empty?"""
    return any(len(r.seq) > 0 for r in records)
