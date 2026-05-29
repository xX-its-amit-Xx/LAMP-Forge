"""MAFFT wrapper for multiple sequence alignment.

We shell out to MAFFT rather than using Biopython's AlignIO front-end so the
user can swap in `--auto`, `--einsi`, or `--retree` strategies without touching
Python. MAFFT is bundled in the Docker image; for local installs we surface a
clear error if the binary is missing.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

logger = logging.getLogger(__name__)


class MafftNotFoundError(RuntimeError):
    """Raised when the ``mafft`` binary cannot be located on PATH."""


class MafftError(RuntimeError):
    """Raised when MAFFT exits non-zero. Message includes captured stderr."""


def _resolve_mafft() -> str:
    """Locate the ``mafft`` binary or raise a friendly error."""
    binary = shutil.which("mafft")
    if not binary:
        raise MafftNotFoundError(
            "MAFFT not found on PATH. Install with:\n"
            "  Ubuntu/Debian:  sudo apt install mafft\n"
            "  macOS (brew):   brew install mafft\n"
            "  Or use the LAMP-Forge Docker image, which bundles MAFFT."
        )
    return binary


def run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    *,
    strategy: str = "--auto",
    threads: int = 1,
    extra_args: list[str] | None = None,
) -> Path:
    """Run MAFFT on ``input_fasta`` and write the alignment to ``output_fasta``.

    Args:
        input_fasta: Path to unaligned FASTA file.
        output_fasta: Path to write the aligned FASTA.
        strategy: MAFFT strategy flag. Defaults to ``--auto`` which picks
            FFT-NS-2 / FFT-NS-i / L-INS-i based on input size. Override with
            ``--einsi`` for higher-quality alignments of <200 sequences with
            local similarity (good for marker genes), or ``--retree 2`` for
            speed on large inputs.
        threads: How many CPU threads MAFFT may use. ``--thread -1`` uses all.
        extra_args: Additional CLI args appended verbatim. Use for
            ``--reorder``, ``--adjustdirection``, etc.

    Returns:
        Path to the alignment (``output_fasta``).

    Raises:
        MafftNotFoundError: If the binary is not on PATH.
        MafftError: If MAFFT exits non-zero.
        FileNotFoundError: If ``input_fasta`` does not exist.
    """
    if not input_fasta.exists():
        raise FileNotFoundError(f"Alignment input not found: {input_fasta}")

    mafft = _resolve_mafft()
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [
        mafft,
        strategy,
        "--thread",
        str(threads),
        "--quiet",
        str(input_fasta),
    ]
    if extra_args:
        # Insert before the input path so they apply.
        cmd = [cmd[0], *cmd[1:-1], *extra_args, cmd[-1]]

    logger.info("Running MAFFT: %s", " ".join(cmd))
    with output_fasta.open("w") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=False)

    if proc.returncode != 0:
        stderr = proc.stderr.decode("utf-8", errors="replace")
        raise MafftError(
            f"MAFFT failed with exit code {proc.returncode}.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr:\n{stderr}"
        )

    logger.info("Alignment written to %s", output_fasta)
    return output_fasta


def load_alignment(path: Path) -> MultipleSeqAlignment:
    """Read an aligned FASTA file into a Biopython MultipleSeqAlignment.

    Args:
        path: Path to aligned FASTA.

    Returns:
        Biopython :class:`MultipleSeqAlignment` ready for entropy computation.

    Raises:
        ValueError: If the file has fewer than 2 sequences or contains
            sequences of unequal length (i.e. not an alignment).
    """
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) < 2:
        raise ValueError(
            f"Alignment file {path} has only {len(records)} sequence(s); "
            "need ≥2 for conservation analysis"
        )

    lengths = {len(r.seq) for r in records}
    if len(lengths) > 1:
        raise ValueError(
            f"Alignment file {path} contains sequences of unequal length "
            f"({sorted(lengths)}); not a valid MSA"
        )
    return MultipleSeqAlignment(records)


def align_records(
    input_fasta: Path,
    output_fasta: Path,
    *,
    strategy: str = "--auto",
    threads: int = 1,
) -> MultipleSeqAlignment:
    """Convenience: run MAFFT and load the resulting alignment in one call.

    Equivalent to :func:`run_mafft` followed by :func:`load_alignment`.
    """
    run_mafft(input_fasta, output_fasta, strategy=strategy, threads=threads)
    return load_alignment(output_fasta)
