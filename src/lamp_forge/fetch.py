"""NCBI Entrez retrieval with on-disk caching.

Why we cache: Entrez throttles to 3 requests/sec without an API key (10/sec
with one). A typical run pulls 10-50 sequences, and during development the
same query is hit dozens of times. Caching by query hash means the pipeline
is reproducible offline once warmed.

Why not pyentrezpy / refgenconf: Biopython's ``Entrez`` is the canonical
binding, and pulling in another wrapper adds surface for the wet-lab user
to debug. The tradeoff is some hand-rolled pagination logic here.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from pathlib import Path
from typing import Any

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from lamp_forge.types import PipelineConfig

logger = logging.getLogger(__name__)


def _query_hash(*parts: str | int | None) -> str:
    """Stable 12-hex-char hash for a query — used as cache filename."""
    blob = "|".join(str(p) for p in parts if p is not None)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


def _cache_path(cache_dir: Path, key: str) -> Path:
    """Resolve a cache key to a FASTA file path under ``cache_dir``."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"fetch_{key}.fasta"


def _configure_entrez(email: str) -> None:
    """Set the module-level Entrez configuration NCBI requires."""
    Entrez.email = email
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key
        logger.debug("Using NCBI API key (rate limit: 10 req/sec)")
    else:
        logger.debug("No NCBI_API_KEY set (rate limit: 3 req/sec)")


def search_taxon(
    taxon_id: int,
    gene: str | None,
    max_sequences: int,
    email: str,
) -> list[str]:
    """Look up nucleotide accessions for ``taxon_id``.

    If ``gene`` is provided, we constrain the search to that gene (e.g. "rpoB",
    "16S ribosomal RNA"); otherwise we pull whole-genome assemblies up to
    ``max_sequences``.

    Args:
        taxon_id: NCBI taxonomy ID (e.g. 1773 for *M. tuberculosis*).
        gene: Optional gene name. None retrieves full assemblies.
        max_sequences: Hard cap on returned accessions.
        email: Required by NCBI Entrez for courtesy headers.

    Returns:
        List of nucleotide accession.version strings, up to ``max_sequences`` long.
    """
    _configure_entrez(email)

    if gene:
        term = f"txid{taxon_id}[Organism] AND {gene}[Gene Name]"
    else:
        term = f"txid{taxon_id}[Organism] AND complete genome[Title]"

    logger.info("Entrez esearch: %s (max %d)", term, max_sequences)

    handle = Entrez.esearch(
        db="nucleotide",
        term=term,
        retmax=max_sequences,
        sort="relevance",
    )
    try:
        result = Entrez.read(handle)
    finally:
        handle.close()

    ids: list[str] = list(result.get("IdList", []))
    logger.info("Entrez returned %d UIDs", len(ids))
    return ids


def fetch_sequences(
    accessions: list[str],
    email: str,
    cache_dir: Path,
    batch_size: int = 10,
) -> list[SeqRecord]:
    """Fetch FASTA records for a list of accessions, with caching.

    Args:
        accessions: NCBI nucleotide accession.version strings or UIDs.
        email: NCBI courtesy email.
        cache_dir: Where to persist FASTA records between runs.
        batch_size: How many records to pull per efetch call. NCBI accepts up to
            10 000 but we keep batches small so a partial failure loses less work.

    Returns:
        List of :class:`Bio.SeqRecord.SeqRecord`. Order matches ``accessions``
        where possible; missing accessions are dropped (logged at WARNING).
    """
    if not accessions:
        return []

    _configure_entrez(email)
    cache_key = _query_hash("fetch", *sorted(accessions))
    cache_file = _cache_path(cache_dir, cache_key)

    if cache_file.exists():
        logger.info("Cache hit: %s (%d accessions)", cache_file, len(accessions))
        return list(SeqIO.parse(cache_file, "fasta"))

    records: list[SeqRecord] = []
    for batch_start in range(0, len(accessions), batch_size):
        batch = accessions[batch_start : batch_start + batch_size]
        logger.info("Entrez efetch batch %d-%d", batch_start, batch_start + len(batch))
        handle = Entrez.efetch(
            db="nucleotide",
            id=",".join(batch),
            rettype="fasta",
            retmode="text",
        )
        try:
            records.extend(SeqIO.parse(handle, "fasta"))
        finally:
            handle.close()
        # Polite pause: NCBI guideline is ≤3 req/sec without API key.
        if "NCBI_API_KEY" not in os.environ:
            time.sleep(0.35)

    SeqIO.write(records, cache_file, "fasta")
    logger.info("Cached %d records to %s", len(records), cache_file)
    return records


def fetch_for_config(config: PipelineConfig) -> list[SeqRecord]:
    """High-level entry: resolve config → list of SeqRecord, with caching.

    Resolution order:
      1. If ``config.accessions`` is non-empty, fetch those directly.
      2. Otherwise look up ``config.taxon_id`` (+ optional gene) via esearch,
         then efetch the resulting UIDs.

    Args:
        config: Parsed pipeline config.

    Returns:
        List of SeqRecord objects ready to feed into the alignment stage.

    Raises:
        RuntimeError: If neither accessions nor taxon-id-driven search yields
            any sequences (the pipeline cannot proceed without input).
    """
    if config.accessions:
        ids = config.accessions
    else:
        assert config.taxon_id is not None, "config validation should have caught this"
        ids = search_taxon(
            taxon_id=config.taxon_id,
            gene=config.gene,
            max_sequences=config.max_sequences,
            email=config.email,
        )

    if not ids:
        raise RuntimeError(
            f"No sequences returned for taxon_id={config.taxon_id}, gene={config.gene}. "
            "Try widening the gene filter or supplying explicit accessions."
        )

    records = fetch_sequences(ids, config.email, config.cache_dir)
    if not records:
        raise RuntimeError("Entrez returned UIDs but efetch yielded no FASTA records")
    return records


def manifest(records: list[SeqRecord]) -> dict[str, Any]:
    """Build a reproducibility manifest entry for fetched sequences.

    Returns a dict with the count and per-record id/length/md5. Used by
    ``run_manifest.json`` so a downstream consumer can detect when the input
    set changed without re-fetching.
    """
    entries = []
    for r in records:
        md5 = hashlib.md5(str(r.seq).encode()).hexdigest()
        entries.append({"id": r.id, "length": len(r.seq), "md5": md5})
    return {"count": len(records), "records": entries}


def write_fasta(records: list[SeqRecord], path: Path) -> None:
    """Write ``records`` to ``path`` as FASTA, creating parent dirs as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "fasta")
    logger.info("Wrote %d records to %s", len(records), path)


def load_fasta(path: Path) -> list[SeqRecord]:
    """Load a FASTA file from disk, returning a list of SeqRecords."""
    return list(SeqIO.parse(path, "fasta"))


def cached_or_fetch(path: Path, config: PipelineConfig) -> list[SeqRecord]:
    """Return records from ``path`` if it exists, else fetch and write them.

    This is the Snakemake-friendly entry point: if the target FASTA already
    exists (e.g. from a previous run), reuse it; otherwise hit Entrez.
    """
    if path.exists():
        return load_fasta(path)
    records = fetch_for_config(config)
    write_fasta(records, path)
    return records
