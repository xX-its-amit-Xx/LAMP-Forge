# LAMP-Forge workflow — the science behind each step

This document is the long-form companion to [README.md](../README.md). The
README tells you *how* to run LAMP-Forge; this tells you *why* each step
exists, what tradeoffs are baked in, and where you'd most plausibly want
to adjust the defaults.

## Why six primers? A two-minute LAMP refresher

A standard PCR reaction uses two primers, exponentially amplifying a template
through repeated thermal cycling. LAMP works at a single temperature (60–65°C)
and uses **six primers** binding to **eight regions** because that's what's
needed to self-prime a strand-displacing polymerase into making the
characteristic cauliflower-like double-loop DNA structure that gives the
reaction its exponential character — without cycling.

```
Template (3' → 5'):

  ────F3c────  ────F2c────  ────F1c────             ────B1────  ────B2────  ────B3────
       ↑           ↑              ↑                      ↓           ↓           ↓
      F3          F2  F1 (RC, in FIP)              B1c(RC,in BIP)    B2          B3
                                       ─────loop─────
                                          LF / LB
```

- **F3 and B3** are the outer displacement primers. They start the reaction.
- **FIP** (F1c + F2) folds onto itself once extended, kicking off the
  loop architecture.
- **BIP** (B1c + B2) does the same on the other side.
- **LF and LB** sit in the single-stranded loop regions formed during the
  reaction and accelerate amplification ~2-3×. They're optional but
  recommended.

If you've never seen a LAMP diagram, the [Eiken animation](https://loopamp.eiken.co.jp/e/lamp/principle.html)
is the clearest 5-minute primer.

## The pipeline, step by step

### Step 1 — Fetch ([src/lamp_forge/fetch.py](../src/lamp_forge/fetch.py))

We hit NCBI Entrez through Biopython. Two modes:

1. **Taxonomy-driven:** Give us a `taxon_id` and optionally a `gene` name,
   we esearch + efetch the top `max_sequences` by relevance.
2. **Explicit accessions:** Give us a list of accession.version strings, we
   efetch directly. Use this for full reproducibility — `relevance` ordering
   can shift as NCBI updates.

**Why we cache.** Without an `NCBI_API_KEY` the rate limit is 3 requests/sec.
A whole-pipeline dev loop hits the same query dozens of times. The cache key
is a SHA256 of the sorted accession list, so the cache hit is deterministic.

**The big tradeoff.** When `gene` is set, Entrez does keyword matching, not
sequence-based matching. If your gene of interest is annotated inconsistently
(e.g. some isolates call it `rpoB`, others `RpoB`, others by locus tag), you
will miss sequences. For canonical marker genes this is fine; for
less-studied genes, prefer the `accessions` route and curate the list yourself.

### Step 2 — Align ([src/lamp_forge/align.py](../src/lamp_forge/align.py))

We shell out to MAFFT in `--auto` mode. MAFFT picks between FFT-NS-2,
FFT-NS-i, and L-INS-i based on input size and similarity. For 10-50
marker-gene sequences this lands on FFT-NS-i, which is a good balance of
speed and accuracy.

**When to override.** For highly divergent sequences (cross-genus
alignments), use `--einsi`. For very large inputs (>500 sequences), `--retree 2`
trades alignment quality for runtime. Both are CLI overrides:

```bash
lamp-forge align --input raw.fa --output aligned.fa --strategy --einsi
```

**What we don't do.** No alignment trimming with Gblocks/TrimAl. The
conservation step's coverage filter (drop columns with >20% gaps) handles the
"edges are noisy" problem without throwing away potentially conserved
regions. If your sequences are wildly variable lengths and you want to trim
first, do it as a pre-processing step and pass the trimmed FASTA in via
`target.accessions: [path/to/trimmed.fa]` — fetch will detect existing
sequences and skip Entrez.

### Step 3 — Conserve ([src/lamp_forge/conserve.py](../src/lamp_forge/conserve.py))

This is where the choice of "conserved" gets pinned down.

For each alignment column we compute Shannon entropy:

H(col) = − Σ p(b) · log₂ p(b)   for b ∈ {A, C, G, T}

with gaps and ambiguous bases excluded from the sum. Then we smooth across a
window the size of a typical primer (default 30 nt) — this turns
single-column flukes into smooth tracks where "is there a primer-sized
patch around here that's conserved?" is a one-glance answer.

**Why Shannon entropy specifically.** Three reasons. (1) Scales naturally to
>2 sequences. (2) Treats partial conservation correctly: a column that's
80% A / 20% G is meaningfully less conserved than one that's 99% A / 1% G,
and entropy captures that linearly. (3) The marker-gene literature reports
conservation in bits (Schneider & Stephens 1990 / Schneider's "sequence logo"
work), so the threshold values translate directly from published surveys.

**The defaults.** `entropy_threshold: 0.25` bits is on the strict side —
appropriate for marker genes (rpoB, 16S rRNA V regions) and bacterial
housekeeping genes. For viral genomes or hypervariable surface antigens,
loosen to 0.5 bits or drop the window size to 20 nt to find smaller
conserved patches. Tighten to 0.1 bits if you want to be very confident
that all isolates in your training set will amplify.

**Coverage filter.** Columns where >20% of the rows are gap or ambiguous get
excluded. This catches MSA edge artefacts (one short sequence with no flanking
data forces gaps in N-1 rows for the column). The threshold is configurable
via `min_coverage` in the function signature but not yet surfaced in the YAML
schema — open an issue if you need it.

### Step 4 — Primer design ([src/lamp_forge/primer_design.py](../src/lamp_forge/primer_design.py))

Geometric constraints come straight from Notomi 2000 and Eiken's "A Guide
to LAMP Primer Design":

| Constraint | Default | Why |
|---|---|---|
| F3/B3 length | 18-22 nt | Standard PCR primer band; outer displacement |
| F2/B2 length | 18-22 nt | Inner forward/reverse; must clear a loop |
| F1c/B1c length | 20-22 nt | The "looped-back" arm of FIP/BIP |
| Loop primer length | 18-22 nt | Bound to single-stranded loop regions |
| F3-F2 / B2-B3 gap | 0-20 bp | Lets F3 displace F2 product cleanly |
| F2-F1 / B1-B2 gap | 40-60 bp | Drives the loop size — too tight = no loop |
| F2-B2 inner amplicon | 120-160 bp | Total inner product; >200bp slows reaction |
| Tm | 60-65 °C | Bst polymerase optimum |
| Tm match across set | ≤ 2 °C spread | Otherwise primers fall out of step |
| GC | 40-65% | Standard primer chemistry |
| Hairpin ΔG | ≥ -2 kcal/mol | Stable hairpins reduce effective primer pool |
| Homodimer ΔG | ≥ -5 kcal/mol | Same logic; dimers compete with template |

We compute Tm and ΔG via [primer3-py](https://pypi.org/project/primer3-py/),
which wraps Untergasser's primer3 C library — the same engine used by Primer-BLAST
and most synthesis-house design tools. The salt corrections default to LAMP
buffer-typical values (50 mM Na+, 8 mM Mg2+, 0.8 mM dNTP). If you're using
an exotic buffer, adjust the constants in `primer_design.tm()`.

**Why we generate-then-filter instead of using primer3's `designPrimers`.**
That orchestrator is PCR-specific. It doesn't understand FIP/BIP chimeras or
the F2-F1 / B1-B2 spacing constraints. So we generate every plausible
candidate ourselves, score with primer3's thermo primitives, and assemble
sets with explicit geometry checks. The combinatorics are bounded by the
length of the conserved region (a 250 bp region produces ~50-200 candidates
per role, well within reach).

### Step 5 — Specificity ([src/lamp_forge/specificity.py](../src/lamp_forge/specificity.py))

We build a local BLAST nucleotide database from user-supplied off-target
FASTAs, then run `blastn-short` (word size 7, optimised for ≤50 nt queries)
for each primer.

**Why local BLAST, not NCBI Web BLAST.** At our query volume (six primers ×
dozens of sets × dozens of off-targets) hitting NCBI's server is slow,
rate-limited, and arguably impolite. The local DB approach also forces the
user to think about which off-targets actually matter for their application
— a deliberate design choice, not a limitation.

**Scoring.** A hit is **flagged** if `percent_identity ≥ min_identity_threshold`
AND `coverage ≥ min_coverage_threshold` (defaults 80% × 80%, after Wang et al.
2018's empirical floor for PCR cross-amplification risk). The set-level
specificity score = `1 - (flagged_pairs / max_possible_pairs)`. A set with no
flagged hits at any primer scores 1.0; a set whose every primer hits every
off-target genome scores 0.0.

**Curating the off-target panel.** For a clinical diagnostic, include:

- Close relatives within the same genus (the species you most need to
  *not* amplify).
- Sister species commonly mistaken for your target in your sample type.
- The host genome (or a representative fragment) if the sample is host-derived.
- Common contaminants of your sample matrix (skin commensals for blood,
  reagent contaminants like *Pseudomonas* for low-biomass samples).

Don't dump the entire `nt` database in — it'll dilute your scoring and
flag hits to organisms that don't exist in your sample type.

### Step 6 — Report ([src/lamp_forge/report.py](../src/lamp_forge/report.py))

Three outputs:

- **HTML report:** Self-contained file (inline matplotlib PNGs, no JS). Open
  it locally, email it, attach to a notebook. Show your PI without spinning
  up a server.
- **JSON:** Complete primer set list with all metadata. For downstream
  automation, e.g. piping into an ordering API.
- **CSV:** One row per primer. The columns map directly to what IDT, Twist,
  and Eurofins ordering forms ask for.

**Composite score weighting.** The ranking is
`0.35 × conservation + 0.45 × specificity + 0.20 × thermodynamic_balance`.
Specificity dominates because false positives are the failure mode that kills
clinical utility. Conservation comes second because if the primer doesn't
match your isolate, it can't amplify. Thermo is necessary but the hard
filters already enforce the floor, so weight is lower. Override these weights
by editing `_composite_score` in [primer_design.py](../src/lamp_forge/primer_design.py)
if your use case differs (e.g. environmental surveillance might weight
conservation higher than specificity).

## Reproducibility

Each run writes a `run_manifest.json` containing:

- LAMP-Forge version
- Full config snapshot
- Input record count + per-record MD5
- Counts of regions found and sets designed

Re-running with the same config + cached FASTAs produces byte-identical
output. For long-term archival of a particular design, copy the cached FASTA
out of `cache/` and reference it via `target.accessions: [<file path>]` —
your design will then be pinned to a specific input even if NCBI's
relevance ranking shifts.

## See also

- [lamp_primer_design.md](lamp_primer_design.md) — deeper dive on LAMP
  primer geometry with diagrams.
- [cookbook.md](cookbook.md) — real-world recipes for specific organisms
  and sample types.
