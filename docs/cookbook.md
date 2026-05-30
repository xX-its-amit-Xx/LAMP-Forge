# LAMP-Forge cookbook — real-world recipes

Worked recipes for designing LAMP assays against real targets, with the
config decisions explained and the off-target panel choices justified.
Run any of these by copying the YAML block into `config/<name>.yaml`,
dropping the listed off-target FASTAs into `input/off_targets/`, and
invoking:

```bash
docker compose run --rm lamp-forge run --config /work/config/<name>.yaml
```

## Recipe 1 — *Mycobacterium tuberculosis* (TB) diagnostic

**Goal.** A LAMP assay that detects *M. tuberculosis complex* (MTBC) DNA
from a sputum sample, without amplifying environmental or commensal
mycobacteria.

**Why it's interesting.** TB is a top-five global cause of death. Existing
LAMP-based TB diagnostics (e.g. Eiken's TB-LAMP) target the gyrB or IS6110
loci. We'll target rpoB — the gene that also harbours rifampicin-resistance
mutations, so the assay output positions you for a follow-up resistance
test on the same amplicon.

**Target sequences.** NCBI taxon ID 1773 (*M. tuberculosis*), gene `rpoB`,
max 20 sequences pulls a diverse clinical isolate set including reference
strains H37Rv, CDC1551, and Erdman, plus contemporary outbreak isolates.

**Off-target panel.** Drop these FASTAs in `input/off_targets/`:

| File | Source | Why |
|---|---|---|
| `m_smegmatis.fasta` | NC_008596.1 | Environmental mycobacterium, lab contaminant |
| `m_avium.fasta` | NC_002944.2 | Non-tuberculous mycobacterium (NTM), most common clinical confounder |
| `m_intracellulare.fasta` | NC_016946.1 | Second-most-common NTM |
| `m_bovis.fasta` | NC_002945.4 | Within MTBC — tests whether your assay can distinguish at the species level (most TB-LAMP assays accept this as on-target) |
| `human_chr_fragment.fasta` | Subset of GRCh38 chr1 | Host DNA contamination from sputum |

**Config.**

```yaml
target:
  name: mtb_rpoB
  taxon_id: 1773
  gene: rpoB
  max_sequences: 25
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.20    # stricter: clinical sensitivity matters
  min_region_length: 250

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 50.0               # rpoB is GC-rich in mycobacteria
  gc_max: 70.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 130
    f2_b2_max: 160

output:
  dir: results/mtb_rpoB
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Expected output.** 6-12 conserved regions ≥250 bp (rpoB is well-conserved
across MTBC), 5-15 ranked primer sets per region. The top-scoring sets
should have 0 flagged off-target hits to *M. smegmatis*, *M. avium*, and
*M. intracellulare*. *M. bovis* will likely show some hits — that's
expected; MTBC is the conventional clinical target.

**Wet-lab validation note.** Sputum is a difficult matrix. Pre-treatment with
*N*-acetyl-L-cysteine or SDS before LAMP is standard; the assay should be
validated against both treated and untreated sputum.

---

## Recipe 2 — SARS-CoV-2 *N*-gene LAMP

**Goal.** A colorimetric LAMP assay for SARS-CoV-2 from nasal swab RNA,
robust across Omicron and post-Omicron lineages.

**Why it's interesting.** This is the assay archetype of the COVID-era
"open-source LAMP" community (Color Genomics, NEB, Zhang lab). Multiple
published primer sets exist; we're checking whether LAMP-Forge can
re-derive equivalent designs from current sequence data.

**Target sequences.** Use `accessions` here, not taxon-driven search — viral
genomes evolve fast and Entrez relevance ranking shifts. Pull a curated set
of recent Omicron subvariants (BA.5, XBB, JN.1, KP.3) and a pre-Omicron
anchor.

```yaml
target:
  name: sars_cov_2_N
  accessions:
    - OQ518516.1   # pre-Omicron WA-CDC reference (post-edit)
    - OR351212.1   # BA.5
    - OQ858867.1   # XBB.1.5
    - OZ016291.1   # JN.1
    - PP389033.1   # KP.3
    # ...add more recent lineages as they appear
  gene: null      # we're giving accessions directly; gene filter unused
  max_sequences: 25
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tighter — coronaviruses cross-react easily
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.30        # looser — RNA viruses have higher background
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0                   # SARS-CoV-2 is AT-leaning
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/sars_cov_2_N
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** Other respiratory coronaviruses and viruses that
co-circulate:

| File | Source | Why |
|---|---|---|
| `hcov_oc43.fasta` | NC_006213.1 | Common-cold coronavirus |
| `hcov_nl63.fasta` | NC_005831.2 | Common-cold coronavirus |
| `mers_cov.fasta` | NC_019843.3 | Same betacoronavirus genus |
| `sars_cov_1.fasta` | NC_004718.3 | Closest sister to SARS-CoV-2 — the specificity acid test |
| `influenza_a_h1n1.fasta` | NC_026438.1 (PB2) | Co-circulating respiratory virus |
| `rsv_a.fasta` | NC_038235.1 | Co-circulating respiratory virus |

**Expected output.** With tight thresholds, expect 3-6 viable sets. The top
set's *N*-gene primers should land in the highly conserved N-gene
3' region (the same region NEB's published colorimetric LAMP kit targets).
*Some* identity hits to SARS-CoV-1 are expected — flag them and decide based on
clinical context whether SARS-CoV-1 cross-reactivity is acceptable (in 2026, it
is).

---

## Recipe 3 — *Salmonella enterica* food-safety screen

**Goal.** Detect any *Salmonella enterica* serovar in food sample DNA
without amplifying common food-microbiome bacteria.

**Why it's interesting.** Foodborne *Salmonella* outbreaks happen
constantly. The challenge is broad inclusivity (~2500 serovars) without
flagging *E. coli*, *Klebsiella*, and other Enterobacteriaceae that share
food-microbiome niches.

**Target.** *invA* gene — the canonical *Salmonella* marker for PCR/LAMP
assays. Highly conserved across serovars, broadly absent in close relatives.

```yaml
target:
  name: salmonella_invA
  taxon_id: 28901   # Salmonella enterica
  gene: invA
  max_sequences: 30  # capture serovar diversity
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.20
  min_region_length: 250

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0
  gc_max: 65.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/salmonella_invA
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** Enterobacteriaceae that share food matrices:

| File | Source | Why |
|---|---|---|
| `e_coli_k12.fasta` | NC_000913.3 | Closest relative; food-microbiome staple |
| `e_coli_o157.fasta` | NC_002695.2 | Pathogenic E. coli — must not co-flag |
| `klebsiella_pneumoniae.fasta` | NC_016845.1 | Enterobacteriaceae, raw produce |
| `enterobacter_cloacae.fasta` | NC_014121.1 | Common food contaminant |
| `citrobacter_freundii.fasta` | NC_021907.1 | Enterobacteriaceae |
| `shigella_flexneri.fasta` | NC_004337.2 | Very close to E. coli, often co-tested |

**Expected output.** *invA* should yield clean separation — the canonical
hit zone for published invA-LAMP primers (e.g. Hara-Kudo et al. 2005)
appears as a long conserved region in the report. Top sets should have 0
flagged off-target hits even at the 80% × 80% threshold.

**Wet-lab note.** Food matrices contain inhibitors. Validate with a spiked
enrichment broth protocol (BAM Ch. 5 standard) and check LOD at 10
CFU/25g — that's the FDA action threshold for ready-to-eat foods.

---

## Recipe 4 — *Plasmodium falciparum* malaria diagnostic

**Goal.** Detect *P. falciparum* from finger-prick blood at field-deployable
sensitivity, without amplifying *P. vivax* or other Plasmodium species
(they require different treatment).

**Why it's interesting.** Malaria LAMP is one of the few diagnostic LAMP
applications that genuinely outperforms RDTs in field deployment (Patel
et al. 2014). Species-level specificity is critical: *P. vivax* causes
relapsing malaria requiring primaquine treatment, while *P. falciparum*
requires artemisinin combination therapy.

**Target.** *Pf*18S rRNA gene. Already used in WHO's MALACHITE LAMP kit;
we'll independently re-derive primers.

```yaml
target:
  name: pfalciparum_18S
  taxon_id: 5833   # Plasmodium falciparum
  gene: 18S        # case-sensitive; alternative: "small subunit ribosomal RNA"
  max_sequences: 20
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # 18S is highly conserved across Plasmodium
  min_coverage_threshold: 0.85

conservation:
  window_size: 30
  entropy_threshold: 0.15        # 18S is extremely conserved within species
  min_region_length: 230

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 30.0                   # Plasmodium is AT-rich
  gc_max: 55.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/pfalciparum_18S
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** Other *Plasmodium* species + human host:

| File | Source | Why |
|---|---|---|
| `p_vivax_18S.fasta` | NC_009906.1 region | Relapsing-malaria species, *must* differentiate |
| `p_ovale_18S.fasta` | XR_002677062.1 | Less common but co-endemic |
| `p_malariae_18S.fasta` | M54897.1 | Less common but co-endemic |
| `p_knowlesi_18S.fasta` | NC_011906.1 region | Zoonotic, Southeast Asia |
| `human_18S.fasta` | NR_003286.4 | Eukaryotic 18S background — major confounder |

**Expected output.** The *P. falciparum* 18S has a published "P.f.-specific
window" (~580-820 of the 18S coordinate system in PlasmoDB). Top LAMP-Forge
sets should localise to this window. *P. vivax* identity should fall to
<80% in the flagged region — if it doesn't, raise `min_identity_threshold`
and re-run.

---

## Recipe 5 — Antibiotic resistance: KPC-type carbapenemase

**Goal.** Detect KPC carbapenemase genes (any variant) in clinical
*Klebsiella* / *E. coli* isolates for AMR surveillance.

**Why it's interesting.** Resistance genes spread across genus boundaries
on plasmids, so you cannot use organism-level specificity — you need
*gene-level* specificity. The off-target panel here is unusual: it's
other β-lactamase families that aren't carbapenemases.

**Target.** `blaKPC` family (KPC-1 through KPC-50+ are all close variants).

```yaml
target:
  name: blaKPC_family
  taxon_id: 0                    # placeholder; not used because accessions are explicit
  accessions:
    - KP119087.1   # KPC-2 (the most common variant)
    - KU902437.1   # KPC-3
    - JN982240.1   # KPC-4
    - KX833570.1   # KPC-5
    # ...
  gene: null
  max_sequences: 30
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.85   # stricter coverage — point mutants matter

conservation:
  window_size: 25
  entropy_threshold: 0.20
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0
  gc_max: 65.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/blaKPC_family
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** Related but *non-carbapenemase* β-lactamases —
the assay must distinguish KPC from ESBLs that look superficially similar:

| File | Source | Why |
|---|---|---|
| `blaCTX-M-15.fasta` | DQ302097.1 | Most common ESBL, ≠ carbapenemase |
| `blaTEM-1.fasta` | AY458016.1 | Classic ESBL family |
| `blaSHV-12.fasta` | X98101.1 | ESBL family |
| `blaNDM-1.fasta` | FN396876.1 | A *different* carbapenemase — may want a separate assay |
| `blaVIM-2.fasta` | AF302086.1 | Another carbapenemase — separate target |
| `blaOXA-48.fasta` | JN626286.1 | OXA carbapenemase — separate target |

**Expected output.** Top sets should have 0 flagged hits against blaCTX-M,
blaTEM, blaSHV. Hits against blaNDM/VIM/OXA at <50% identity are expected
(they're more distantly related). The composite score will hover around
0.85-0.95 for well-designed sets.

**Caveats for AMR LAMP.** Resistance gene assays in clinical use are usually
*multiplexed* (KPC + NDM + VIM + OXA in one tube). LAMP-Forge's current
version designs one target at a time. Run the pipeline five times (once per
carbapenemase family) and validate the primer sets together for multiplex
compatibility (heterodimer ΔG across families) before deployment.

---

---

# Field-deployable recipes (portable, on-site, multiplex)

The recipes below target the settings LAMP was built for: results read off at a
single temperature, no thermocycler, no cold chain. They cover the three
environments where a portable isothermal box earns its keep — **industrial
asset integrity** (oil & gas), **on-farm animal biosecurity**, and
**multiplexed panels** that interrogate one sample for several organisms at
once. The multiplex recipe (Recipe 8) uses the `lamp-forge panel` command to
check that independently-designed assays can share a tube.

> **Accession note.** Where a recipe lists accessions marked *(representative)*,
> treat them as a starting point and verify against current NCBI before
> ordering — viral taxonomy and RefSeq accessions drift. Taxon-driven recipes
> (`taxon_id` + `gene`) are more robust to this and are preferred where the
> marker gene is annotated consistently.

## Recipe 6 — Sulfate-reducing bacteria (oilfield souring / MIC monitoring)

**Goal.** Detect sulfate-reducing bacteria (SRB) broadly across genera in
produced water, pipeline biofilm, and pigging solids — the primary drivers of
reservoir souring (H₂S generation) and microbiologically influenced corrosion
(MIC) in oil & gas infrastructure.

**Why it's interesting.** SRB are *polyphyletic* — they span
Deltaproteobacteria (*Desulfovibrio*, *Desulfobacter*), Firmicutes
(*Desulfotomaculum*), Thermodesulfobacteria, and even Archaea
(*Archaeoglobus*). You cannot get broad SRB coverage from an organism-level
(16S, single-genus) assay. The field-standard handle is the **functional gene
`dsrAB`** (dissimilatory sulfite reductase) — present in essentially all
sulfate reducers and the basis of the qPCR assays already used for oilfield
souring surveillance. This is a *gene-level inclusivity* problem, the mirror
image of the AMR recipe (Recipe 5): we want to catch every SRB, not exclude
relatives.

**The specificity trap.** Sulfur-*oxidizing* bacteria (e.g. *Chlorobium*,
some Proteobacteria) carry a **reverse-type `dsrAB`** that is homologous to the
SRB enzyme. A naive assay will light up on these and over-report corrosion
risk. They are the key off-target panel members below.

**Target.** `dsrB` (beta subunit — slightly more conserved primary structure
than `dsrA` for primer placement). Broad SRB coverage needs a *curated
multi-genus* reference set; the cleanest source is the FunGene `dsrB`
repository (fungene.cme.msu.edu) or an NCBI protein-guided nucleotide pull
across the genera below. Start from a dominant oilfield genus and widen:

```yaml
target:
  name: srb_dsrB_oilfield
  taxon_id: 872          # Desulfovibrio (anchor genus); widen with accessions
  accessions:
    # Curate a multi-genus dsrB set for true broad coverage — pull dsrB CDS for
    # Desulfovibrio, Desulfobacter, Desulfotomaculum, Desulfobulbus,
    # Archaeoglobus from NCBI/FunGene and list them here.
    []
  gene: dsrB
  max_sequences: 40      # SRB diversity is wide — sample it
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.35   # looser — dsrB is functionally, not sequence, conserved
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0
  gc_max: 65.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/srb_dsrB
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** Reverse-`dsr` organisms + non-SRB oilfield community:

| File | Source | Why |
|---|---|---|
| `chlorobium_tepidum_rdsr.fasta` | NC_002932.3 (*dsrAB* region) | Sulfur **oxidizer**, reverse-type *dsrAB* — the false-positive trap |
| `allochromatium_vinosum_rdsr.fasta` | reverse-*dsr* region | Purple sulfur bacterium, reverse-type *dsr* |
| `pseudomonas_aeruginosa.fasta` | NC_002516.2 | Ubiquitous non-SRB oilfield biofilm former |
| `thauera_nitrate_reducer.fasta` | representative NRB | Nitrate-reducing bacteria used to *out-compete* SRB in souring control — must not co-flag |
| `crude_oil_microbiome_mixed.fasta` | your site metagenome | Real background; best off-target you can supply |

**Expected output.** `dsrB` is a protein-coding gene under strong functional
constraint but with wobble-position variability, so expect conserved windows
that track codon structure. Top sets should show **no flagged hits** to the
reverse-`dsr` oxidizers at 80% × 80% — if they do, raise
`min_identity_threshold` and prefer windows in the most conserved catalytic
motifs. Broad SRB inclusivity is only as good as your input accession set;
widen `accessions` until added genera stop changing the conserved windows.

**Field note.** Produced water and crude carry severe inhibitors (high salt,
hydrocarbons, sulfide). Validate the chosen set on a spiked produced-water
matrix and pair with a simple DNA cleanup (e.g. magnetic-bead capture) before
the LAMP step. The win for a portable assay here is turnaround: souring
decisions (biocide dosing, pigging schedule) are made on-site, where shipping
samples to a lab costs days.

---

## Recipe 7 — African swine fever virus (on-farm biosecurity)

**Goal.** Detect African swine fever virus (ASFV) from whole blood, serum,
tissue, or oral fluid on the farm or at a border post — fast enough to inform
a quarantine decision before the animal moves.

**Why it's interesting.** ASFV is the highest-consequence pig pathogen in the
world: ~100% case fatality for the virulent genotype II strains driving the
current Eurasian/global panzootic, no widely deployed vaccine, and control that
depends entirely on early detection and culling. It is a large **double-stranded
DNA** virus (*Asfarviridae*) — so, unlike most farm-virus LAMP targets, it needs
**no reverse transcription**, making it an unusually clean isothermal,
field-deployable target. The clinical picture (high fever, hemorrhage) overlaps
classical swine fever, so the assay's job is as much *differential diagnosis*
as detection.

**Target.** `B646L`, the gene encoding the **p72** major capsid protein — the
WOAH/OIE-referenced PCR target and the basis of published ASFV LAMP assays. It
is well-conserved across genotypes, which is what you want for broad detection.

```yaml
target:
  name: asfv_p72
  taxon_id: 10497        # African swine fever virus
  gene: B646L            # p72 capsid; alt annotation "p72" / "major capsid protein"
  max_sequences: 25      # span genotypes I/II at minimum
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tighter — must distinguish from other swine viruses
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.25
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/asfv_p72
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** The swine viral/host differential:

| File | Source | Why |
|---|---|---|
| `csfv.fasta` | classical swine fever virus *(representative)* | The #1 clinical confounder — pestivirus, RNA, must not co-flag |
| `prrsv.fasta` | PRRS virus *(representative)* | Endemic respiratory/reproductive pathogen, co-circulates |
| `pcv2.fasta` | porcine circovirus 2 *(representative)* | Ubiquitous swine DNA virus — closest "DNA background" |
| `porcine_parvovirus.fasta` | PPV *(representative)* | Common swine DNA virus |
| `sus_scrofa_fragment.fasta` | subset of *Sus scrofa* genome | Host DNA — the dominant background in blood/tissue |

**Expected output.** p72 yields several broadly-conserved windows. Top sets
should have zero flagged hits across the entire panel — ASFV is genetically
distant from the pestivirus/arterivirus confounders, so clean separation is the
expected (and required) result. If any porcine DNA-virus hits appear, tighten
to 90% identity.

**Field note.** This is a textbook BioVind-style deployment: a portable box at
a farm or checkpoint, blood or oral-fluid in, quarantine decision out in under
an hour. Because ASFV is a reportable disease, a positive must be confirmed by
the national reference lab — position this assay as a **triage/screening** tool,
not the confirmatory test.

---

## Recipe 8 — Multiplex panel: one sample, several organisms (`lamp-forge panel`)

**Goal.** Build a small multiplexed isothermal panel — several LAMP assays that
run together in one reaction on one sample — and verify the primer sets won't
cross-react (heterodimerise) with each other before you commit to ordering.

**Why it's interesting.** Portable diagnostics earn their value when one swab or
one water sample answers several questions at once. But LAMP is primer-dense
(six primers per target), so a 3-plex is 18 oligos sharing a tube and a 5-plex
is 30 — every inter-assay pair is a chance for a cross-dimer that quietly kills
sensitivity. LAMP-Forge designs *one* target at a time; the **`panel`** command
closes the loop by screening combinations of independently-designed sets for
inter-assay compatibility (minimum cross-dimer ΔG) and proposing the most
compatible one-set-per-target combination.

**Worked example — oilfield souring panel.** Three functional guilds that
together describe corrosion/souring risk in a produced-water sample:

1. **SRB** via `dsrB` (Recipe 6) — sulfide generation / MIC.
2. **Methanogens** via `mcrA` (methyl-coenzyme M reductase) — methanogenic
   corrosion and reservoir gas.
3. **General bacterial load** via a conserved `16S rRNA` window — total burden /
   sample adequacy control.

Run the pipeline once per guild, then combine:

```bash
# 1. Design each assay independently (three configs, three output dirs).
docker compose run --rm lamp-forge run --config /work/config/srb_dsrB.yaml
docker compose run --rm lamp-forge run --config /work/config/methanogen_mcrA.yaml
docker compose run --rm lamp-forge run --config /work/config/general_16s.yaml

# 2. Check cross-assay compatibility and propose a panel.
docker compose run --rm lamp-forge panel \
  --set SRB=/work/results/srb_dsrB/primer_sets.json \
  --set MCR=/work/results/methanogen_mcrA/primer_sets.json \
  --set 16S=/work/results/general_16s/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out /work/results/souring_panel
```

**What `panel` does.** For each target it considers the top *N* ranked sets
(`--top-per-target`), then searches combinations (one set per target) for the
one whose **worst inter-assay heterodimer ΔG** is least negative (i.e. the most
thermodynamically independent combination). It writes:

| File | What it is |
|---|---|
| `panel_report.html` | Cross-dimer heatmap + the recommended compatible combination |
| `panel.json` | Full pairwise ΔG matrix and ranked combinations, for downstream tooling |
| `panel_primers.csv` | Order-ready flat list for the recommended panel, target-tagged |

**Reading the output.** Any inter-assay primer pair with ΔG below
`--dimer-dg-threshold` is flagged. A clean panel has *no* flagged pairs; if every
combination flags, your targets are too primer-dense to co-multiplex and should
be split across two tubes (or one target re-designed in a different conserved
window). This is exactly the manual "validate the sets together for multiplex
compatibility" step that Recipe 5 told you to do by hand — now automated.

**Scaling note.** BioVind-class platforms read out up to ~18 targets per sample;
the same workflow scales — design each target, then run one `panel` call across
all of them. Combinatorial search is capped by `--top-per-target`, so keep that
small (3–5) for large panels.

---

## Recipe 9 — Estimating assay LOD before ordering primers (`lamp-forge lod`)

**Goal.** Before committing to ordering and validating a primer set, estimate
the analytical limit of detection (LOD) the assay can achieve given your
specific sample type and DNA/RNA extraction protocol.

**Why it matters for field deployment.** BioVind-style portable platforms
process small sample volumes and use rapid extraction kits optimised for speed
over yield — the effective sample volume reaching the reaction is a small
fraction of what was collected. LOD_95 (the copy number at which 95% of
reactions are positive) is the standard regulatory benchmark for a *detection*
assay. Knowing your LOD before ordering lets you decide:

- Is 300 copies/mL sufficient sensitivity for the clinical or industrial
  threshold? (TB in sputum: ~100–1000 GE/mL; souring SRB in produced water:
  often >10³ cells/mL.)
- Would increasing sample input volume or improving extraction efficiency push
  the LOD below the required threshold without needing a more sensitive assay?

**Model.** The `lod` command uses the **Poisson single-molecule** model:

```
P(detect | mean copies = λ) = 1 - exp(-λ)
LOD_95: λ = -ln(0.05) ≈ 2.996 copies / reaction
```

Back-calculation to the original sample:

```
LOD [copies/mL] = λ_LOD / (V_sample × η × V_rxn_input / V_eluate) × 1000
```

**Usage.** The `lamp-forge lod` command takes extraction-chain parameters and
prints an LOD table. All volumes are in microlitres (uL).

```bash
# Produced-water SRB screen (Recipe 6):
# 1 mL produced water, 40% bead-extraction efficiency,
# eluted in 100 uL buffer, 5 uL added to 25 uL LAMP reaction.
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.40 \
  --eluate-volume 100 \
  --reaction-input 5 \
  --out-csv results/srb_lod.csv
```

Expected output:

```
Extraction chain: 1000 uL sample, 40% efficiency, 100 uL eluate, 5.0 uL to reaction
Effective sample volume per reaction: 20.00 uL

   P(detect)   lambda (copies/rxn)   LOD (copies/mL)
-----------------------------------------------------------
       0.900               2.303                115.1
       0.950               2.996                149.8
       0.990               4.605                230.3
       0.999               6.908                345.4
```

A 150 copies/mL LOD_95 is comfortably below the SRB threshold of concern for
oilfield souring (~10³ cells/mL) — the assay design is sufficient.

**Improving LOD without changing the primer set.** The same command lets you
model protocol changes:

| Change | Effect |
|---|---|
| Double sample volume (1 mL → 2 mL) | halves LOD in copies/mL |
| Improve extraction efficiency (40% → 70%) | LOD × 0.57 |
| Add more eluate to reaction (5 uL → 10 uL) | halves LOD, but check inhibitors |
| Concentrate eluate (100 uL → 50 uL) | halves LOD |

**Clinical point-of-care example (respiratory LAMP panel).**
Nasal swab in 400 uL VTM, 50% RNA extraction, 50 uL eluate, 5 uL to reaction:

```bash
lamp-forge lod \
  --sample-volume 400 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 400 × 0.5 × (5/50) = 20 uL per reaction → LOD_95 ≈ 150
copies/mL of VTM, equivalent to ~60 copies/swab — within the range of
published colorimetric LAMP clinical performance for SARS-CoV-2.

**Script usage.** The underlying model is importable for use in a notebook or
downstream tooling:

```python
from lamp_forge.lod import ExtractionParams, lod_table

params = ExtractionParams(
    sample_volume_ul=1000.0,
    extraction_efficiency=0.40,
    eluate_volume_ul=100.0,
    reaction_input_ul=5.0,
)
for estimate in lod_table(params):
    print(
        f"LOD_{estimate.detection_probability*100:.0f}: "
        f"{estimate.lod_genome_eq_per_ml:.1f} GE/mL"
    )
```

---

## Recipe 10 — Methanogens via *mcrA* (oil & gas souring panel, complete)

**Goal.** Detect methanogens broadly in produced water and pipeline biofilm —
the second functional guild in the oilfield souring panel described in Recipe 8.
Methanogens drive reservoir souring through a different pathway than SRB
(methanogenesis rather than sulfidogenesis) and are an independent corrosion
and gas-quality concern.

**Why it's interesting.** Methanogens are a deeply polyphyletic archaeal group
spanning Methanobacteriales, Methanomicrobiales, Methanosarcinales, and
Methanopyrales.  The universal methanogen marker is the **functional gene
`mcrA`** (methyl-coenzyme M reductase alpha subunit), which is present in
all methanogens and absent from non-methanogens (Luton et al. 2002).  Paired
with an SRB `dsrB` assay (Recipe 6) in a multiplexed panel (Recipe 8), it gives
a complete picture of the sulfide/methane risk balance in a produced-water
sample.

**Target.** `mcrA` — pull a taxonomically broad set anchored on the dominant
oil-reservoir methanogens (*Methanobacterium*, *Methanothermobacter*,
*Methanococcus*, *Methanosaeta*) and the aceticlastic lineage
(*Methanosarcina*):

```yaml
target:
  name: methanogen_mcrA
  taxon_id: 2157         # Archaea — anchor; widen with accessions below
  accessions:
    # Curate mcrA CDS across genera for oilfield coverage.
    # Start from NCBI Protein-guided nucleotide pulls for:
    # Methanobacterium thermoautotrophicum, Methanosarcina mazei,
    # Methanosaeta concilii, Methanococcus jannaschii, Methanobrevibacter.
    []
  gene: mcrA
  max_sequences: 35       # span all five major orders
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.35   # mcrA is functionally conserved but sequence-diverse
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0              # Archaea span a wide GC range
  gc_max: 65.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/methanogen_mcrA
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.** Non-methanogenic archaea and common oilfield bacteria:

| File | Source | Why |
|---|---|---|
| `sulfolobus_acidocaldarius.fasta` | NC_007181.1 | Thermoacidophilic crenarchaeon — archaea without mcrA |
| `archaeoglobus_fulgidus.fasta` | NC_000917.1 | Sulfate-reducing archaeon; shares oilfield niche with methanogens |
| `desulfovibrio_vulgaris.fasta` | NC_002937.3 | Dominant SRB in oilfields — must not co-flag in mcrA assay |
| `pseudomonas_aeruginosa.fasta` | NC_002516.2 | Ubiquitous oilfield bacterium |
| `crude_oil_microbiome_mixed.fasta` | your site metagenome | Real produced-water background |

**Expected output.** `mcrA` is a protein-coding gene under strong catalytic
constraint; conserved windows map to the alpha-helical ARS/methane-binding
domain.  Top sets should have 0 flagged hits to *Archaeoglobus* (a sulfate
reducer, no mcrA) and all bacterial off-targets.

**Completing the souring panel.** Use Recipe 8 to combine this assay with
the SRB `dsrB` set (Recipe 6) and a universal 16S control:

```bash
lamp-forge panel \
  --set SRB=results/srb_dsrB/primer_sets.json \
  --set MCR=results/methanogen_mcrA/primer_sets.json \
  --set 16S=results/general_16s/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/souring_panel
```

**Field note.** For oilfield deployment, co-quantification of SRB (Recipe 6)
and methanogens (this recipe) in a single tube provides a risk index:
high SRB + low methanogens → sulfidogenesis risk; high methanogens → gas quality
and acetate-corrosion concern.  A portable panel that resolves both guilds in
under an hour justifies the biocide and pigging decisions made on-site.

---

## Recipe 11 — PRRSV via RT-LAMP (on-farm biosecurity, RNA virus)

**Goal.** Detect porcine reproductive and respiratory syndrome virus (PRRSV)
from nasal swab, oral fluid, or serum on-farm — fast enough to inform a
quarantine or movement-restriction decision before animals are transported.

**Why it's interesting.** PRRSV is the most economically significant swine
pathogen in North America and Europe, costing the US industry alone an estimated
USD 664 million per year (Holtkamp et al. 2013).  It is an **RNA virus**
(*Arterivirus*, positive-sense ssRNA), so detection requires **RT-LAMP** — a
reverse-transcription step before the isothermal amplification.  One-step
RT-LAMP (Bst 2.0 WarmStart + NEB RTx or WarmStart RTx) runs at **63-65 degC**
in a single tube, with a positive read-out in 30 minutes, making it ideal for
the BioVind-style portable platform.

PRRSV exists as two major genotypes: **Type 1** (European, Lelystad lineage)
and **Type 2** (North American; the dominant genotype in the US and Asia).
They share ~60% nucleotide identity.  A pan-PRRSV assay must capture both, so
the primer design needs a diverse, genotype-spanning input set.

**Target.** `ORF7` (nucleocapsid protein N) — the most conserved ORF across
PRRSV genotypes and the target of most published PRRSV RT-PCR / RT-LAMP assays:

```yaml
target:
  name: prrsv_ORF7
  taxon_id: 28344   # Porcine reproductive and respiratory syndrome virus
  gene: ORF7        # nucleocapsid protein N; alt annotation: "N protein" / "N"
  max_sequences: 30 # span Type 1 + Type 2 + high-diversity strains
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 25
  entropy_threshold: 0.35   # RNA viruses are more variable — loosen slightly
  min_region_length: 200

primer:
  tm_min: 63.0              # RT-LAMP optimal: >= 63 degC for one-step RT
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0              # PRRSV is GC-moderate (~50-55%)
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/prrsv_ORF7
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config note.** `primer.tm_min: 63.0` — not the default 60.0.  One-step
RT-LAMP runs at 63-65 degC; primers with Tm below 63 degC may reduce
reverse-transcriptase co-activity.  See `lamp-forge rt-check` below.

**Off-target panel.** Porcine viruses that co-circulate and the host background:

| File | Source | Why |
|---|---|---|
| `porcine_circovirus_2.fasta` | NC_005148.1 | Ubiquitous swine DNA virus — no RT needed; must not co-flag |
| `asfv_p72.fasta` | NC_001659.1 region | Large DNA virus; the clinical confounder for respiratory cases is CSF, not ASFV, but include for completeness |
| `csfv.fasta` | classical swine fever virus *(representative)* | RNA pestivirus; most likely RNA co-detection false positive |
| `porcine_parvovirus.fasta` | NC_001718.1 | Common swine DNA virus |
| `sus_scrofa_fragment.fasta` | subset of *Sus scrofa* genome | Host RNA/DNA background from nasal swab |

**Expected output.** ORF7 is ~375 nt, so expect 1-3 conserved windows.
Top-scoring sets should have zero flagged hits across the panel; PRRSV is
genetically distant from all listed off-targets.  If you add Type 1 and Type 2
genotype sequences, a good pan-PRRSV window will appear in the N-terminal
half of ORF7 (consistent with published primer positions).

**Verify RT-LAMP readiness.** After the run, confirm all primer sets meet the
one-step RT-LAMP Tm floor before ordering:

```bash
lamp-forge rt-check \
  --input results/prrsv_ORF7/primer_sets.json \
  --na-type rna \
  --out-csv results/prrsv_ORF7/rt_check.csv
```

Example output:

```
RT-LAMP compatibility check -- target: RNA
Parameters: core primers 63.0-65.0 degC, loop primers >= 60.0 degC

Set ID                      | N | In-range | Core-low | Status
----------------------------+---+----------+----------+--------
region_01_set_001           | 6 |        6 |        0 | OK
region_01_set_002           | 6 |        5 |        1 | NOT OPTIMIZED
...

Summary: 7 of 10 set(s) RT-LAMP optimized.
Tip: tighten primer.tm_min to >= 63.0 in your config and re-run for suboptimal sets.
```

Sets marked **NOT OPTIMIZED** have at least one core primer (F3/B3/FIP/BIP)
with Tm below 63.0 degC.  Re-run with `primer.tm_min: 63.0` (already set in
the config above) and `primer.tm_max: 65.0` to push all primers into the
RT-LAMP-compatible window.

**Wet-lab note.** Oral fluid is the preferred non-invasive sample type for
PRRSV surveillance in group-housed pigs (rope sampling).  The fluid has
moderate RT-PCR inhibitors; validate with a spiked oral-fluid extraction
(e.g. MagMAX 96 Viral RNA kit) and confirm LOD at a clinically relevant
copy number (typically 10^2 - 10^3 GE/mL of oral fluid).  Use `lamp-forge lod`
to estimate the achievable LOD given your extraction protocol before ordering.

---

## Appendix — choosing your parameters for a new target

When you start from a target not covered above, set parameters in this order:

1. **Off-target panel.** Closest relatives + sample-matrix contaminants.
2. **Identity / coverage thresholds.** 80% × 80% is fine for most bacteria;
   tighten to 85% × 85% for highly conserved targets (viruses, ribosomal).
3. **Conservation threshold.** Start at 0.25 bits and adjust based on how
   many regions you get. <2 regions → loosen; >50 → tighten.
4. **Primer Tm window.** Default 60-65 °C is right for almost everyone.
5. **GC range.** Adjust if your target is GC-extreme (mycobacteria > 65%,
   *Plasmodium* < 30%).

When in doubt, run twice with different thresholds and compare HTML reports.
The pipeline is fast enough that A/B sweeps are practical.

---

## Recipe 12 — Avian Influenza A via RT-LAMP (on-farm / wild-bird biosecurity)

**Goal.** Detect highly pathogenic avian influenza A (HPAI, H5N1/H5N2/H7N9
and related subtypes) from cloacal or oropharyngeal swabs on-farm, without
a laboratory cold chain, in under 60 minutes.

**Why it's interesting.** HPAI H5N1 clade 2.3.4.4b caused catastrophic losses
in poultry flocks across North America and Europe from 2021 onward, with
ongoing spread into dairy cattle herds.  Rapid on-site detection using a
BioVind-style portable platform is exactly the workflow that allows
quarantine and culling decisions to be made before the virus spreads to
adjacent houses.  This recipe targets the influenza A **M-gene** (matrix
protein M1/M2) — the most conserved genomic segment across all influenza A
subtypes and the universal target of WHO-recommended influenza A molecular
assays.

Influenza A is a **negative-sense ssRNA** virus, so detection requires
**RT-LAMP**.  Design primers against the cDNA sequence; the pipeline handles
this automatically because sequences are deposited as cDNA in NCBI.  After
the run, validate Tm suitability with `lamp-forge rt-check` before ordering.

**Target.** Segment 7 (M-gene) — universal influenza A marker; does not
discriminate H5 from H3 or H1, but that is desirable for a first-pass
pan-influenza-A screen.  For subtype resolution, run a second assay targeting
the HA gene with H5- or H7-specific primers (separate design run).

```yaml
target:
  name: influenza_a_M_gene
  taxon_id: 11520   # Influenza A virus
  gene: M           # matrix protein; alt annotation: "M segment", "segment 7"
  max_sequences: 40 # span H5N1, H5N2, H7N9, H3N2, H1N1 for broad conservation
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight — influenza B is the closest non-target
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.35        # RNA viruses vary more than bacteria
  min_region_length: 200

primer:
  tm_min: 63.0                   # RT-LAMP one-step: co-activity with RTx enzyme
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/influenza_a_M
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config note.** `primer.tm_min: 63.0` — one-step RT-LAMP at 63-65 degC;
same rationale as Recipe 11 (PRRSV).  Run `lamp-forge rt-check` after to
confirm all sets fall in the optimised Tm window.

**Off-target panel.** Other respiratory pathogens and host background from
oropharyngeal swabs:

| File | Source | Why |
|---|---|---|
| `influenza_b_M.fasta` | NC_002204.1 | Closest non-target: IBV segment 7; ~40% aa identity to IAV M1 |
| `newcastle_disease_virus.fasta` | NC_002617.1 | Co-circulating avian paramyxovirus |
| `infectious_bronchitis_virus.fasta` | NC_001451.1 | Avian coronavirus, respiratory niche |
| `avian_metapneumovirus.fasta` | NC_007652.1 | AMPV, respiratory, co-circulates with HPAI |
| `gallus_gallus_fragment.fasta` | Subset of GalGal6 chr1 | Host DNA from cloacal swab |

**Expected output.** The IAV M-gene has a well-known conserved central
region (~positions 300-900 nt) that anchors all published RT-PCR primer
pairs.  Top-scoring LAMP-Forge sets should localise there, with zero flagged
hits against IBV, NDV, or AMPV.

**Verify RT-LAMP readiness.**

```bash
lamp-forge rt-check \
  --input results/influenza_a_M/primer_sets.json \
  --na-type rna \
  --out-csv results/influenza_a_M/rt_check.csv
```

**Multiplex with other farm-biosecurity targets.** HPAI often co-circulates
with PRRSV, ASFV, and other notifiable pathogens.  Design each target
independently, then combine:

```bash
lamp-forge panel \
  --set IAV=results/influenza_a_M/primer_sets.json \
  --set PRRSV=results/prrsv_ORF7/primer_sets.json \
  --set ASFV=results/asfv_p72/primer_sets.json \
  --top-per-target 5 \
  --out results/farm_biosecurity_panel
```

**Then generate the pooling sheet** to tell the wet-lab how to combine all
18 (3 x 6) primer stocks into a single 10x working mix:

```bash
lamp-forge pool \
  --panel results/farm_biosecurity_panel/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/farm_biosecurity_panel/pool_sheet.csv
```

The `pool_sheet.csv` tells you exactly how many microlitres of each primer
synthesis tube to pipette.  For a 3-target panel with 200 uM stocks the
water fraction is (1 - 3x44/200) * 500 = 170 uL, leaving plenty of room.

**Field note.** HPAI is a Tier 1/2 select agent in some jurisdictions.  Use
inactivated positive controls (e.g. heat-killed virus at a certified BSL-2
facility) for initial wet-lab validation.  The in silico assay design this
pipeline produces does not require BSL-3 handling.

---

## Recipe 13 — Pooling a multiplex working stock (`lamp-forge pool`)

**Goal.** After designing and screening a multi-target panel with
`lamp-forge panel`, produce the pipetting sheet that tells the wet-lab
*exactly* how to combine the individual primer synthesis tubes into a single
10x working-stock pool ready for the portable device.

**Why it's needed.** Primer synthesis vendors (IDT, Twist) ship each oligo
in its own tube at 100 uM.  For a multiplex LAMP panel with N targets you
have 6N separate tubes.  Pipetting the correct volume from each into a common
vessel — respecting the role-specific LAMP stoichiometry (FIP/BIP at 8x the
molar concentration of outer primers) — is error-prone by hand and grows
quadratically in complexity with panel size.  `lamp-forge pool` calculates it
in one command.

**LAMP 10x working-stock stoichiometry** (Li et al. 2018 / NEB E1700):

| Role | Pool conc (uM) | Final in 1x reaction (uM) |
|---|---|---|
| FIP, BIP | 16 | 1.6 |
| LF, LB | 4 | 0.4 |
| F3, B3 | 2 | 0.2 |

**Stock-concentration feasibility.**  The minimum required stock concentration
equals the sum of all primer pool concentrations across all targets.  For a
6-primer set that is 44 uM per target:

| Targets (N) | Min stock (uM) | Available from IDT/Twist? |
|---|---|---|
| 1 | 44 | Yes — standard 100 uM resuspension |
| 2 | 88 | Yes — standard 100 uM |
| 3 | 132 | Request 200 uM resuspension |
| 6 | 264 | Request 500 uM resuspension |
| 11 | 484 | Request 500 uM or vendor pre-mix |
| 18 | 792 | Vendor pre-mix recommended |

If `lamp-forge pool` raises a stock-concentration error, increase
`--stock-conc` to match the minimum shown in the error message.

**Usage.**

```bash
# After lamp-forge panel has produced panel.json:

# 2-target panel at standard 100 uM stocks, 500 uL pool:
lamp-forge pool \
  --panel results/souring_panel/panel.json \
  --stock-conc 100 \
  --total-volume 500 \
  --out results/souring_panel/pool_sheet.csv

# 3-target oilfield souring panel (SRB + MCR + 16S), 200 uM stocks:
lamp-forge pool \
  --panel results/souring_panel/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/souring_panel/pool_sheet.csv
```

**Reading the output.**  `pool_sheet.csv` has one row per primer plus a final
WATER row.  Pipette in any order; the last step is topping up to the total
volume with nuclease-free water.  Each target's primers are grouped together
so you can process one synthesis plate at a time.

```
target  role  primer_name  sequence  stock_conc_um  target_conc_um  vol_stock_ul
SRB     F3    SRB_F3       ACGT...   100.0          2.0             10.000
SRB     B3    SRB_B3       TGCA...   100.0          2.0             10.000
SRB     FIP   SRB_FIP      ACGT...   100.0         16.0             80.000
SRB     BIP   SRB_BIP      TGCA...   100.0         16.0             80.000
SRB     LF    SRB_LF       GGCC...   100.0          4.0             20.000
SRB     LB    SRB_LB       CCGG...   100.0          4.0             20.000
MCR     F3    MCR_F3       ...
...
WATER   ...   Nuclease-free water              ...             280.000
```

**Script usage.** The pool calculator is also importable:

```python
from lamp_forge.pool import PoolingParams, build_pool_plan, write_pool_csv
import json
from pathlib import Path

panel = json.loads(Path("results/souring_panel/panel.json").read_text())
params = PoolingParams(stock_conc_um=200.0, total_pool_volume_ul=500.0)
plan = build_pool_plan(panel["selection"], params)
write_pool_csv(plan, Path("results/souring_panel/pool_sheet.csv"))
print(f"{plan.n_primers} primers, {plan.water_volume_ul:.1f} uL water")
```

---

## Recipe 14 — Foot-and-mouth disease virus (FMDV) via RT-LAMP (on-farm / checkpoint biosecurity)

**Goal.** Detect FMDV broadly across all seven serotypes (O, A, C, SAT1, SAT2, SAT3, Asia1)
from epithelial vesicular fluid, blood, nasal swab, or oral fluid — fast enough to justify
immediate quarantine before suspect animals move.

**Why it's interesting.** FMDV is the most economically damaging livestock disease on the
planet: the 2001 UK outbreak cost an estimated GBP 8 billion and triggered the slaughter of
>10 million animals.  It is a **positive-sense ssRNA** virus (*Picornaviridae*, *Aphthovirus*)
— unlike ASFV, which is a large DNA virus (Recipe 7) — so detection requires **RT-LAMP**.
One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx) at 63-65 degC delivers a yes/no result in
under 30 minutes, compatible with a BioVind-style portable device at the farm gate or a
border checkpoint.

FMDV has **seven antigenically distinct serotypes** (O, A, C, SAT1, SAT2, SAT3, Asia1)
sharing as little as 30% VP1 nucleotide identity across serotypes.  A pan-FMDV assay must
target a conserved region.  The **3D RNA-dependent RNA polymerase** gene (~1800 nt) is the
most conserved coding region across all serotypes and is the target of published RT-PCR and
RT-LAMP FMDV assays (de Vries et al. 2010; Soltan et al. 2018).

**Target.** 3D gene (RNA polymerase), pulling sequences spanning all seven serotypes.
Use the ready-made config at `config/fmdv_3dpol.yaml`, or paste the block below:

```yaml
target:
  name: fmdv_3Dpol
  taxon_id: 12110          # Foot-and-mouth disease virus
  gene: 3D                 # RNA polymerase; alt annotation: "3Dpol" / "RdRp"
  max_sequences: 40        # span O, A, C, SAT1, SAT2, SAT3, Asia1
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- SVD is a picornavirus; cross-reactivity risk
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.30        # RNA-virus drift: slightly looser than bacteria
  min_region_length: 200

primer:
  tm_min: 63.0                   # RT-LAMP one-step: co-activity floor for NEB RTx / AMV-RT
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/fmdv_3Dpol
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config note.** `primer.tm_min: 63.0` — same one-step RT-LAMP rationale as Recipes 11
and 12 (PRRSV and AIV).  The `conservation.entropy_threshold: 0.30` is slightly looser than
for bacterial targets — the 3D gene is under strong catalytic constraint but still subject to
RNA-virus synonymous drift.

**Off-target panel.** The critical differential is **vesicular disease** — four distinct
diseases produce indistinguishable clinical signs (vesicles on feet, tongue, snout):

| File | Source | Why |
|---|---|---|
| `svd.fasta` | swine vesicular disease virus *(representative)* | Vesicular lesions identical to FMDV in pigs; closest picornavirus in panel |
| `vs_virus.fasta` | vesicular stomatitis virus *(representative)* | Same clinical picture; rhabdovirus (cattle, horses, pigs) |
| `csfv.fasta` | classical swine fever virus *(representative)* | Co-circulates in pigs; RNA pestivirus |
| `bvdv.fasta` | bovine viral diarrhea virus *(representative)* | Common bovine RNA pestivirus |
| `bos_taurus_fragment.fasta` | Subset of ARS-UCD2.0 chr1 | Host DNA from cattle blood or epithelial swab |
| `sus_scrofa_fragment.fasta` | Subset of *Sus scrofa* genome | Host DNA from pig samples |

**Expected output.** The 3D gene has a well-documented conserved central block (~780-1200 nt
of the 3D ORF relative to the O/Kaufbeuren reference).  Top-scoring sets should localise there,
with zero flagged off-target hits — FMDV is genetically distant from CSFV and BVD (pestiviruses),
VSV (rhabdovirus), and SVD (Enterovirus E).

**Verify RT-LAMP readiness.**

```bash
lamp-forge rt-check \
  --input results/fmdv_3Dpol/primer_sets.json \
  --na-type rna \
  --out-csv results/fmdv_3Dpol/rt_check.csv
```

Sets marked **NOT OPTIMIZED** have at least one core primer below 63.0 degC.  Because
`primer.tm_min: 63.0` is already set in this config, re-running is usually not necessary
— but verify before ordering.

**Multiplex with other farm-biosecurity targets.** FMDV, ASFV, PRRSV, and AIV can co-occur
in high-density livestock areas.  Design each target independently, then combine:

```bash
lamp-forge panel \
  --set FMDV=results/fmdv_3Dpol/primer_sets.json \
  --set ASFV=results/asfv_p72/primer_sets.json \
  --set PRRSV=results/prrsv_ORF7/primer_sets.json \
  --set IAV=results/influenza_a_M/primer_sets.json \
  --top-per-target 5 \
  --out results/farm_biosecurity_4plex
```

Then pool and export:

```bash
lamp-forge pool \
  --panel results/farm_biosecurity_4plex/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/farm_biosecurity_4plex/pool_sheet.csv

lamp-forge export \
  --input results/fmdv_3Dpol/primer_sets.json \
  --format idt \
  --target-label FMDV_3D \
  --out orders/fmdv_idt_order.csv
```

**Estimate LOD before ordering** (nasal swab in 2 mL VTM, 50% RNA extraction, 100 uL eluate,
5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 2000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5
```

Effective sample = 2000 x 0.5 x (5/100) = 50 uL per reaction -> LOD_95 approx 60 copies/mL
of VTM.  FMDV viral loads in epithelial vesicular fluid during the acute febrile phase can
reach 10^8 copies/mL, so a 60 copies/mL LOD provides orders-of-magnitude headroom.

**Field note.** FMDV is a WOAH-listed Tier 1 notifiable disease.  A positive result should
trigger immediate notification to the national veterinary authority and confirmatory testing
at a WOAH Reference Laboratory (e.g. Pirbright Institute, FAO World Reference Laboratory for
FMD).  Position this assay as a **triage/screening** tool, not the confirmatory test.  Wet-lab
validation must be performed at a certified facility using inactivated positive controls or
synthetic RNA standards — the in silico design produced by this pipeline does not require
BSL-2+ handling.

**References.**

- de Vries AAF et al. (2010) Real-time RT-LAMP for rapid, sensitive detection of FMDV
  serotypes O, A and Asia 1. *J Virol Methods* 163:303-311.
  doi:10.1016/j.jviromet.2009.09.029
- Soltan MA et al. (2018) Loop-mediated isothermal amplification (LAMP) for rapid detection
  of FMDV serotype O. *J Virol Methods* 251:6-10. doi:10.1016/j.jviromet.2017.10.003
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA. *Nucleic Acids Res*
  28(12):e63. doi:10.1093/nar/28.12.e63
