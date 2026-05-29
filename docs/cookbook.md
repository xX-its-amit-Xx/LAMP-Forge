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
