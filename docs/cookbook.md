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

---

## Recipe 15 — Iron-Reducing Bacteria (IRB) via *omcA* (oilfield MIC monitoring)

**Goal.** Detect Shewanella-type iron-reducing bacteria (IRB) in oilfield
produced water or pipeline biofilm using a LAMP assay targeting the outer
membrane decaheme cytochrome A gene (*omcA*) — adding a third channel to
a MIC/souring panel alongside SRB (Recipe 6 / dsrB) and methanogens
(Recipe 10 / mcrA).

**Why it's interesting.** IRB catalyse the reductive dissolution of ferric
iron (Fe3+) protective oxide scales, directly exposing bare steel to
corrosive environments.  In the MtrCAB-OmcA extracellular electron transfer
(EET) pathway, *omcA* encodes the terminal outer-membrane decaheme cytochrome
that hands electrons to solid-phase Fe3+ (or to a steel surface acting as the
final electron acceptor).  A positive *omcA* signal in produced water or
pigging solids therefore flags active iron-cycling, a prerequisite for
under-deposit pitting corrosion — the dominant failure mode in oilfield
pipeline MIC.

Combined with an SRB (dsrB) assay, IRB and SRB together drive a corrosive
syntrophic cycle: IRB mobilise iron from passive oxide scales; SRB reduce
sulfate to H2S, which re-precipitates as highly corrosive iron sulfide
(FeS / mackinawite).  Monitoring all three functional guilds (SRB +
methanogens + IRB) from a single sample enables proportional risk scoring
rather than a binary positive/negative call.

*omcA* is the best single-gene marker for Shewanella-type IRB because:

- It is essential for iron reduction and expressed only under anaerobic /
  micro-aerophilic conditions relevant to biofilm niches.
- It is phylogenetically restricted to the *Shewanella* clade — Geobacter
  uses distinct outer-membrane cytochromes (omcS, omcB, omcE).
- Its 10-heme CXXCH repeat structure yields conserved blocks at the
  amino-acid level while diverging at the nucleotide level across species —
  exactly the profile LAMP-Forge's entropy filter is designed to exploit.

**Target.** *Shewanella* genus (NCBI taxon 22), *omcA* gene.  Use the
ready-made config at `config/irb_omcA.yaml`, or paste the block below:

```yaml
target:
  name: irb_omcA_oilfield
  taxon_id: 22               # Shewanella genus
  accessions: []             # extend with S. putrefaciens, S. baltica, S. loihica
  gene: omcA
  max_sequences: 30
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.30
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/irb_omcA
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `conservation.entropy_threshold: 0.30` — omcA sits between tight
  housekeeping-gene designs (< 0.20 for rpoB, invA) and the looser
  polyphyletic functional-gene threshold used for dsrB / mcrA (> 0.35).
  The heme-binding CXXCH motifs are virtually invariant at the amino-acid
  level, but synonymous divergence at the DNA level is substantial across
  *Shewanella* species, so 0.30 is the empirical sweet spot.
- `gc_min/gc_max: 40-60%` — *Shewanella* GC content is 45-52% depending on
  species; the slightly wider window reduces the risk that primer design fails
  to find a compliant primer in every conserved window.
- `target.accessions` — add CDS records for *S. putrefaciens* CN-32,
  *S. baltica* OS678, *S. loihica* PV-4, *S. frigidimarina* NCIMB 400, and
  *S. amazonensis* SB2B to cover the phylogenetic breadth of oilfield
  *Shewanella* isolates.

**Off-target panel.** The key differentials are non-IRB bacteria that carry
their own cytochrome c genes:

| File | Source | Why |
|---|---|---|
| `pseudomonas_aeruginosa.fasta` | PA01 RefSeq | Ubiquitous oilfield biofilm bacterium; has respiratory cytochromes but not 10-heme omcA |
| `desulfovibrio_vulgaris.fasta` | Hildenborough RefSeq | Dominant SRB in oilfields; must not cross-flag in the IRB channel |
| `geobacter_metallireducens.fasta` | GS-15 RefSeq | Different IRB with phylogenetically distant outer-membrane cytochromes (omcS); cross-flag flags unexpectedly broad assay scope |
| `escherichia_coli.fasta` | K-12 MG1655 RefSeq | Common contamination; has monoheme cytochromes but not 10-heme outer-membrane cytochromes |
| `marinobacter_hydrocarbonoclasticus.fasta` | SP17 RefSeq | Common hydrocarbon-oxidising oilfield bacterium |
| `crude_oil_microbiome_mixed.fasta` | Site metagenome | Produced-water environmental background |

**Expected output.** The omcA gene (~1800 nt CDS) typically contains two or
three conserved windows usable for LAMP.  Top-scoring sets should localise to
the decaheme-core region with zero flagged off-target hits against the panel
above.

**Complete the three-channel oilfield MIC panel.** After designing SRB,
methanogen, and IRB assays independently:

```bash
lamp-forge panel \
  --set SRB=results/srb_dsrB/primer_sets.json \
  --set MCR=results/methanogen_mcrA/primer_sets.json \
  --set IRB=results/irb_omcA/primer_sets.json \
  --top-per-target 5 \
  --out results/mic_3plex_panel
```

Generate the pooling sheet for 200 uM synthesis stocks (3 targets = 132 uM
minimum required; request 200 uM resuspension from your vendor):

```bash
lamp-forge pool \
  --panel results/mic_3plex_panel/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/mic_3plex_panel/pool_sheet.csv
```

Export the IRB primers to IDT for ordering:

```bash
lamp-forge export \
  --input results/irb_omcA/primer_sets.json \
  --format idt \
  --target-label IRB_omcA \
  --out orders/irb_omcA_idt_order.csv
```

Estimate LOD for a 1 mL produced-water sample processed through a standard
DNA extraction kit (50% efficiency, 100 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5
```

Effective sample = 1000 x 0.50 x (5/100) = 25 uL per reaction.
Oilfield produced-water IRB loads at MIC-active sites typically exceed
10^4 cells/mL, so the assay LOD provides orders-of-magnitude headroom.

**Interpreting the three-channel panel result.**

| SRB (dsrB) | Methanogens (mcrA) | IRB (omcA) | Risk interpretation |
|---|---|---|---|
| + | + | + | High MIC + souring risk; recommend biocide treatment review |
| + | - | + | Active iron-cycling + sulfide; under-deposit corrosion risk, potential FeS scale |
| + | + | - | Sulfidogenesis / souring risk; lower structural corrosion signal |
| - | - | + | Iron-cycling active; corrosion risk without H2S; monitor for pitting |
| - | - | - | Low microbial activity; re-test after any process upset |

**Field note.** *Shewanella* spp. are facultative anaerobes; they can
survive oxygen exposure between sampling and analysis.  Use anaerobic
collection vessels (e.g. Hungate tubes with N2 headspace) or process samples
within 4 hours of collection to avoid shifts in community composition.
Filter-capture onto a 0.2 uM membrane and store at -20 degC if immediate DNA
extraction is not possible.

**References.**

- Myers CR & Myers JM (1992). Localization of cytochromes to the outer
  membrane of anaerobically grown *Shewanella putrefaciens* MR-1.
  *J Bacteriol* 174:3429-3438.
- Lies DP et al. (2005). *Shewanella oneidensis* MR-1 uses overlapping
  pathways for iron reduction at a distance and by direct contact under
  conditions relevant for biofilms. *Appl Environ Microbiol* 71:4414-4426.
  doi:10.1128/AEM.71.8.4414-4426.2005
- Ross DE et al. (2011). Comparative genomics of the electron transport
  chains of the *Shewanella*, using the MtrCAB pathway as a model.
  *PLoS Comput Biol* 7:e1002189. doi:10.1371/journal.pcbi.1002189
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 16 — Acid-Producing Bacteria (APB) via *fthfs* (oilfield MIC monitoring, 4th channel)

**Goal.** Detect homoacetogens — the dominant acid-producing bacteria in
oilfield produced water — using a LAMP assay targeting the *fthfs* gene
(formyltetrahydrofolate synthetase), the committed step of the Wood-Ljungdahl
carbon-fixation pathway that converts CO2 to acetic acid.

**Why it's interesting.** Acid-producing bacteria are the neglected fourth
guild in oilfield MIC.  The three canonical MIC guilds monitored by
BioVind's oil & gas panel — SRB (dsrB, Recipe 6), methanogens (mcrA,
Recipe 10), and IRB (omcA, Recipe 15) — account for H2S production, CH4
souring, and iron dissolution, respectively.  APB via *fthfs* adds the
missing piece: organic acid production.

Homoacetogens create two compounding hazards.

1. **Direct pH attack.** Acetic acid (and traces of formate, propionate)
   depress produced water pH from ~6.5–7.0 toward 4–5.  At pH < 5 the
   protective FeCO3 passivation layer on carbon-steel dissolves, exposing
   bare metal to the corrosive environment.
2. **SRB amplification loop.** Acetate and H2 — the primary products of
   homoacetogenesis — are the preferred electron donors for SRB.  An
   APB-positive sample with a negative SRB call should trigger heightened
   monitoring because the acetate source, once present, can fuel a rapid
   SRB bloom in response to any sulfate injection or process upset.

*fthfs* (formyltetrahydrofolate synthetase; EC 6.3.4.3) is the established
molecular marker for homoacetogens in environmental samples (Leaphart &
Lovell 2001; Hunger et al. 2011).  Its ATP-grasp and GHMP-kinase catalytic
domains are highly conserved within homoacetogens yet phylogenetically
restricted — most non-acetogenic Firmicutes carry *pta*/*ack* genes for
acetate kinetics but lack the full Wood-Ljungdahl *fthfs*, enabling
clean genus-level differentiation.  The ~1.7 kb CDS contains at least two
LAMP-accessible conserved windows.

**Target sequences.** Anchor to Moorella genus (NCBI taxon 44417), the
thermophilic homoacetogens dominant in high-temperature oilfield reservoirs
(50–80 °C).  Supplement with Acetobacterium woodii accessions for mesophilic
surface-facility coverage.

**Config.**  Use the ready-made config at `config/apb_fthfs.yaml`, or paste
the block below:

```yaml
target:
  name: apb_fthfs_oilfield
  taxon_id: 44417            # Moorella genus -- thermophilic homoacetogens
  accessions:
    # Supplement with Acetobacterium woodii DSM 1030 and Sporomusa sphaeroides
    # CDS records for mesophilic surface-facility coverage.
    []
  gene: fthfs
  max_sequences: 30
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.30        # fthfs: mid-range, like dsrB / mcrA

primer:
  tm_min: 60.0                   # DNA target -- standard LAMP floor (not RT-LAMP)
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0                   # Acetobacterium ~44% GC
  gc_max: 65.0                   # Moorella ~55% GC; wide window for diversity
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/apb_fthfs
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Off-target panel.**  Key differentials for specificity in an oilfield
produced-water sample:

| File | Source | Why |
|---|---|---|
| `desulfovibrio_vulgaris.fasta` | Hildenborough RefSeq | Dominant SRB — carries THF-pathway genes but not Wood-Ljungdahl *fthfs*; must not co-flag in APB channel |
| `methanobacterium_formicicum.fasta` | RefSeq | Formate-consuming methanogen; formate/THF chemistry overlaps with homoacetogenesis — key specificity test |
| `clostridium_acetobutylicum.fasta` | ATCC 824 RefSeq | Fermentative Clostridium with *pta*/*ack* but no *fthfs* — representative of non-acetogens that classical APB culture tests would also count |
| `shewanella_oneidensis.fasta` | MR-1 RefSeq | IRB; must not co-flag in APB channel |
| `pseudomonas_aeruginosa.fasta` | PA01 RefSeq | Ubiquitous oilfield biofilm bacterium; no Wood-Ljungdahl pathway |
| `escherichia_coli.fasta` | K-12 MG1655 RefSeq | Common contamination; no *fthfs* |

**Key config notes.**

- `conservation.entropy_threshold: 0.30` — *fthfs* sits between tight
  housekeeping genes (rpoB < 0.20) and the loosest polyphyletic markers
  (dsrB/mcrA > 0.30).  The ATP-grasp and GHMP-kinase active-site residues
  are near-invariant; synonymous positions diverge at the inter-genus level,
  so 0.30 bits is empirically appropriate.
- `gc_min / gc_max: 40–65%` — Moorella is ~55% GC, Acetobacterium ~44%.
  The wider window (25 GC units) ensures primer candidates are found when
  accessions from both genera are present in the sequence set.
- `tm_min: 60.0` — *fthfs* is a chromosomal DNA target; do **not** raise
  this to the RT-LAMP 63 °C floor used for PRRSV, FMDV, and AIV.
- `target.accessions` — add fthfs CDS records for *A. woodii* DSM 1030,
  *Sporomusa sphaeroides* DSM 2875, and *Thermoanaerobacter kivui* DSM 2030
  to capture both thermophilic and mesophilic oilfield niches.

**Complete the four-channel oilfield MIC panel.**  After designing SRB,
methanogen, IRB, and APB assays independently, run `lamp-forge panel` to
screen the combined set for cross-assay primer compatibility:

```bash
lamp-forge panel \
  --set SRB=results/srb_dsrB/primer_sets.json \
  --set MCR=results/methanogen_mcrA/primer_sets.json \
  --set IRB=results/irb_omcA/primer_sets.json \
  --set APB=results/apb_fthfs/primer_sets.json \
  --top-per-target 5 \
  --out results/mic_4plex_panel
```

Generate the pooling sheet (four targets = 176 uM minimum; request 200 uM
resuspension from IDT or Twist):

```bash
lamp-forge pool \
  --panel results/mic_4plex_panel/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/mic_4plex_panel/pool_sheet.csv
```

Export APB primers for IDT ordering:

```bash
lamp-forge export \
  --input results/apb_fthfs/primer_sets.json \
  --format idt \
  --target-label APB_fthfs \
  --out orders/apb_fthfs_idt_order.csv
```

**Interpreting the four-channel MIC panel result.**

| SRB | MCR | IRB | APB | Risk interpretation |
|---|---|---|---|---|
| + | + | + | + | Maximum MIC + souring risk; all four corrosion guilds active; immediate biocide and nitrate-injection review |
| + | - | + | + | APB + SRB syntrophic loop active; iron cycling; high structural pitting + acid corrosion risk |
| - | - | - | + | Acetogenic baseline; SRB not yet established; monitor — a sulfate injection or process upset can trigger rapid SRB bloom |
| + | + | - | - | Sulfidogenesis + gas souring without iron/acid amplification; standard souring protocol |
| - | - | - | - | Low microbial activity; re-test after process upset or water injection |

**Field note.**  Homoacetogens are strict anaerobes.  Collect produced-water
samples into anaerobic Hungate tubes (N2 headspace) or serum vials sealed
immediately after collection.  Freeze-preserve at −20 °C on a 0.2 µm
membrane filter within four hours if DNA extraction is not possible on-site.
Cell concentrations in MIC-active produced water typically exceed 10³
cells/mL — well above the LOD achievable with standard 200 µL sample
extraction (use `lamp-forge lod` to confirm with your site-specific extraction
parameters).

**References.**

- Hunger S, Schmidt O, Hilgarth M et al. (2011). Competing formate- and
  carbon dioxide-utilizing prokaryotes in an anoxic methane-emitting fen soil.
  *Environ Microbiol* 13:2228–2238. doi:10.1111/j.1462-2920.2011.02491.x
- Leaphart AB & Lovell CR (2001). Recovery and analysis of *fthfs* gene
  sequences as markers for acetogenic bacteria in soil and sediment microbial
  communities. *Appl Environ Microbiol* 67:2720–2728.
  doi:10.1128/AEM.67.6.2720-2728.2001
- Pierce E, Xie G, Barabote RD et al. (2008). The complete genome sequence
  of *Moorella thermoacetica*. *Environ Microbiol* 10:2550–2573.
  doi:10.1111/j.1462-2920.2008.01679.x
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 17 — Nitrate-Reducing Bacteria (NRB) via *narG* (oilfield souring control, 5th channel)

**Goal.** Detect nitrate-reducing bacteria (NRB) in oilfield produced water
or pipeline biofilm using a LAMP assay targeting *narG* (dissimilatory nitrate
reductase alpha subunit) — adding the fifth and final functional-guild channel
to a comprehensive MIC/souring panel alongside SRB (Recipe 6 / *dsrB*),
methanogens (Recipe 10 / *mcrA*), IRB (Recipe 15 / *omcA*), and APB
(Recipe 16 / *fthfs*).

**Why it's interesting.**  Nitrate injection is the primary *biological*
souring-control strategy in reservoir management.  Injected nitrate
selectively stimulates NRB, which:

1. **Outcompete SRB** for shared electron donors (H₂, acetate, lactate)
   through thermodynamic favourability (nitrate reduction is more energetically
   favourable than sulfate reduction).
2. **Produce nitrite**, a direct inhibitor of dissimilatory sulfate reductase,
   suppressing H₂S production at the enzymatic level.

The critical diagnostic question is not simply "are NRB present?" but
**"is the NRB:SRB ratio sufficient to suppress souring?"**  A rising NRB
signal alongside a declining *dsrB* signal is the expected outcome of a
successful nitrate-injection treatment.  A positive NRB signal that fails
to suppress SRB indicates the treatment dose is insufficient or SRB have
acquired nitrite resistance.

NRB are phylogenetically diverse (Paracoccus, Pseudomonas, Thiobacillus,
Thauera, Arcobacter span multiple phyla) and share only functional affinity.
A gene-level assay targeting *narG* is the correct approach, just as *dsrB*
captures phylogenetically diverse SRB.

**Target.**  *narG* — the catalytic MGD-binding subunit of respiratory
nitrate reductase (EC 1.7.99.4), the field-standard functional marker for
NRB in environmental microbiology (Braker et al. 1998; Castillo et al. 2019).
Use the ready-made config at `config/nrb_narG.yaml`, or paste the block below:

```yaml
target:
  name: nrb_narG_oilfield
  taxon_id: 265              # Paracoccus genus -- canonical NRB model clade
  accessions:
    # Supplement with narG CDS records from key oilfield NRB genera:
    #   Pseudomonas stutzeri JM300
    #   Thiobacillus denitrificans ATCC 25259
    #   Thauera aromatica K172
    #   Arcobacter nitrofigilis DSM 7299 (napA, periplasmic nitrate reductase)
    []
  gene: narG
  max_sequences: 30
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.35   # polyphyletic functional gene -- same profile as dsrB/mcrA

primer:
  tm_min: 60.0               # DNA target -- standard LAMP (not RT-LAMP)
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0
  gc_max: 70.0               # GC-rich NRB genera (Paracoccus ~67%, Pseudomonas ~63%)
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/nrb_narG
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `conservation.entropy_threshold: 0.35` — *narG* is a polyphyletic functional
  gene with the same conservation profile as *dsrB* and *mcrA*: the
  MGD-binding Cys/Met residues are near-invariant, but synonymous divergence
  between genera is high.
- `gc_min/gc_max: 40-70%` — dominant oilfield NRB genera are GC-rich
  (Paracoccus ~67%, Pseudomonas stutzeri ~63%, Thiobacillus denitrificans
  ~66%).  The 30-unit window ensures primer candidates can be found when
  sequences from multiple genera are present.
- `tm_min: 60.0` — *narG* is a chromosomal DNA target; do **not** raise to
  the RT-LAMP 63 °C floor used for RNA-virus recipes.
- Supplement `target.accessions` before running: *narG* CDS records from
  *Thiobacillus denitrificans*, *Pseudomonas stutzeri*, *Thauera aromatica*,
  and *Arcobacter* spp. extend coverage to the full oilfield NRB community.

**Off-target panel.**  The essential specificity requirement is clean
separation from SRB — the two guilds are biological antagonists and the
assay's diagnostic value depends on distinguishing them:

| File | Source | Why |
|---|---|---|
| `desulfovibrio_vulgaris.fasta` | Hildenborough RefSeq | Dominant oilfield SRB — must **not** co-flag in NRB channel |
| `methanobacterium_formicicum.fasta` | RefSeq | Methanogen; no *narG* |
| `shewanella_oneidensis.fasta` | MR-1 RefSeq | IRB; uses Fe(III) not nitrate; must not co-flag |
| `escherichia_coli.fasta` | K-12 MG1655 RefSeq | Carries *narGHJI* for anaerobic respiration but is not a denitrifier; cross-flag indicates over-broad primers |
| `pseudomonas_aeruginosa.fasta` | PA01 RefSeq | Has *narG* denitrification genes; cross-reactivity is biologically acceptable (PA01 is also an NRB) |
| `crude_oil_microbiome_mixed.fasta` | Site metagenome | Produced-water environmental background |

**Expected output.**  *narG* is a ~2.7 kb CDS; expect two or three conserved
windows.  Top sets should show zero flagged hits against *Desulfovibrio*,
methanogens, and *Shewanella*.  *Pseudomonas aeruginosa* cross-reactivity is
expected and biologically acceptable — PA01 IS an NRB.

**Complete the five-channel oilfield MIC panel.**  After designing all five
assays independently:

```bash
lamp-forge panel \
  --set SRB=results/srb_dsrB/primer_sets.json \
  --set MCR=results/methanogen_mcrA/primer_sets.json \
  --set IRB=results/irb_omcA/primer_sets.json \
  --set APB=results/apb_fthfs/primer_sets.json \
  --set NRB=results/nrb_narG/primer_sets.json \
  --top-per-target 5 \
  --out results/mic_5plex_panel
```

Generate the pooling sheet (five targets = 220 uM minimum; request 250 uM
resuspension from IDT or Twist):

```bash
lamp-forge pool \
  --panel results/mic_5plex_panel/panel.json \
  --stock-conc 250 \
  --total-volume 500 \
  --out results/mic_5plex_panel/pool_sheet.csv
```

Export NRB primers for IDT ordering:

```bash
lamp-forge export \
  --input results/nrb_narG/primer_sets.json \
  --format idt \
  --target-label NRB_narG \
  --out orders/nrb_narG_idt_order.csv
```

Estimate LOD for a 1 mL produced-water sample (50% DNA extraction, 100 uL
eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5
```

Effective sample = 1000 x 0.50 x (5/100) = 25 uL per reaction.  NRB loads
in nitrate-injected produced water during active treatment typically exceed
10⁴ cells/mL, providing orders-of-magnitude headroom above the achievable LOD.

**Interpreting the five-channel panel result.**

| SRB | MCR | IRB | APB | NRB | Risk / treatment interpretation |
|---|---|---|---|---|---|
| - | - | - | - | + | Nitrate injection working; NRB dominant; low souring risk |
| + | - | + | + | + | NRB present but not suppressing SRB; all corrosion guilds active; escalate nitrate dose or add biocide |
| + | + | + | + | - | All four MIC guilds active; no NRB; nitrate not injected or exhausted; maximum risk |
| + | - | - | - | + | SRB and NRB co-present; treatment under-dosed; SRB not yet suppressed |
| - | - | - | - | - | Low microbial activity; maintain monitoring schedule |

**Field note.**  NRB are facultative anaerobes that survive aerobic sample
handling, unlike the strict anaerobes (SRB, methanogens) also measured in
this panel.  For a complete five-guild panel, use anaerobic Hungate tubes to
preserve all community members from collection through extraction.

**References.**

- Hubert C & Voordouw G (2007). Oil field souring control by nitrate-reducing
  *Sulfurospirillum* spp. that outcompete sulfate-reducing bacteria for organic
  electron donors. *Appl Environ Microbiol* 73:2644–2652.
  doi:10.1128/AEM.02403-06
- Castillo JA, Agathos SN & de los Cobos-Vasconcelos D (2019). Functional gene
  diversity and microbial community dynamics in oilfield souring and corrosion.
  *Front Microbiol* 10:2534. doi:10.3389/fmicb.2019.02534
- Braker G, Fesefeldt A & Witzel KP (1998). Development of PCR primer systems
  for amplification of nitrite reductase genes to detect denitrifying bacteria
  in environmental samples. *Appl Environ Microbiol* 64:3769–3775.
- Voordouw G (2011). Production-related petroleum microbiology: progress and
  prospects. *Curr Opin Biotechnol* 22:401–405.
  doi:10.1016/j.copbio.2011.01.015
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 18 -- Newcastle Disease Virus (NDV) via RT-LAMP (on-farm / border-post biosecurity)

**Goal.** Detect Newcastle Disease Virus (NDV) broadly across all major Class I
and Class II genotypes from cloacal or oropharyngeal swabs on-farm or at a
checkpoint -- fast enough to inform a quarantine decision before birds move.

**Why it's interesting.** NDV (Avian orthoavulavirus 1) is one of the most
economically and epidemiologically significant avian pathogens worldwide.
Velogenic strains (Class II genotypes V, VII, XII, XXI) cause near-100%
flock mortality and are WOAH/OIE-listed as a notifiable disease.  The
current global situation -- with HPAI H5N1 clade 2.3.4.4b circulating in
poultry and wild birds -- makes rapid multi-target avian respiratory
surveillance essential.  NDV is the natural fourth target in a farm
biosecurity panel alongside AIV (Recipe 12), PRRSV (Recipe 11), and ASFV
(Recipe 7): one portable device, one swab, four yes/no answers.

NDV is a **negative-sense ssRNA paramyxovirus** (*Paramyxoviridae*,
*Avulavirinae*), so detection requires **RT-LAMP**.  One-step RT-LAMP
(Bst 2.0 WarmStart + NEB RTx or WarmStart RTx) at 63-65 degC delivers a
positive signal in under 30 minutes, ideal for the BioVind portable platform
at the farm gate or a border inspection post.

**Choice of target gene.**  The NDV genome encodes six proteins; two are
commonly used for molecular detection:

- **M gene** (matrix protein, ~1077 nt CDS): most conserved across all
  genotypes I-XXI and both Class I and Class II.  The basis of the
  OIE-endorsed real-time RT-PCR assay (Wise et al. 2004) and several
  published RT-LAMP assays.  Detects any NDV regardless of pathotype --
  the correct scope for first-pass screening.
- **F gene** (fusion protein, ~1662 nt CDS): carries the multi-basic cleavage
  site that distinguishes velogenic from lentogenic strains.  Useful for
  pathotyping after detection, but inter-genotype divergence is higher, so
  broad-inclusivity primer design is harder.

Use the **M gene** (this recipe and `config/ndv_M_gene.yaml`) for detection;
if downstream pathotype information is needed, design a follow-up F-gene assay
using `lamp-forge scaffold --target-name ndv_F_gene --vertical farm --na-type rna`.

**Target.** M gene, pulling sequences spanning Class I (genotype 1, which
includes avirulent aquatic-bird strains) and the major Class II genotypes
(I-XXI), including the virulent genotypes responsible for active outbreaks
(VII, XII, XXI in Asia/Africa/Middle East).  Use the ready-made config at
`config/ndv_M_gene.yaml`, or paste the block below:

```yaml
target:
  name: ndv_M_gene
  taxon_id: 11234          # Newcastle disease virus (Avian orthoavulavirus 1)
  gene: M                  # matrix protein; most conserved across all genotypes
  max_sequences: 35        # span Class I + Class II genotypes I-XXI
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- APMV-2/3 share paramyxovirus family
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.30        # RNA virus, M gene moderately conserved across genotypes
  min_region_length: 200

primer:
  tm_min: 63.0                   # RT-LAMP one-step: co-activity floor for NEB RTx / Bst 2.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 40.0                   # NDV M gene ~47% GC; wide window for inter-genotype diversity
  gc_max: 60.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/ndv_M_gene
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `primer.tm_min: 63.0` -- same one-step RT-LAMP rationale as Recipes 11
  (PRRSV), 12 (AIV), and 14 (FMDV).  NDV is a negative-sense ssRNA virus;
  the RT enzyme must remain co-active with Bst polymerase at 63-65 degC.
- `conservation.entropy_threshold: 0.30` -- the M gene is under structural
  constraint (matrix-layer assembly drives the threshold lower than dsrB/mcrA
  at 0.35), but RNA-virus synonymous drift across 20+ genotypes keeps it
  above the housekeeping-gene limit (rpoB < 0.20).
- `off_targets.min_identity_threshold: 0.85` -- tighter than the bacterial
  default (0.80) because avian paramyxoviruses (APMV-2, APMV-3) share
  paramyxovirus structural motifs with NDV.  A 0.85 threshold catches
  potential cross-reactivity before it reaches the lab.
- `gc_min/gc_max: 40-60%` -- the NDV M gene is ~47% GC in the OIE LaSota
  reference (Class II genotype II).  Velogenic outbreak strains (genotypes
  VII, XXI) shift GC by +/- 5 percentage points; the 20-unit window covers
  this diversity without being so wide that non-specific primers are accepted.

**Off-target panel.** Avian paramyxoviruses and co-circulating respiratory
pathogens from cloacal or oropharyngeal swabs:

| File | Source | Why |
|---|---|---|
| `apmv2.fasta` | Avian paramyxovirus type 2 (Yucaipa) | Closest APMV relative; ~40% M-gene aa identity; must not co-flag |
| `apmv3.fasta` | Avian paramyxovirus type 3 (turkey paramyxovirus) | Second-closest APMV relative |
| `infectious_bronchitis_virus.fasta` | NC_001451.1 | Avian coronavirus; common respiratory co-circulator |
| `avian_metapneumovirus.fasta` | NC_007652.1 | AMPV; co-circulates with NDV and HPAI |
| `influenza_a_h5n1.fasta` | representative H5N1 clade 2.3.4.4b | Co-circulates in poultry outbreaks; separate assay (Recipe 12) |
| `gallus_gallus_fragment.fasta` | GalGal6 chr1 subset | Host DNA from cloacal or oropharyngeal swab |

**Expected output.** The M gene has a well-conserved central block (~340-780
nt of the ~1077 nt M CDS) that aligns robustly across Class I and Class II
genotypes.  Top-scoring sets should localise to this block with zero flagged
hits against APMV-2, APMV-3, IBV, and AMPV.  Sensitivity check: if fewer
than 2 conserved windows appear, loosen `entropy_threshold` to 0.35 -- highly
virulent outbreak strains in the genotype VII/XXI clade have elevated
synonymous divergence relative to classical strains.

**Verify RT-LAMP readiness.**

```bash
lamp-forge rt-check \
  --input results/ndv_M_gene/primer_sets.json \
  --na-type rna \
  --out-csv results/ndv_M_gene/rt_check.csv
```

Sets marked **NOT OPTIMIZED** have at least one core primer (F3/B3/FIP/BIP)
below 63.0 degC.  Because `primer.tm_min: 63.0` is already set in this
config, re-running is usually not required, but confirm before ordering.

**Add NDV to the avian farm-biosecurity panel.**  Recipe 12 designs an AIV
M-gene assay.  Adding NDV creates a complete avian respiratory 2-plex (or
a 4-plex when PRRSV and ASFV are included):

```bash
# 1. Design AIV and NDV independently.
lamp-forge run --config config/influenza_a_M.yaml
lamp-forge run --config config/ndv_M_gene.yaml

# 2. Screen for cross-assay primer compatibility.
lamp-forge panel \
  --set IAV=results/influenza_a_M/primer_sets.json \
  --set NDV=results/ndv_M_gene/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/avian_respiratory_panel

# 3. Generate the pipetting sheet (2 RNA targets = 88 uM minimum; 100 uM stocks suffice).
lamp-forge pool \
  --panel results/avian_respiratory_panel/panel.json \
  --stock-conc 100 \
  --total-volume 500 \
  --out results/avian_respiratory_panel/pool_sheet.csv

# 4. Export to IDT for ordering.
lamp-forge export \
  --input results/ndv_M_gene/primer_sets.json \
  --format idt \
  --target-label NDV_M \
  --out orders/ndv_M_idt_order.csv
```

**Estimate LOD before ordering** (cloacal swab in 1 mL PBS, 50% RNA
extraction, 100 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5
```

Effective sample = 1000 x 0.50 x (5/100) = 25 uL per reaction ->
LOD_95 approx 120 copies/mL of PBS.  NDV shedding in acutely infected
chickens reaches 10^6-10^9 EID50/mL in cloacal swabs -- orders-of-magnitude
above the achievable LOD.  For sentinel surveillance of subclinical shedding,
increase sample volume to 2 mL or double the reaction input (10 uL) to push
LOD_95 below 60 copies/mL.

**Biosafety note.**  Velogenic NDV (vNDV, Class II genotype V and VII) is a
USDA-regulated select agent (7 CFR Part 331) in the United States.  Wet-lab
validation of an NDV RT-LAMP assay must use inactivated positive controls
(heat-inactivated virus or synthetic RNA standards) unless performed in a
certified facility with appropriate select-agent registration.  The in silico
primer design produced by this pipeline does not require select-agent
handling.

A positive result in the field should trigger immediate reporting to the
national veterinary authority and confirmatory testing at a WOAH Reference
Laboratory (e.g. USDA APHIS NVSL, Plum Island Animal Disease Center).
Position this assay as a **triage/screening** tool, not the confirmatory test.

**Multiplex with ASFV, PRRSV, and AIV** for a complete mixed-species farm
biosecurity panel:

```bash
lamp-forge panel \
  --set IAV=results/influenza_a_M/primer_sets.json \
  --set NDV=results/ndv_M_gene/primer_sets.json \
  --set PRRSV=results/prrsv_ORF7/primer_sets.json \
  --set ASFV=results/asfv_p72/primer_sets.json \
  --top-per-target 5 \
  --out results/farm_biosecurity_4plex
```

Note: PRRSV and NDV are both RNA viruses but from different families
(Arterivirus vs. Paramyxovirus), and ASFV is a DNA virus -- the pooled
reaction uses a one-step RT + LAMP master mix (e.g. NEB WarmStart LAMP Kit
2.0 + WarmStart RTx), which amplifies both RNA and DNA targets in a single
tube at 63-65 degC.

**References.**

- Wise MG, Suarez DL, Seal BS et al. (2004). Development of a real-time
  reverse-transcription PCR for detection of Newcastle disease virus RNA in
  clinical samples. *J Clin Microbiol* 42:329-338.
  doi:10.1128/JCM.42.1.329-338.2004
- Tsai HJ, Chang KH, Tseng CH et al. (2016). Loop-mediated isothermal
  amplification assay for detection of Newcastle disease virus. *J Virol
  Methods* 235:88-95. doi:10.1016/j.jviromet.2016.05.013
- Dimitrov KM, Afonso CL, Yu Q & Miller PJ (2017). Newcastle disease vaccines
  -- A solved problem or a continuous challenge? *Vet Microbiol* 206:126-136.
  doi:10.1016/j.vetmic.2016.12.019
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 19 -- Universal 16S rRNA sample-adequacy control (all panels)

**Goal.** Design a LAMP assay that amplifies any bacterium present in the
sample -- serving as an internal positive control that confirms DNA extraction
succeeded and sufficient microbial nucleic acid reached the reaction.

**Why it's needed.** Every BioVind functional-guild panel (oilfield souring,
farm biosecurity, point-of-care) can return an all-negative result for two
completely different reasons: (1) the target organisms are genuinely absent,
or (2) DNA extraction failed.  Without an internal control you cannot
distinguish these cases and will mis-report extraction failures as clean
samples.  A 16S positive control resolves the ambiguity:

| 16S result | Functional-guild channels | Interpretation |
|---|---|---|
| + | at least one + | Extraction OK; target guilds detected |
| + | all - | Extraction OK; target guilds genuinely absent |
| - | all - | **Extraction failure** -- repeat extraction; do not report |
| - | at least one + | Internal inconsistency; investigate (partial inhibition?) |

In the oilfield souring context, this is the sixth channel that completes
the SRB + methanogen + IRB + APB + NRB five-guild panel (Recipes 6, 10,
15, 16, 17).

**Why 16S rRNA.**  The 16S ribosomal RNA gene is present in all bacteria in
multiple copies per cell (4-10 in most species), making it the most
sensitive possible DNA target for confirming bacterial nucleic acid is
present.  Its mosaic structure -- highly conserved flanking blocks flanking
nine variable regions (V1-V9) -- is precisely what makes it possible to
design six primers all landing in invariant sequences while the amplicon
spans a variable region (making the amplicon unambiguously microbial).

**Key design inversion.**  This is the only LAMP-Forge assay where the
design goal is **inclusivity**, not exclusivity.  All other assays set a
low entropy threshold to find the conserved marker gene; here the low
entropy threshold selects positions that match **all** bacteria.  The
off-target panel does not check for false-positives against organisms in
the same guild -- it checks for false-positives against host DNA and Archaea,
which would inflate the bacterial load count.

**Target sequences.**  Use the ready-made config at `config/general_16s.yaml`,
or paste the block below.  The `taxon_id: 2` anchor (Bacteria) retrieves a
phylogenetically diverse set; supplement with explicit accessions spanning
Firmicutes, Actinobacteria, and Bacteroidetes if auto-retrieval is dominated
by Proteobacteria:

```yaml
target:
  name: universal_16S_control
  taxon_id: 2                      # Bacteria -- broadest anchor for universal coverage
  accessions: []                   # extend with representatives of major phyla if needed
  gene: "16S ribosomal RNA"
  max_sequences: 50                # sample major phyla for cross-phylum conservation check
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85    # tighter: bacteria and archaea share rRNA motifs
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.10         # strict: only truly invariant universal positions
  min_region_length: 220

primer:
  tm_min: 60.0                    # DNA target -- standard LAMP (not RT-LAMP)
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
  dir: results/general_16s
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `conservation.entropy_threshold: 0.10` -- the most distinctive parameter in
  this config.  All other LAMP-Forge configs set 0.20-0.35 bits; this assay
  requires truly invariant positions (entropy < 0.05 bits in large 16S databases)
  to guarantee cross-phylum amplification.  The invariant blocks of 16S rRNA
  (Lane et al. 1985) are the same regions targeted by universal PCR primers
  (27F/1492R, 515F/806R).
- `off_targets.min_identity_threshold: 0.85` -- slightly tighter than the
  standard 0.80 because archaeal 16S rRNA shares secondary-structure motifs
  with bacterial 16S.  In an **oilfield** panel, methanogens and Archaeoglobus
  share the sample matrix; in a **bacterial-only** POC panel, archaeal
  cross-reactivity would be a false-positive.
- `tm_min: 60.0` -- ribosomal DNA is a chromosomal DNA target; do **not** raise
  to the RT-LAMP 63 degC floor used for RNA-virus targets.
- `target.taxon_id: 2` -- unique among all LAMP-Forge configs; every other assay
  uses a genus- or species-level anchor.

**Off-target panel.**  Non-bacterial nucleic acids that co-extract from
oilfield and clinical samples:

| File | Source | Why |
|---|---|---|
| `archaeoglobus_fulgidus.fasta` | NC_000917.1 | Sulfate-reducing archaeon; key oilfield off-target; cross-reactivity would miscount archaeal DNA as bacterial load |
| `methanobacterium_formicicum.fasta` | representative RefSeq | Methanogen archaeon sharing oilfield sample matrix with SRB/IRB |
| `sulfolobus_acidocaldarius.fasta` | NC_007181.1 | Thermoacidophilic crenarchaeon; very distant from Bacteria; strong negative control |
| `human_18S_rRNA.fasta` | NR_003286.4 | Eukaryotic 18S rRNA; must NOT amplify in clinical or PoC contexts |
| `human_chr_fragment.fasta` | Subset of GRCh38 chr1 | Human gDNA background; reservoir inspector or patient sample contamination |

**Expected output.**  With `entropy_threshold: 0.10`, expect 2-5 conserved
windows -- fewer than for specific functional-gene markers, because the
invariant 16S blocks are short (~50-100 nt each).  Top-scoring sets should
show zero flagged hits against human and eukaryotic DNA.  Archaeal 16S may
flag at high identity (>85%) -- document this in the run report and decide
whether archaeal cross-reactivity is acceptable for your panel's context.

**Run.**

```bash
docker compose run --rm lamp-forge run --config /work/config/general_16s.yaml
```

**Complete the six-channel oilfield MIC panel.**  After designing all six
assays independently:

```bash
lamp-forge panel \
  --set SRB=results/srb_dsrB/primer_sets.json \
  --set MCR=results/methanogen_mcrA/primer_sets.json \
  --set IRB=results/irb_omcA/primer_sets.json \
  --set APB=results/apb_fthfs/primer_sets.json \
  --set NRB=results/nrb_narG/primer_sets.json \
  --set 16S=results/general_16s/primer_sets.json \
  --top-per-target 5 \
  --out results/mic_6plex_panel
```

Generate the pooling sheet (six targets = 264 uM minimum; request 300 uM
resuspension from IDT or Twist, or use vendor pre-mix):

```bash
lamp-forge pool \
  --panel results/mic_6plex_panel/panel.json \
  --stock-conc 300 \
  --total-volume 500 \
  --out results/mic_6plex_panel/pool_sheet.csv
```

Export the 16S control primers for IDT ordering:

```bash
lamp-forge export \
  --input results/general_16s/primer_sets.json \
  --format idt \
  --target-label 16S_ctrl \
  --out orders/16S_ctrl_idt_order.csv
```

**Six-channel result interpretation table.**

| SRB | MCR | IRB | APB | NRB | 16S | Interpretation |
|---|---|---|---|---|---|---|
| - | - | - | - | - | - | **Extraction failure** -- repeat extraction; do not report |
| - | - | - | - | - | + | Clean sample; no MIC guilds detected; low risk |
| + | + | + | + | + | + | All five MIC guilds active; maximum corrosion/souring risk |
| + | - | + | + | + | + | APB-SRB syntrophic loop active; iron cycling; escalate biocide |
| - | - | - | - | + | + | NRB dominant; nitrate injection suppressing SRB; treatment working |
| + | + | + | + | + | - | Inconsistent -- possible assay inhibition; investigate |

**LOD context for the 16S control.**  16S rRNA is present in 4-10 copies per
cell, so the effective LOD is 2-5x lower than for single-copy functional genes.
Use `lamp-forge lod` with the same extraction parameters as the functional
channels to model the worst case:

```bash
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5
```

For single-copy genes, LOD_95 approx 150 copies/mL.  With 4-10 copies per
cell for 16S, the control channel calls positive at 15-38 cells/mL -- well
below the 10^3 cells/mL threshold at which oilfield MIC risk is actionable.

**Wet-lab notes.**

- Run no-template controls (NTCs) alongside every batch; 16S amplicon carryover
  is the primary false-positive risk because the gene is environmentally abundant.
- For oilfield samples use bead-capture DNA extraction (e.g. PowerSoil Pro)
  rather than simple boiling; the 16S channel is the canary for extraction
  failures and must be the most sensitive assay in the panel.
- For clinical contexts, validate that the chosen primer set does not amplify
  human 18S or mitochondrial 16S rRNA under standard LAMP conditions.

**References.**

- Lane DJ, Pace B, Olsen GJ et al. (1985). Rapid determination of 16S
  ribosomal RNA sequences for phylogenetic analyses. *Proc Natl Acad Sci USA*
  82:6955-6959. doi:10.1073/pnas.82.20.6955
- Parada AE, Needham DM & Fuhrman JA (2016). Every base matters: assessing
  small subunit rRNA primers for marine microbiomes with mock communities,
  time series and global field samples. *Environ Microbiol* 18:1403-1414.
  doi:10.1111/1462-2920.13023
- Woese CR & Fox GE (1977). Phylogenetic structure of the prokaryotic domain:
  the primary kingdoms. *Proc Natl Acad Sci USA* 74:5088-5090.
  doi:10.1073/pnas.74.11.5088
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 20 — Time-to-positive estimation (`lamp-forge ttp`)

**Goal.** Before ordering primers, confirm that a designed LAMP assay will
turn positive within the device run window (30-60 min on BioVind BioID) for
the expected copy count in a real field or clinical sample.

**Why it matters for BioVind.** BioVind's BioID device has a fixed 30-60 min
run window.  An assay that reads out in 62 min at the diagnostic threshold is
clinically useless on that device — even if it is thermodynamically perfect.
`lamp-forge ttp` surfaces this risk in silico, before a single primer is
ordered.

The model is the empirically validated linear log10 relationship between
threshold time and initial copy count:

    TTP(N) = TTP_1copy - slope × log10(N)

Default parameters come from Bst 2.0 WarmStart at 65 °C (Tanner 2012,
NEB E1700 manual):

    TTP_1copy = 55 min   (conservative 95th-percentile single-molecule TTP)
    slope     = 6 min/decade of copies

---

### 20a — SRB dsrB assay (oil and gas, DNA-LAMP)

**Context.** The SRB assay from Recipe 6 targets dsrB, a single-copy gene.
Produced-water samples contain ~10^3 to 10^6 SRB cells/mL.  After the
standard extraction (1 mL sample, 50% efficiency, 100 uL eluate, 5 uL to
reaction) this maps to 25-25 000 copies per reaction.  Verify all these
copy counts yield TTP < 60 min:

```bash
lamp-forge ttp \
  --preset dna-lamp \
  --device-window 60 \
  --copies-min 1 \
  --copies-max 100000 \
  --out-csv results/srb_dsrB/ttp_dna.csv
```

**Expected output (abridged):**

```
TTP model: preset=dna-lamp, TTP@1cp=55.0 min, slope=6.0 min/decade, min_TTP=10.0 min
Device window: 60 min

    Copies/rxn   TTP (min)   vs window
------------------------------------------
           1.0        55.0        PASS
          10.0        49.0        PASS
         100.0        43.0        PASS
       1 000.0        37.0        PASS
      10 000.0        31.0        PASS
     100 000.0        25.0        PASS
```

All target copy counts are within the 60-min window.  The assay is safe to
order for BioVind deployment.

**LOD + TTP combined check (the standard pre-order workflow):**

```bash
# Step 1: LOD (what is the minimum detectable concentration in sample?)
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5

# Step 2: TTP (will the assay read out in time at that concentration?)
lamp-forge ttp --preset dna-lamp --device-window 60
```

---

### 20b — PRRSV ORF7 assay (farm biosecurity, RT-LAMP)

**Context.** A PRRSV assay on an RNA target uses one-step RT-LAMP with NEB
RTx reverse transcriptase and Bst 2.0 WarmStart at 63-65 °C.  The reverse-
transcription lag adds ~5 min to TTP at low copy counts.  Use the `rt-lamp`
preset:

```bash
lamp-forge ttp \
  --preset rt-lamp \
  --device-window 60 \
  --copies-min 1 \
  --copies-max 1000000 \
  --out-csv results/prrsv_orf7/ttp_rtlamp.csv
```

**Expected output (abridged):**

```
TTP model: preset=rt-lamp, TTP@1cp=60.0 min, slope=6.0 min/decade, min_TTP=12.0 min
Device window: 60 min

    Copies/rxn   TTP (min)   vs window
------------------------------------------
           1.0        60.0        PASS
          10.0        54.0        PASS
         100.0        48.0        PASS
```

At 1 copy/reaction TTP equals the 60-min window exactly (borderline PASS).
If the platform real run time is 55 min, tighten the primer Tm floor to
63.5 °C in the config to speed the RT step and recheck.

---

### 20c — Custom chemistry parameters

If you have measured TTP on your own kit from a calibration curve with
synthetic template, override the defaults:

```bash
lamp-forge ttp \
  --ttp-one-copy 48 \
  --slope 5.5 \
  --device-window 45 \
  --copies-min 10 \
  --copies-max 1000000 \
  --out-csv results/custom_ttp.csv
```

---

**References.**

- Tanner NA, Zhang Y & Evans TC Jr (2012). Visual detection of isothermal
  nucleic acid amplification using pH-sensitive dyes. *BioTechniques*
  53(2):81-89. doi:10.2144/0000113902
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63
- Dao Thi VL et al. (2020). A colorimetric RT-LAMP assay and LAMP-sequencing
  for in-field diagnosis of SARS-CoV-2. *Science* 370(6518):914-917.
  doi:10.1126/science.abc7075
- NEB WarmStart(R) LAMP Kit protocol (E1700). New England Biolabs, 2023.

---

## Recipe 21 -- Group A Streptococcus (GAS) via *speB* (human point-of-care pharyngitis)

**Goal.** Detect *Streptococcus pyogenes* (Group A Streptococcus, GAS) from a
throat swab at the point of care -- fast enough to drive an antibiotic-
prescribing decision before the patient leaves the clinic or teleconsultation.

**Why it matters for BioVind.** GAS pharyngitis is the most common bacterial
cause of sore throat worldwide and the #1 reason for antibiotic prescription
in primary care globally.  Crucially, untreated GAS carries sequelae
(rheumatic fever, rheumatic heart disease, post-streptococcal
glomerulonephritis) that are preventable with a 10-day penicillin course --
but only if the diagnosis is made correctly.  The current standard of care
is the rapid antigen detection test (RADT), which has ~86% sensitivity;
missed cases proceed to 48-hour throat culture.

A 30-minute LAMP assay on a throat swab closes this gap:

- **Sensitivity >= 95%** in comparable throat-swab LAMP studies
  (Kitagawa et al. 2011; Li et al. 2019).
- **DNA target**: no RT step, no cold chain for the LAMP master mix.
- **Single-tube read-out** by turbidity, pH dye, or fluorescence --
  compatible with BioVind BioID portable device.
- **Decision time < 30 min** from swab to result: well within the
  urgent-care or teleconsultation workflow.

This recipe is the first in the human POC vertical that targets a bacterium
rather than an RNA virus.  BioVind's third vertical (rural urgent care,
telemedicine) is where portable bacterial pharyngitis testing earns its
commercial case: shipping throat-swab swabs to a remote lab for culture costs
48 hours and $50-100 per test; a portable LAMP assay on the same swab in
real time supports immediate antibiotic stewardship decisions.

**Why *speB*.** *speB* (streptococcal cysteine protease, also called
*streptococcal pyrogenic exotoxin B* or SCP) encodes a ~38 kDa zymogen
that is:

- Present in **all** *S. pyogenes* strains (no known natural deletion).
- **Absent from other streptococcal species**, including Group B
  (*S. agalactiae*), Group C/G (*S. dysgalactiae*), and *S. pneumoniae*
  -- the intra-genus off-targets most likely to appear on a throat swab.
- Under strong pathogenicity-island selection: the catalytic Cys/His dyad
  and prodomain cleavage site are near-invariant across M-types, yielding
  at least two LAMP-accessible conserved windows in the ~1.4 kb CDS.
- The basis of multiple published LAMP assays, confirming empirical
  primer-design feasibility.

**Target sequences.** NCBI taxon 1314 (*S. pyogenes*), gene `speB`,
25 sequences spans the globally dominant M-types (M1, M3, M12, M28, M89,
M77, M4) needed for a conserved-region analysis that reflects clinical
diversity.

Use the ready-made config at `config/gas_speB.yaml`, or paste the block
below:

```yaml
target:
  name: gas_speB_poc
  taxon_id: 1314          # Streptococcus pyogenes (Group A Streptococcus)
  gene: speB              # cysteine protease exotoxin B; GAS-specific marker
  max_sequences: 25       # span key M-types: M1, M3, M12, M28, M89, M77, M4
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- intra-genus streptococci share GC content
  min_coverage_threshold: 0.85

conservation:
  window_size: 30
  entropy_threshold: 0.20        # speB is well-conserved within GAS (similar to rpoB)
  min_region_length: 220

primer:
  tm_min: 60.0                   # DNA target -- standard LAMP (no RT step needed)
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 30.0                   # S. pyogenes is AT-rich (~38.5% GC); floor at 30%
  gc_max: 55.0                   # cap below GC-rich streptococci / staphylococci
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/gas_speB
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `gc_min: 30.0` -- S. pyogenes is one of the most AT-rich beta-haemolytic
  streptococci (~38.5% genomic GC; speB CDS is ~37% GC).  Most other
  LAMP-Forge configs use gc_min >= 35; lowering to 30 is required here to
  find primer candidates in the lower-GC codon positions without missing
  conserved windows.
- `gc_max: 55.0` -- caps primers well below GC-rich common off-targets
  (S. pneumoniae ~40%, S. aureus ~33%, human DNA ~41%).  A 55% ceiling
  is deliberately lower than the default 65% to reduce the chance of
  designing primers that also match GC-richer respiratory flora.
- `entropy_threshold: 0.20` -- speB is under strong pathogenicity-island
  selection comparable to single-copy housekeeping genes (rpoB, invA).
  This is tighter than the functional-gene markers used for polyphyletic
  targets (dsrB/mcrA at 0.35); the speB catalytic domains are near-
  invariant across M-types and the threshold exploits this.
- `min_identity_threshold: 0.85` and `min_coverage_threshold: 0.85` --
  intra-genus streptococcal off-targets share some conserved gene blocks;
  the stricter 0.85/0.85 thresholds (versus the bacterial default 0.80/0.80)
  are required to surface any cross-reactive primer candidates before
  wet-lab validation.
- `tm_min: 60.0` -- DNA target; do **not** raise to the RT-LAMP 63 degC
  floor used for PRRSV, FMDV, AIV, and NDV.

**Off-target panel.** The key differentials are intra-genus streptococci and
common pharyngeal flora:

| File | Source | Why |
|---|---|---|
| `strep_agalactiae.fasta` | NEM316 RefSeq (NC_004116.1) | Group B Streptococcus; same genus; no speB; must not cross-flag |
| `strep_pneumoniae.fasta` | D39 RefSeq (NC_008533.2) | Pneumococcus; common pharyngeal/respiratory coloniser |
| `strep_salivarius.fasta` | CCHSS3 RefSeq | Oral commensal; abundant on throat swabs |
| `strep_mutans.fasta` | UA159 RefSeq (NC_004350.2) | Oral streptococcus; GC content (~36.8%) close to GAS |
| `staph_aureus.fasta` | MRSA252 RefSeq (NC_002952.2) | Common pharyngeal pathogen; mimics GAS on clinical exam |
| `human_chr_fragment.fasta` | GRCh38 chr1 subset | Host DNA from throat swab; dominant background nucleic acid |

**Expected output.** The speB CDS (~1.4 kb) should yield 2-4 conserved
windows at the 0.20-bit threshold.  Top-scoring sets should localise to the
pro-domain and zymogen regions, with zero flagged hits against all off-target
streptococci at the 0.85 x 0.85 threshold.  *S. aureus* cross-reactivity at
any identity is biologically unacceptable and would be a false-positive in
the clinical context; confirm it is zero before ordering.

**Run.**

```bash
docker compose run --rm lamp-forge run --config /work/config/gas_speB.yaml
```

**Estimate LOD before ordering** (throat swab in 200 uL lysis buffer, 50%
DNA extraction, 50 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 200 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 200 x 0.50 x (5/50) = 10 uL per reaction.
LOD_95 approx 300 copies/mL of lysis buffer, equivalent to approximately
60 GAS cells per swab -- consistent with published LAMP sensitivity
(Kitagawa et al. 2011 reported LOD of 10 CFU/reaction on synthetic template).
GAS carriage on throat swabs during active pharyngitis typically exceeds
10^4 CFU/mL; the assay LOD provides two orders-of-magnitude headroom.

**Export primers for ordering.**

```bash
lamp-forge export \
  --input results/gas_speB/primer_sets.json \
  --format idt \
  --target-label GAS_speB \
  --out orders/gas_speB_idt_order.csv
```

**Verify config before running.**

```bash
lamp-forge validate \
  --config config/gas_speB.yaml \
  --no-check-dirs
```

**Wet-lab notes.**

- Throat swab (posterior pharynx + tonsillar pillars) is the standard
  sample type.  Swab directly into lysis buffer or a universal transport
  medium (UTM) for same-day processing.
- The dominant PCR inhibitor in throat swabs is mucinous glycoproteins.
  A simple bead-beating or CTAB extraction removes most inhibition; avoid
  boiling alone as it leaves protein aggregates that can inhibit Bst.
- Validate the primer set against a mixed-organism throat-swab matrix
  spiked with GAS reference strain ATCC 12344 (M6) at 10, 100, and 1000
  CFU/swab to confirm LOD in the clinical matrix.
- Run a no-template control (NTC) and a positive extraction control
  (purified GAS DNA) alongside every batch.
- A positive result should be reported as a presumptive positive pending
  culture confirmation per local clinical guideline.  For paediatric
  patients (<15 years) with pharyngitis score >= 3 (Centor / McIsaac),
  treat empirically without waiting for culture confirmation per IDSA 2012.

**References.**

- Kitagawa Y, Ueno M, Shinozuka N et al. (2011). Loop-mediated isothermal
  amplification for rapid and sensitive detection of *Streptococcus pyogenes*.
  *J Infect Chemother* 17(4):486-493. doi:10.1007/s10156-010-0189-3
- Li J, Macdonald J & von Stetten F (2019). Review: a comprehensive summary of
  a decade development of the recombinase polymerase amplification. *Analyst*
  144(1):31-67. doi:10.1039/C8AN01621F
- Shulman ST, Bisno AL, Clegg HW et al. (2012). Clinical practice guideline
  for the diagnosis and management of group A streptococcal pharyngitis: 2012
  update by the IDSA. *Clin Infect Dis* 55(10):e86-e102.
  doi:10.1093/cid/cis629
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 22 -- Clostridioides difficile via *tcdB* (hospital point-of-care CDI)

**Goal.** A LAMP assay that detects toxigenic *Clostridioides difficile* DNA
from a stool sample, enabling a cohort isolation and antibiotic-prescribing
decision at the point of care without shipping the sample to a centralised
laboratory.

**Why it matters.** *C. difficile* infection (CDI) is the leading cause of
hospital-acquired infectious diarrhoea in high-income countries (~230,000
hospitalisations per year in the US; CDC HAI data 2021).  The standard
diagnostic pathway suffers a critical delay: enzyme immunoassay (EIA) for
toxin A/B is rapid but only ~75% sensitive; nucleic-acid amplification tests
(NAATs) achieve >95% sensitivity but require a centralised laboratory and
1-4 hours.  That delay forces hospitals to manage symptomatic patients under
precautionary contact precautions for hours while the result is awaited -- a
ward-management bottleneck that drives infection transmission and bed occupancy
costs.  A 30-minute LAMP assay from a liquid stool sample collapses that delay
to within a clinical encounter, enabling immediate cohort isolation or
clearance of contact precautions -- a direct fit for BioVind's human
point-of-care vertical targeting urgent care, rural hospitals, and
telemedicine sample-collection sites.

**Why *tcdB* and not *tcdA*.**
The IDSA/SHEA 2017 CDI guidelines (McDonald et al. 2018) and ESCMID 2021
update both recommend *tcdB*-targeted molecular tests because:

- A subset of epidemic ribotypes (RT017, some RT033) carry *tcdB* but have a
  natural deletion of *tcdA* (the so-called "A-B+" phenotype).  These strains
  produce disease that is indistinguishable clinically from A+B+ strains; a
  *tcdA*-only assay misses them.
- *tcdB* encodes the major cytotoxin responsible for intestinal epithelial
  damage and inflammation; its presence is the clinically actionable
  correlate of toxigenic *C. difficile* carriage.
- Published *tcdB* LAMP assays achieve LOD_95 <= 25 copies/reaction in
  stool-matrix studies (Shin et al. 2016).

**AT-rich genome challenge.**
*C. difficile* has one of the most AT-rich genomes among Gram-positive
pathogens: ~28.6% GC genome-wide (Sebaihia et al. 2006), lower even than
*S. pyogenes* (~38.5%).  The *tcdB* gene and the surrounding pathogenicity
locus (PaLoc) are ~26-30% GC.  Standard LAMP primer GC floors of 40-45%
used for GC-richer organisms would yield no candidates.  This recipe lowers
`gc_min` to 25%, consistent with published LAMP designs for other low-GC
clostridia (Nakagawa et al. 2010 for *C. perfringens*), while capping
`gc_max` at 50% to avoid cross-matching GC-richer gut-flora off-targets
(*Bacteroides* spp. ~43-48% GC, *E. coli* K-12 ~51% GC).

**Conserved region choice.**
The N-terminal glucosyltransferase catalytic domain of TcdB (~nt 1-1100 of
the 7.1 kb CDS) is under strong purifying selection at the catalytic DXD
motif (Asp286, Asp288) and tryptophan-rich substrate-binding loops.  This
~1.1 kb window is broadly conserved across all major clinical ribotypes and
is the target used in published LAMP (Shin 2016) and RT-PCR (Crobach 2016)
assays.  Set `max_sequences: 25` to retrieve sequences from reference strains
630 (RT012), R20291 (RT027/NAP1), and diverse clinical isolates spanning
RT001, RT078, RT014, RT017.

**Target sequences.** NCBI taxon ID 1496 (*Clostridioides difficile*), gene
`tcdB`, max 25 sequences retrieves a ribotype-diverse clinical isolate set.

**Off-target panel.** Drop these FASTAs into `input/off_targets/`:

| File | Source | Why |
|---|---|---|
| `c_sordellii.fasta` | ATCC 9714T RefSeq | Closest functional relative; lethal toxin (LT) with ~40% TcdB protein identity |
| `c_perfringens.fasta` | ATCC 13124 RefSeq (NC_003366.1) | Common gut Clostridium; no tcdB |
| `c_botulinum.fasta` | ATCC 19397 RefSeq (NC_009699.1) | Anaerobic Gram-positive; no tcdB |
| `bacteroides_fragilis.fasta` | NCTC 9343 RefSeq (NC_003228.3) | Dominant stool anaerobe; BFAG enterotoxin is unrelated |
| `e_coli_k12.fasta` | MG1655 RefSeq (NC_000913.3) | Abundant stool commensal; high background in stool LAMP |
| `human_chr_fragment.fasta` | GRCh38 chr1 subset | Host colonocyte DNA shed in stool; dominant nucleic acid in stool samples |

**Config.**

Use the ready-made config at `config/cdiff_tcdB.yaml`, or paste the block
below:

```yaml
target:
  name: cdiff_tcdB_poc
  taxon_id: 1496          # Clostridioides difficile
  gene: tcdB              # toxin B glucosyltransferase; IDSA-recommended target
  max_sequences: 25       # span ribotypes: RT001, RT027 (NAP1), RT106, RT078, RT014
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- tcdB-like domain in C. sordellii LT
  min_coverage_threshold: 0.85

conservation:
  window_size: 30
  entropy_threshold: 0.25        # catalytic core conserved but more inter-ribotype drift
                                 # than housekeeping genes (rpoB, speB)
  min_region_length: 200

primer:
  tm_min: 60.0                   # DNA target -- standard LAMP (no RT step needed)
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 25.0                   # C. diff genome ~28.6% GC (very AT-rich); tcdB
                                 # PaLoc is ~26-30% GC; must drop well below 40%
  gc_max: 50.0                   # cap below GC-richer gut-flora off-targets
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/cdiff_tcdB
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `gc_min: 25.0` -- The most important deviation from standard configs.
  *C. difficile* is one of the most AT-rich Gram-positive pathogens
  (~28.6% GC).  Every LAMP-Forge config for GC-richer organisms uses
  `gc_min: 30-40`.  Lowering to 25% here is essential: without it the
  primer design stage would find zero candidate windows in the tcdB CDS
  and the run would fail with zero primer sets.  Published LAMP assays
  for low-GC clostridia routinely include primers with GC% as low as
  26-30% (Nakagawa 2010; Shin 2016).

- `gc_max: 50.0` -- Caps primers well below the GC-rich gut-flora
  background: *Bacteroides* spp. ~43-48% GC, *E. coli* K-12 ~51% GC.
  Any primer with GC > 50% is more likely to find a binding site in
  stool background organisms than in the AT-rich *C. diff* target.

- `entropy_threshold: 0.25` -- Slightly relaxed compared to well-
  conserved housekeeping genes (rpoB at 0.20, speB at 0.20).  The
  *tcdB* catalytic core is strongly conserved but epidemic ribotypes
  (RT027 vs. RT078 vs. RT001) have more synonymous divergence than
  single-copy chromosomal markers.

- `min_identity_threshold: 0.85` / `min_coverage_threshold: 0.85` --
  The chief off-target is *Clostridioides sordellii* lethal toxin (LT),
  which shares ~40% protein identity with TcdB.  The 0.85/0.85
  threshold surfaces any residual cross-reactive primer candidates before
  ordering.

- `tm_min: 60.0` -- *C. difficile* tcdB is a chromosomal DNA target; do
  **not** raise to the RT-LAMP 63 degC co-activity floor used for RNA
  targets (PRRSV, FMDV, AIV, NDV).

**Run.**

```bash
docker compose run --rm lamp-forge run --config /work/config/cdiff_tcdB.yaml
```

**Estimate LOD before ordering** (liquid stool 200 uL in 800 uL lysis
buffer, ~40% DNA extraction efficiency, 100 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 200 \
  --efficiency 0.40 \
  --eluate-volume 100 \
  --reaction-input 5
```

Effective sample = 200 x 0.40 x (5/100) = 4 uL per reaction.
LOD_95 approx 750 copies/mL of stool input, comfortably within the range
reported by Shin et al. 2016 (~25 copies/reaction on extracted DNA).
Symptomatic CDI patients typically shed 10^5-10^8 *C. diff* spores/mL;
the assay LOD provides three or more orders-of-magnitude headroom.

**Check TTP for the 60-min BioVind device window.**

```bash
lamp-forge ttp \
  --copies 750 \
  --preset lamp
```

At 750 copies/reaction the predicted TTP is ~34 min -- well inside the
60-min BioVind BioID device window.

**Export primers for ordering.**

```bash
lamp-forge export \
  --input results/cdiff_tcdB/primer_sets.json \
  --format idt \
  --target-label CDiff_tcdB \
  --out orders/cdiff_tcdB_idt_order.csv
```

**Verify config before running.**

```bash
lamp-forge validate \
  --config config/cdiff_tcdB.yaml \
  --no-check-dirs
```

**Wet-lab notes.**

- Sample type is liquid (unformed) stool from symptomatic patients.
  Do not submit formed stool -- IDSA guidelines recommend molecular
  testing only on diarrhoeal specimens to reduce asymptomatic carriage
  detection.
- Major LAMP inhibitors in stool: bile salts, fatty acids, and complex
  polysaccharides.  A commercial stool DNA kit (e.g. Qiagen QIAamp DNA
  Stool Mini Kit) with bead-beating removes most inhibition; a 1:10
  dilution of crude lysate in molecular-grade water is the minimum
  effective clean-up when kit processing is not feasible in the field.
- *C. difficile* spores survive drying and are the main infective form.
  Process stool samples in a BSL-2 cabinet; treat all positive extracts
  and amplicons as potentially hazardous.
- Run a no-template control (NTC) and a positive extraction control
  (purified *C. difficile* reference-strain DNA, e.g. ATCC 9689 / strain
  VPI 10463) in every batch.
- Include a 16S rRNA internal control (Recipe 19) to confirm the
  extraction step succeeded on each stool sample; an all-negative result
  may indicate inhibition rather than true absence of *C. diff*.
- A positive tcdB LAMP result is a presumptive positive for toxigenic
  CDI.  Do not use LAMP for test-of-cure: NAAT positivity can persist
  for weeks after resolution of symptoms due to residual DNA from dead
  organisms.

**References.**

- McDonald LC, Gerding DN, Johnson S et al. (2018). Clinical practice
  guidelines for *Clostridium difficile* infection in adults and children:
  2017 update by IDSA and SHEA. *Clin Infect Dis* 66(7):e1-e48.
  doi:10.1093/cid/cix1085
- Shin HB, Yoon J, Lee Y et al. (2016). Evaluation of a loop-mediated
  isothermal amplification assay for rapid detection of *Clostridium
  difficile* toxin B gene in stool specimens. *BMC Infect Dis* 16:371.
  doi:10.1186/s12879-016-1700-z
- Sebaihia M, Wren BW, Mullany P et al. (2006). The multidrug-resistant
  human pathogen *Clostridium difficile* has a highly mobile, mosaic
  genome. *Nat Genet* 38(7):779-786. doi:10.1038/ng1830
- Crobach MJT, Planche T, Eckert C et al. (2016). ESCMID: update of the
  diagnostic guidance document for *Clostridium difficile* infection.
  *Clin Microbiol Infect* 22 Suppl 4:S63-81. doi:10.1016/j.cmi.2016.03.010
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 23 -- Bovine Respiratory Disease (BRD) five-target panel (`lamp-forge bov-risk`)

**Goal.** Build a five-channel multiplex isothermal panel that detects the major
viral and bacterial pathogens driving bovine respiratory disease (BRD) from a
single nasal swab or bronchoalveolar lavage sample -- fast enough to influence
treatment decisions before antimicrobials are administered.

**Why it's interesting.** BRD ("shipping fever") is the leading cause of morbidity
and mortality in beef feedlots and a major cause of dairy production loss worldwide,
with estimated annual losses exceeding USD 900 million in North America alone.  The
disease follows a predictable two-hit pattern: a primary viral pathogen (BRSV, BCoV,
BVDV, or IBR) depresses respiratory immunity, allowing *Mannheimia haemolytica* to
invade the lower airways and cause the acute fibrinous pneumonia that kills animals.
Early identification of the viral trigger and bacterial co-infection shapes the
treatment decision: viral-only BRD may not require antimicrobials; viral + bacterial
co-infection requires immediate antimicrobial intervention.

A BioVind-style portable platform delivering a five-plex result at the pen-side in
under 60 minutes enables this treatment-decision window in a way that laboratory
culture (24-72 hours) cannot.

**Panel targets.**

| Target | Gene | NA type | Config | bov-risk flag |
|---|---|---|---|---|
| BRSV | N (nucleoprotein) | RNA | `config/brsv_N_gene.yaml` | `--brsv` |
| BCoV | N (nucleoprotein) | RNA | `config/bcov_N_gene.yaml` | `--bcov` |
| BVDV (types 1 + 2) | 5'-UTR | RNA | `config/bvdv_5utr.yaml` | `--bvdv` |
| IBR / BoHV-1 | gB (UL27) | DNA | `config/ibr_gB.yaml` | `--ibr` |
| *Mannheimia haemolytica* | lktA | DNA | `config/mhae_lktA.yaml` | `--mhae` |

Three targets are RNA viruses (BRSV, BCoV, BVDV) requiring one-step RT-LAMP at
63-65 degC (NEB RTx + Bst 2.0 WarmStart).  Two targets are DNA (IBR, MHAE) and
run at the standard 60-65 degC LAMP window.  Use a single reaction temperature of
63-65 degC for all five when running as a single multiplexed tube: the DNA-LAMP
targets function well at this temperature.

**Step 1 -- Design each assay independently.**

```bash
# RNA virus targets (RT-LAMP: tm_min 63 degC in all three configs)
docker compose run --rm lamp-forge run --config /work/config/brsv_N_gene.yaml
docker compose run --rm lamp-forge run --config /work/config/bcov_N_gene.yaml
docker compose run --rm lamp-forge run --config /work/config/bvdv_5utr.yaml

# DNA targets (standard DNA-LAMP: tm_min 60 degC)
docker compose run --rm lamp-forge run --config /work/config/ibr_gB.yaml
docker compose run --rm lamp-forge run --config /work/config/mhae_lktA.yaml
```

**Step 2 -- Verify RT-LAMP readiness for the three RNA targets.**

```bash
lamp-forge rt-check \
  --input results/brsv_N_gene/primer_sets.json --na-type rna \
  --out-csv results/brsv_N_gene/rt_check.csv

lamp-forge rt-check \
  --input results/bcov_N_gene/primer_sets.json --na-type rna

lamp-forge rt-check \
  --input results/bvdv_5utr/primer_sets.json --na-type rna
```

Sets marked **NOT OPTIMIZED** have core primers below 63 degC.  Re-run with
`primer.tm_min: 63.0` (already set in all three RNA configs) and re-check.

**Step 3 -- Screen for cross-assay primer compatibility.**

```bash
lamp-forge panel \
  --set BRSV=results/brsv_N_gene/primer_sets.json \
  --set BCOV=results/bcov_N_gene/primer_sets.json \
  --set BVDV=results/bvdv_5utr/primer_sets.json \
  --set IBR=results/ibr_gB/primer_sets.json \
  --set MHAE=results/mhae_lktA/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/brd_5plex_panel
```

**Step 4 -- Calculate the pooling sheet.**

Five targets = 30 primers total.  Minimum required stock concentration =
5 x 44 uM = 220 uM; request 250 uM resuspension from IDT or Twist.

```bash
lamp-forge pool \
  --panel results/brd_5plex_panel/panel.json \
  --stock-conc 250 \
  --total-volume 500 \
  --out results/brd_5plex_panel/pool_sheet.csv
```

**Step 5 -- Export to vendor order.**

```bash
lamp-forge panel-export \
  --panel results/brd_5plex_panel/panel.json \
  --format idt \
  --out orders/brd_5plex_idt_order.csv
```

**Step 6 -- Estimate assay LOD before ordering.**

Nasal swab in 500 uL PBS, 50% extraction, 50 uL eluate, 5 uL to reaction:

```bash
lamp-forge lod \
  --sample-volume 500 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 500 x 0.5 x (5/50) = 25 uL per reaction -> LOD_95 approx
120 copies/mL.  BRSV loads in nasal secretions during acute febrile phase reach
10^5-10^7 copies/mL, providing orders-of-magnitude headroom.

**Step 7 -- Interpret BRD panel results.**

```bash
# BRSV + M. haemolytica co-detected (highest-mortality BRD pattern):
lamp-forge bov-risk --brsv --mhae

# BVDV detected alone (screen herd for persistently infected animals):
lamp-forge bov-risk --bvdv

# Full panel from a JSON flags file (e.g. output from a portable device):
lamp-forge bov-risk --input-json results/brd_flags.json \
  --out-json results/brd_assessment.json
```

Key clinical patterns interpreted by `lamp-forge bov-risk`:

| Flags | Alert | Interpretation |
|---|---|---|
| BRSV + MHAE | CRITICAL | Viral-bacterial co-infection; initiate antimicrobials immediately |
| BVDV alone | HIGH | Screen herd for PI animals; movement restriction; notify veterinarian |
| IBR alone | MODERATE | Mandatory notification in EU/UK/Scandinavia; check programme rules |
| BCoV alone | MODERATE | Supportive care; monitor for Mannheimia secondary infection |

**Target biology notes.**

*BRSV.* Bovine orthopneumovirus (Family Pneumoviridae, negative-sense ssRNA) has
two subgroups (A and B; ~80-85% N-gene nt identity).  The N-gene config includes
representative accessions from both subgroups; verify from NCBI and expand to
>= 20 diverse field strains before running.

*BCoV.* Bovine coronavirus (Betacoronavirus 1, positive-sense ssRNA) is closely
related to human HCoV-OC43 (~95-96% N-gene nt identity).  Review the specificity
heatmap for OC43 hits and prefer primer windows in BCoV-specific N-gene regions.

*BVDV.* Pestivirus bovis (positive-sense ssRNA; Family Flaviviridae) has two types
(BVDV-1: NCBI TaxID 11099; BVDV-2: NCBI TaxID 97012) sharing ~75-80% 5'-UTR identity.
Add BVDV-2 complete genome accessions to `bvdv_5utr.yaml` before running for
pan-BVDV coverage.

*IBR / BoHV-1.* Bovine alphaherpesvirus 1 (dsDNA, ~72% GC) causes latent infection.
The gB config sets gc_min=50% / gc_max=75% for the GC-rich UL27 gene.  Standard
DNA-LAMP (no RT); IBR detection may trigger mandatory reporting in EU/UK/Scandinavia.

*Mannheimia haemolytica.* Gram-negative Pasteurellaceae (~41% GC); the only bacterial
channel.  lktA encodes the RTX leukotoxin; it is absent from *Pasteurella multocida*,
making it a highly specific M. haemolytica marker.  Standard DNA-LAMP.

**Wet-lab notes.**

- Use a combined RNA/DNA extraction kit to capture both RNA virus targets and the
  bacterial DNA target from a single nasal swab eluate.
- Validate the five-plex at 63 degC (the RT-LAMP floor) across all five channels
  simultaneously using synthetic RNA/DNA standards.
- Include a 16S rRNA internal control (Recipe 19) to confirm extraction success.
- IBR detection requires confirmatory serology at an accredited laboratory before
  mandatory reporting obligations are triggered.

**References.**

- USDA-NAHMS (2011) Feedlot 2011, Part I: Baseline reference of feedlot management.
  USDA APHIS Veterinary Services.
- Valarcher JF et al. (1999) Evolution of bovine respiratory syncytial virus.
  *J Virol* 74:10714-10728. doi:10.1128/JVI.74.22.10714-10728.2000
- Saif LJ (2004) Animal coronaviruses: what can they teach us about SARS?
  *Rev Sci Tech* 23:643-660.
- Ridpath JF (2010) Bovine viral diarrhea virus: global status. *Vet Clin North
  Am Food Anim Pract* 26:105-121. doi:10.1016/j.cvfa.2009.10.007
- Highlander SK et al. (1989) Secretion and virulence of *Pasteurella haemolytica*
  leukotoxin. *J Bacteriol* 171:1862-1872.
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 24 -- *Mycoplasma bovis* via *uvrC* (BRD supplemental 6th channel)

**Goal.** Detect *Mycoplasma bovis* in nasal swab, bronchoalveolar lavage
(BAL), or joint fluid from cattle and calves, adding a mycoplasmal sixth
channel to the five-target BRD panel (Recipe 23) on a BioVind-style
portable platform.

**Why it's interesting.**  *Mycoplasma bovis* is the most pathogenic
*Mycoplasma* in cattle and an increasingly urgent veterinary problem.
It causes:

- **Chronic fibrinous pneumonia** in feedlot calves (responsible for
  30-50% of feedlot pneumonia cases that fail standard antimicrobial therapy)
- **Middle-ear disease** (otitis media) in young calves, a leading
  cause of calf mortality in dairy operations
- **Septic arthritis** -- often poly-articular, refractory to treatment
- **Sub-clinical and clinical mastitis** in dairy cattle

The critical clinical distinction from the five BRD panel targets in
Recipe 23 is therapeutic: *M. bovis* is intrinsically resistant to
beta-lactams and aminoglycosides (it lacks a cell wall), and pan-resistant
strains failing all frontline antimicrobials (enrofloxacin, tulathromycin,
florfenicol) are now routinely isolated in North America and Europe.
Detecting *M. bovis* specifically -- rather than assuming the causative agent
is *M. haemolytica* -- allows a directed decision: isolate the animal,
perform susceptibility testing before committing to therapy, and apply
enhanced biosecurity.

*M. bovis* is a **DNA target** (no reverse transcription required).
Standard DNA-LAMP (Bst 2.0 at 60-65 degC) delivers a result in under
30 minutes, fully compatible with a BioVind-style portable platform.

**Target gene: *uvrC*.** The UvrC subunit of the nucleotide excision repair
complex (excinuclease ABC) is the published LAMP target of choice for *M. bovis*
(Chen et al. 2017, *J Vet Sci*).  It is:

- Present as a single copy in all *M. bovis* genomes sequenced to date
- Sequence-divergent from its orthologues in bovine respiratory *Mycoplasma* spp.
  sharing the same niche (*M. bovirhinis*, *M. dispar*, *M. bovoculi*)
- Under strong purifying selection at the ATP-binding and HhH-GPD domains,
  yielding at least two LAMP-accessible conserved windows in the ~1.8 kb CDS
- LOD_95 of 10 copies per reaction demonstrated on extracted DNA in the
  original Chen et al. 2017 validation

**The GC-content challenge.**  *M. bovis* has one of the lowest genomic GC
contents of any livestock pathogen: **~29.5%** (compare *S. pyogenes* ~38.5%,
*C. difficile* ~28.6%).  The standard LAMP-Forge `gc_min` of 40% would yield
zero primer candidates for this organism.  The config below drops `gc_min` to
25% -- the single most important parameter deviation from standard configs.
If the run returns zero primer sets, lower `gc_min` further to 20% and re-run.

**Config.** Use the ready-made config at `config/mbovis_uvrC.yaml`, or paste
the block below:

```yaml
target:
  name: mbovis_uvrC
  taxon_id: 28903          # Mycoplasma bovis (NCBI TaxID)
  accessions: []           # extend with PG45, HB0801, NP151, Ningxia-1 strains
  gene: uvrC               # DNA repair subunit C; alt annotation: "uvrC protein"
  max_sequences: 25        # spans key M. bovis field strains
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- ruminant Mycoplasma spp. share conserved blocks
  min_coverage_threshold: 0.85

conservation:
  window_size: 30
  entropy_threshold: 0.20        # uvrC is conserved within M. bovis (>97% inter-strain)
  min_region_length: 200

primer:
  tm_min: 60.0                   # DNA target -- standard LAMP; NOT RT-LAMP
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 25.0                   # M. bovis genomic GC ~29.5%; CRITICAL deviation from default
  gc_max: 47.0                   # cap below GC-richer ruminant Mycoplasma off-targets
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/mbovis_uvrC
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `gc_min: 25.0` -- The most critical departure from standard configs.  Every
  other LAMP-Forge config uses `gc_min >= 25-35`; the default is 40%.
  *M. bovis* at ~29.5% GC requires primers in the 25-35% GC range to hit the
  conserved catalytic-domain blocks of uvrC.  Without this adjustment the run
  will fail with zero primer candidates.
- `gc_max: 47.0` -- Deliberately lower than the 65% default.  Most bovine
  respiratory *Mycoplasma* off-targets have slightly higher GC than *M. bovis*
  (*M. agalactiae* ~33%, *M. bovirhinis* ~31%); capping at 47% reduces the
  risk of cross-reactive primers binding these species at the GC-rich end.
- `entropy_threshold: 0.20` -- *uvrC* is under strong functional constraint
  within *M. bovis*; within-species identity > 97%.  This is the same tight
  threshold used for single-species housekeeping genes (*rpoB* at 0.20,
  *speB* at 0.20, *lktA* at 0.20).
- `min_identity_threshold: 0.85` -- Tighter than the 0.80 bacterial default.
  *Mycoplasma agalactiae* shares > 80% 16S identity with *M. bovis*; the
  0.85 threshold is required to surface any residual uvrC cross-reactive
  primer candidates before wet-lab validation.
- `tm_min: 60.0` -- *M. bovis* uvrC is a chromosomal DNA target; do **not**
  raise to the RT-LAMP 63 degC floor used for BRSV, BCoV, and BVDV in
  Recipe 23.

**Off-target panel.**  The critical differentials are bovine respiratory
*Mycoplasma* species and common BRD bacterial co-pathogens:

| File | Source | Why |
|---|---|---|
| `mycoplasma_bovirhinis.fasta` | NCBI TaxID 28901 representative | Bovine respiratory commensal Mycoplasma; healthy-cattle coloniser; must not co-flag |
| `mycoplasma_dispar.fasta` | NCBI TaxID 29556 representative | Mild-calf-pneumonia Mycoplasma; shared niche with M. bovis |
| `mycoplasma_agalactiae.fasta` | PG2 RefSeq (NC_009497.1) | Closest relative (>80% 16S identity); goat/sheep agalactia; must not flag in cattle panel |
| `pasteurella_multocida.fasta` | ATCC 43137 RefSeq | Most common BRD bacterial co-pathogen after M. haemolytica; must not cross-flag |
| `bos_taurus_fragment.fasta` | ARS-UCD2.0 chr1 subset | Host DNA from nasal swab or BAL; dominant background nucleic acid |

**Expected output.**  The *uvrC* CDS (~1.8 kb) should yield 1-3 conserved
windows.  Top-scoring sets should have zero flagged hits against all
off-targets at the 0.85 x 0.85 threshold.  If fewer than two conserved
windows appear, check that the sequence retrieval captured strains from
both the North American (HB0801, NP151) and European (PG45, Ningxia-1)
lineages, which have some inter-lineage synonymous divergence.

**Run.**

```bash
docker compose run --rm lamp-forge run --config /work/config/mbovis_uvrC.yaml
```

**Estimate LOD before ordering** (nasal swab in 500 uL PBS, 50% DNA
extraction, 50 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 500 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 500 x 0.5 x (5/50) = 25 uL per reaction ->
LOD_95 approx 120 copies/mL.  In clinically active *M. bovis* pneumonia
nasal secretions contain 10^4-10^6 genomic equivalents/mL, providing
two or more orders-of-magnitude headroom.

**Check TTP for the 60-min BioVind BioID device window.**

```bash
lamp-forge ttp --preset dna-lamp --device-window 60
```

At 100 copies/reaction (approx LOD_95 for a 25 uL effective sample) the
predicted DNA-LAMP TTP is ~43 min -- well inside the 60-min BioVind
BioID device window.

**Export primers for ordering.**

```bash
lamp-forge export \
  --input results/mbovis_uvrC/primer_sets.json \
  --format idt \
  --target-label MBovis_uvrC \
  --out orders/mbovis_uvrC_idt_order.csv
```

**Verify config before running.**

```bash
lamp-forge validate \
  --config config/mbovis_uvrC.yaml \
  --no-check-dirs
```

**Add M. bovis to the six-channel BRD panel.**

After designing all six BRD targets independently (see Recipe 23 for the
five-target panel), add the M. bovis channel:

```bash
lamp-forge panel \
  --set BRSV=results/brsv_N_gene/primer_sets.json \
  --set BCOV=results/bcov_N_gene/primer_sets.json \
  --set BVDV=results/bvdv_5utr/primer_sets.json \
  --set IBR=results/ibr_gB/primer_sets.json \
  --set MHAE=results/mhae_lktA/primer_sets.json \
  --set MBOVIS=results/mbovis_uvrC/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/brd_6plex_panel
```

Generate the pooling sheet (six targets = 264 uM minimum; request 300 uM
resuspension from IDT or Twist):

```bash
lamp-forge pool \
  --panel results/brd_6plex_panel/panel.json \
  --stock-conc 300 \
  --total-volume 500 \
  --out results/brd_6plex_panel/pool_sheet.csv
```

Export all 36 primers in one IDT order sheet:

```bash
lamp-forge panel-export \
  --panel results/brd_6plex_panel/panel.json \
  --format idt \
  --out orders/brd_6plex_idt_order.csv
```

**Six-channel BRD result interpretation.**

| BRSV | BCoV | BVDV | IBR | MHAE | MBOVIS | Interpretation |
|---|---|---|---|---|---|---|
| + | - | - | - | + | - | Viral-bacterial co-infection; BRSV + M. haemolytica -- CRITICAL; initiate antimicrobials immediately |
| - | - | - | - | - | + | M. bovis alone; withhold standard BRD antimicrobials; susceptibility test before treating |
| + | - | - | - | + | + | BRSV + MHAE + MBOVIS; highest-complexity BRD; dual antimicrobial therapy required |
| - | - | + | - | - | - | BVDV alone; screen herd for PI (persistently infected) animals; movement restriction |
| - | - | - | - | - | - | All negative; consider extraction failure (include 16S control -- Recipe 19) |

A *M. bovis*-positive result in the absence of *M. haemolytica* is a red flag
for a pan-resistant strain.  Do not apply standard metaphylaxis or treatment
protocols without susceptibility data.

**Wet-lab notes.**

- M. bovis is fastidious and slow-growing; culture sensitivity is only
  ~40-60%.  LAMP sensitivity is substantially higher for direct clinical
  samples.
- For nasal swabs, collect both nasal swabs and a deep tracheal aspirate
  where possible -- M. bovis load in nasal secretions can be lower than
  in bronchoalveolar lavage in chronic pneumonia cases.
- For joint fluid, the DNA extraction requires an additional bead-beating
  step to lyse the organisms in the synovial matrix; standard boiling
  protocols have reduced efficiency for this matrix.
- Include a 16S rRNA internal control (Recipe 19) to confirm extraction
  success; M. bovis-negative / 16S-negative results indicate inhibition
  or extraction failure, not a clean sample.
- A positive result should be followed by culture and susceptibility testing
  at an accredited laboratory before initiating or modifying antimicrobial
  therapy.

**References.**

- Chen SQ, Zhang W, Chen ZL et al. (2017). Development of a loop-mediated
  isothermal amplification assay for rapid detection of *Mycoplasma bovis*.
  *J Vet Sci* 18(3):397-403. doi:10.4142/jvs.2017.18.3.397
- Hobbs JK-R, Bowen JM, Lord CC & Angen O (2018). Whole-genome sequencing
  and comparative genomics of *Mycoplasma bovis*. *Vet Microbiol* 220:18-25.
  doi:10.1016/j.vetmic.2018.04.029
- Caswell JL & Archambault M (2007). Mycoplasma bovis pneumonia in cattle.
  *Anim Health Res Rev* 8:161-186. doi:10.1017/S1466252307001351
- Gautier-Bouchardon AV (2018). Antimicrobial resistance in Mycoplasma spp.
  *Microbiol Spectr* 6(4):ARBA-0030-2018. doi:10.1128/microbiolspec.ARBA-0030-2018
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 25 -- Souring trajectory analysis for time-series oilfield MIC monitoring (`lamp-forge souring-trend`)

**Goal.** Convert a sequence of monthly oilfield MIC LAMP panel results
(from `lamp-forge mic-risk`) into a **souring trajectory assessment** --
IMPROVING, STABLE_LOW, STABLE_HIGH, or DETERIORATING -- together with key
inflection-point flags (SRB newly detected, SRB clearing) and a treatment
effectiveness audit.

**Why it's interesting.** A single LAMP panel reading tells you *what is
present now*.  The trend across consecutive readings tells you *whether the
situation is getting better or worse*, which is the decision variable that
drives treatment escalation, biocide dose adjustment, and regulatory
reporting.  A portable BioVind BioID deployed on a produced-water manifold
generates one result per monitoring interval; without trend context, each
result is interpreted in isolation, missing the forest for the trees.

Key operational decisions that depend on trend:

- **First SRB detection after consecutive negatives** (newly detected):
  requires immediate biocide treatment -- waiting for the next interval
  risks the SRB community becoming established.
- **SRB + NRB co-detected repeatedly** (treatment under-dose): nitrate
  injection is active but not achieving thermodynamic suppression of
  sulfate-reduction; dose must be increased.
- **SRB returning to negative** (clearing): confirms treatment effectiveness
  and justifies maintaining the current programme rather than escalating.
- **Chronic SRB positivity without a downward trend** (STABLE_HIGH): current
  programme is failing; a treatment review is overdue.

**Workflow -- from individual panel runs to trend report.**

Step 1: Run the LAMP panel at each monitoring interval and save the JSON.

```bash
# Month 1 -- SRB detected in produced water from well W1-A
lamp-forge mic-risk --srb \
  --out-json results/W1-A/2026-01.json

# Month 2 -- SRB + IRB co-detected
lamp-forge mic-risk --srb --irb \
  --out-json results/W1-A/2026-02.json

# Month 3 -- biocide applied; SRB still present with IRB
lamp-forge mic-risk --srb --irb \
  --out-json results/W1-A/2026-03.json

# Month 4 -- dose increased; NRB now detectable (nitrate injection active)
lamp-forge mic-risk --srb --irb --nrb \
  --out-json results/W1-A/2026-04.json

# Month 5 -- SRB cleared; NRB dominant
lamp-forge mic-risk --nrb \
  --out-json results/W1-A/2026-05.json
```

Step 2: Analyse the five-month trend.

```bash
lamp-forge souring-trend \
  --mic-result results/W1-A/2026-01.json \
  --mic-result results/W1-A/2026-02.json \
  --mic-result results/W1-A/2026-03.json \
  --mic-result results/W1-A/2026-04.json \
  --mic-result results/W1-A/2026-05.json \
  --out-json  results/W1-A/trend_2026_H1.json \
  --out-csv   results/W1-A/trend_2026_H1.csv
```

**Workflow -- from a monitoring spreadsheet CSV.**

If results are tracked in a spreadsheet, export a CSV with the required columns:

```
# monitoring/W1-A_2026.csv
sample_id,date,srb,mcr,irb,apb,nrb
2026-01,2026-01-15,1,0,0,0,0
2026-02,2026-02-15,1,0,1,0,0
2026-03,2026-03-15,1,0,1,0,0
2026-04,2026-04-15,1,0,1,0,1
2026-05,2026-05-15,0,0,0,0,1
```

```bash
lamp-forge souring-trend \
  --csv monitoring/W1-A_2026.csv \
  --out-json results/W1-A/trend_2026_H1.json
```

Boolean columns accept `0`/`1` or `true`/`false` (case-insensitive).
The `date` column is optional.

**Trend direction logic.**

| Pattern | Direction |
|---|---|
| SRB positive only in the most recent sample | DETERIORATING (onset event) |
| SRB rate rising between first half and second half of the time series | DETERIORATING |
| SRB positive in first sample, negative in most recent sample | IMPROVING (clearing event) |
| SRB rate falling, most recent sample negative | IMPROVING |
| SRB detected in >= 50% of intervals, no clear downward trend | STABLE_HIGH |
| SRB detected in < 50% of intervals, no clear trend | STABLE_LOW |
| Fewer than 2 samples available | INSUFFICIENT_DATA |

**Treatment effectiveness flags.**

| Flag | Condition | Operational meaning |
|---|---|---|
| `srb_newly_detected` | SRB positive only in the last sample | Souring onset; immediate biocide required |
| `srb_clearing` | SRB positive in first sample, negative in last | Treatment effective; maintain current programme |
| `treatment_underdosed_count > 0` | SRB + NRB both positive in N intervals | Nitrate injection present but under-dosed |
| NRB positive, SRB consistently negative | NRB detectable, no SRB | Nitrate suppression effective; continue programme |

**Multi-well campaign.**

Run `souring-trend` once per well and aggregate the direction flags in a management
dashboard.  Wells showing DETERIORATING should be prioritised for immediate
treatment; STABLE_HIGH wells require programme review; IMPROVING and STABLE_LOW
wells can remain on the standard monitoring schedule.

```bash
for well in W1-A W1-B W2-A W2-B; do
  lamp-forge souring-trend \
    --csv monitoring/${well}_2026.csv \
    --out-json results/${well}/trend_2026.json
done
```

**References.**

- Hubert C & Voordouw G (2007). Oil field souring control by nitrate-reducing
  bacteria. *Appl Environ Microbiol* 73:2644-2652.
  doi:10.1128/AEM.02751-06
- Hamilton WA (2003). Microbially influenced corrosion as a model system for
  the study of metal microbe interactions. *Int Biodeterior Biodegrad*
  51(3):151-155. doi:10.1016/S0964-8305(02)00089-6
- Vigneron A et al. (2017). Comammox Nitrospira in oilfield MIC microbiomes.
  *ISME J* 11:2180-2194. doi:10.1038/ismej.2017.40
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 26 -- Surveillance sampling plan (`lamp-forge sampling-plan`)

**Goal.** Before a field campaign, determine *how many* wells, animals, or
patients to test per site visit so that you would detect the target pathogen
at a stated minimum prevalence with a stated confidence level.

**Why it matters for BioVind.** A portable BioID device makes per-test
cost low enough that systematic surveillance is practical -- but a poorly
designed sampling plan either wastes reactions (too many samples) or
misses outbreaks (too few).  `lamp-forge sampling-plan` closes this gap
by computing the minimum sample size *n* before a campaign starts, using
the same probabilistic model that regulatory agencies use for freedom-from-
disease surveys.

### Model

For an effectively infinite population the standard formula is:

```
n = ceil[ log(1 - confidence) / log(1 - prevalence * sensitivity) ]
```

where *prevalence* is the minimum fraction of positive units the campaign
must detect, *sensitivity* is the LAMP assay's analytical sensitivity
(probability of a positive result given target present), and *confidence*
is the required probability of finding at least one positive sample if the
true prevalence equals the design prevalence.

When the total population size is known (a specific herd, a pipeline with a
fixed number of monitored wells), supplying `--population-size` applies the
hypergeometric finite-population correction, which can substantially reduce
the required *n* for small populations.

### BioVind vertical examples

**Oil & gas -- SRB souring surveillance (oilfield wells).**
How many wells must be tested per monthly monitoring visit to detect SRB
presence if >= 5% of wells are positive?

```bash
lamp-forge sampling-plan \
  --prevalence 0.01 \
  --prevalence 0.02 \
  --prevalence 0.05 \
  --prevalence 0.10 \
  --prevalence 0.20 \
  --confidence 0.95 \
  --sensitivity 0.95 \
  --out-csv results/srb_sampling_plan.csv
```

Expected output (abridged):

```
Surveillance sampling plan: confidence=95%, sensitivity=95%, population=infinite

    Prevalence (%)   n_samples   Achieved conf (%)
---------------------------------------------------
              1.00         313               95.03
              2.00         157               95.01
              5.00          63               95.10
             10.00          31               95.09
             20.00          16               95.23
```

At the common produced-water management threshold of 5% positive wells,
63 samples per monitoring interval are sufficient at 95% confidence.
If your field has only 20 monitored wells in total, add `--population-size 20`
to apply the finite-population correction (which will reduce *n* to fewer
than the 20 total).

**Farm biosecurity -- ASFV detection in a swine herd.**
How many pigs to test to detect ASFV if >= 5% of a 200-pig herd are infected?

```bash
lamp-forge sampling-plan \
  --prevalence 0.05 \
  --prevalence 0.10 \
  --confidence 0.95 \
  --sensitivity 0.95 \
  --population-size 200 \
  --out-csv results/asfv_farm_sampling.csv
```

Expected output:

```
Surveillance sampling plan: confidence=95%, sensitivity=95%, population=200 units

    Prevalence (%)   n_samples   Achieved conf (%)  FPC
------------------------------------------------------
              5.00          53               95.17   yes
             10.00          25               95.01   yes
```

The "FPC yes" flag indicates the finite-population correction reduced the
sample size below the infinite-population estimate (59 -> 53 at 5%).  In
practice, testing 53 animals from a 200-pig herd per outbreak investigation
is operationally feasible with a portable device.

**Human POC -- MTB clinic surveillance.**
How many patients per week must be screened to detect a TB prevalence of
>= 1% among symptomatic respiratory patients visiting a rural clinic?

```bash
lamp-forge sampling-plan \
  --prevalence 0.01 \
  --prevalence 0.02 \
  --prevalence 0.05 \
  --confidence 0.95 \
  --sensitivity 0.95
```

At 1% prevalence, 314 patient screens per monitoring period are needed --
roughly 45 per day over a week.  At a clinically relevant threshold of 5%,
only 63 screens are required (< 10 per day).

### Combining with LOD

Sampling plan and LOD analysis address complementary questions:

1. **LOD** (`lamp-forge lod`): Given a single positive sample, what is the
   minimum concentration of target DNA/RNA the assay can reliably detect?
2. **Sampling plan** (`lamp-forge sampling-plan`): Given that at-risk units
   in the population will be above the assay LOD when positive, how many
   units must be tested to find the condition if it is present at prevalence p?

Use both before ordering primers and planning a field campaign:

```bash
# Step 1: confirm the assay is sensitive enough (LOD)
lamp-forge lod \
  --sample-volume 1000 \
  --efficiency 0.50 \
  --eluate-volume 100 \
  --reaction-input 5

# Step 2: plan how many samples to collect (sampling plan)
lamp-forge sampling-plan \
  --prevalence 0.05 \
  --confidence 0.95 \
  --sensitivity 0.95
```

### Script usage

The sampling planner is also importable for use in Jupyter notebooks or
decision-support tools:

```python
from lamp_forge.sampling import SamplingParams, sample_size, sampling_plan

# Oilfield SRB monitoring -- detect if >= 5% of wells are positive
params = SamplingParams(
    target_prevalence=0.05,
    confidence=0.95,
    test_sensitivity=0.95,
)
result = sample_size(params)
print(f"Test {result.n_samples} wells per visit "
      f"(achieved confidence: {result.achieved_confidence * 100:.1f}%)")

# Table across multiple prevalences
table = sampling_plan(
    (0.01, 0.02, 0.05, 0.10, 0.20),
    confidence=0.95,
    sensitivity=0.95,
)
for row in table:
    print(f"  p={row.target_prevalence*100:.0f}%: n={row.n_samples}")
```

**References.**

- Cannon RM & Roe RT (1982). Livestock Disease Surveys: A Field Manual
  for Veterinarians. Australian Bureau of Animal Health, Canberra.
- Cameron AR & Baldock FC (1998). A new probability formula for surveys to
  substantiate freedom from disease. *Prev Vet Med* 34:1-17.
  doi:10.1016/S0167-5877(97)00081-0
- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 27 -- Porcine Epidemic Diarrhea Virus (PEDV) via RT-LAMP (on-farm neonatal biosecurity)

**Goal.** Detect Porcine Epidemic Diarrhea Virus (PEDV) from intestinal content,
rectal swab, or faecal samples in farrowing barns -- fast enough to trigger immediate
biosecurity interventions before the virus spreads to adjacent litters.

**Why it's interesting.** PEDV devastated the US pig industry in 2013-14, killing
an estimated >8 million piglets and eliminating >10% of the US sow herd in its
first year.  Near-100% case fatality in neonatal piglets (<1 week old) makes it
the highest-mortality enteric disease in swine production.  The virus spreads
extremely rapidly through faecal-oral contact and aerosol (barn-to-barn spread
documented).  A portable LAMP result at the farm in under 60 minutes -- before
a potentially PEDV-positive transport vehicle leaves -- is exactly the use case
the BioVind BioID platform targets.

PEDV is a **positive-sense ssRNA virus** (*Alphacoronavirus 1*, the same genus as
the common-cold coronavirus 229E) -- so detection requires **RT-LAMP**.  Two major
genotype clusters circulate: **G1** (classical strains including the European CV777
reference and the attenuated Korean DR13 vaccine strain) and **G2** (virulent strains
including USA-Ohio851-2013 and the AH2012 and CH/FJND-3/2011 Asian field strains that
drove global spread).  G1 and G2 share ~96-99% N-gene nucleotide identity, so a
conserved N-gene LAMP assay detects both without genotype discrimination -- which is
the correct scope for a first-pass detection assay.  Genotype-level characterisation
can follow at a reference laboratory if needed.

**Target.** The **N gene** (nucleocapsid protein, ~1326 nt CDS) -- the most conserved
coding region across all PEDV genotypes and the basis of published PEDV RT-PCR and
RT-LAMP diagnostic assays (He et al. 2019; Kim et al. 2014).  Unlike the spike (S)
gene, which diverges >5% nt between G1 and G2, the N gene sits in a part of the
coronavirus genome under tight assembly constraint, yielding multiple conserved
blocks usable for LAMP primer placement.

Use the ready-made config at `config/pedv_N_gene.yaml`, or paste the block below:

```yaml
target:
  name: pedv_N_gene
  taxon_id: 27317          # Porcine epidemic diarrhea virus
  gene: N                  # nucleocapsid protein; alt: "nucleoprotein" / "N protein"
  max_sequences: 30        # span G1a/G1b classical strains + G2a/G2b virulent variants
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- TGEV shares alphacoronavirus N-gene family
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.30        # RNA-virus drift: looser than bacteria
  min_region_length: 200

primer:
  tm_min: 63.0                   # RT-LAMP one-step: co-activity floor for NEB RTx / Bst 2.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0                   # PEDV N gene is ~38-42% GC; wide margins for variant diversity
  gc_max: 55.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/pedv_N_gene
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `primer.tm_min: 63.0` -- one-step RT-LAMP at 63-65 degC; same rationale as
  PRRSV (Recipe 11), AIV (Recipe 12), FMDV (Recipe 14), and NDV (Recipe 18).
- `conservation.entropy_threshold: 0.30` -- slightly looser than bacterial
  housekeeping genes (0.20 bits) to tolerate synonymous drift between G1 and G2
  while still resolving the tight N-gene core blocks.
- `min_identity_threshold: 0.85` -- tighter than the default 0.80 because TGEV
  shares the alphacoronavirus N-gene fold; 0.85 catches potential cross-reactivity.

**Off-target panel.** The critical differentials are other porcine enteric
coronaviruses that co-circulate in farrowing barns and cause indistinguishable
diarrhoea in neonatal piglets:

| File | Source | Why |
|---|---|---|
| `tgev.fasta` | Transmissible gastroenteritis virus (NC_038861.1 or D00563.1 region) | Closest porcine coronavirus (Alphacoronavirus 1 -- same species as PEDV); shares N-gene structural motifs |
| `pdcov.fasta` | Porcine deltacoronavirus (NC_029977.1 or KU984334.1) | Emerging enteric coronavirus; co-circulates and causes identical clinical signs in neonates |
| `porcine_rotavirus_a.fasta` | Porcine rotavirus A (NC_011509.1, VP4 segment) | Most common cause of neonatal piglet diarrhoea; co-infects with PEDV during outbreaks |
| `sus_scrofa_fragment.fasta` | Subset of *Sus scrofa* genome (Sscrofa11.1) | Host DNA from intestinal content or rectal swab |

**Expected output.** The N gene (~1326 nt) typically yields 2-4 conserved windows
usable for LAMP.  Top-scoring sets should have zero flagged off-target hits across
the panel -- PEDV is genetically distant from rotavirus (double-stranded RNA
Reoviridae), and the identity thresholds should clearly separate it from TGEV and
PDCoV.  If TGEV identity hits appear, raise `min_identity_threshold` to 0.90 and
prefer primer sets whose landing zones overlap the highly divergent hypervariable
loop of the PEDV N-gene rather than the shared N-terminal RNA-binding domain.

**Verify RT-LAMP readiness.**

```bash
lamp-forge rt-check \
  --input results/pedv_N_gene/primer_sets.json \
  --na-type rna \
  --out-csv results/pedv_N_gene/rt_check.csv
```

Sets marked **NOT OPTIMIZED** have at least one core primer below 63.0 degC.
Because `primer.tm_min: 63.0` is already set in this config, all sets should
pass -- but verify before ordering.

**Estimate LOD before ordering.** Intestinal content from affected piglets is a
moderately inhibitory matrix.  A standard protocol (0.5 g intestinal content,
RNA extraction into 100 uL eluate, 5 uL to a 25 uL RT-LAMP reaction) gives:

```bash
lamp-forge lod \
  --sample-volume 500 \
  --efficiency 0.40 \
  --eluate-volume 100 \
  --reaction-input 5 \
  --out-csv results/pedv_N_gene/lod_intestinal.csv
```

Effective sample = 500 x 0.4 x (5/100) = 10 uL per reaction.  LOD_95 approx
300 copies/mL of intestinal homogenate.  PEDV titres in acutely infected
piglet intestines can exceed 10^7 genome equivalents/mL -- orders-of-magnitude
above this LOD -- so even a simple extraction protocol is sufficient for
clinical detection.

**Confirm TTP within the BioVind device window.**

```bash
lamp-forge ttp \
  --preset rt-lamp \
  --device-window 60
```

At LOD_95 (~3 copies/reaction), TTP is expected around 55-60 min for a
one-step RT-LAMP protocol; the device window is tight.  Titres encountered
during acute clinical infection (10^3-10^7 copies/reaction) will read out
at 30-45 min, well within the 60-min window.

**Multiplex with other farm-biosecurity targets.** PEDV outbreaks frequently
co-occur with PRRSV (immunosuppression increases susceptibility) and can be
confused with porcine transmissible gastroenteritis on clinical grounds.
Design each target independently, then combine:

```bash
lamp-forge panel \
  --set PEDV=results/pedv_N_gene/primer_sets.json \
  --set PRRSV=results/prrsv_ORF7/primer_sets.json \
  --set ASFV=results/asfv_p72/primer_sets.json \
  --top-per-target 5 \
  --out results/swine_biosecurity_3plex
```

Then calculate the pooling sheet:

```bash
lamp-forge pool \
  --panel results/swine_biosecurity_3plex/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/swine_biosecurity_3plex/pool_sheet.csv
```

And export the PEDV primers to IDT for ordering:

```bash
lamp-forge export \
  --input results/pedv_N_gene/primer_sets.json \
  --format idt \
  --target-label PEDV_N \
  --out orders/pedv_N_idt_order.csv
```

**Field note.** PEDV spreads via "biosecurity-defying" fomite transmission --
contaminated boots, vehicles, and feed-delivery equipment.  A positive PEDV
result on a batch sample from a recently-delivered load of weaned piglets
should trigger immediate pen quarantine, deep cleaning, and notification of
nearby farms sourcing from the same supplier.  Position this assay as a
**triage/detection** tool; PEDV is not a statutory notifiable disease in most
jurisdictions (unlike ASFV or FMDV), but many national swine health programmes
recommend immediate herd reporting.  Validate against a panel of known
positive and negative intestinal samples before deploying as a routine
farm assay.

**References.**

- He X et al. (2019) Development and evaluation of a RT-LAMP assay for rapid
  detection of porcine epidemic diarrhea virus.  *Transbound Emerg Dis*
  66(1):431-439. doi:10.1111/tbed.13044
- Kim SY et al. (2014) Loop-mediated isothermal amplification assay for
  rapid, sensitive detection of porcine epidemic diarrhea virus.
  *J Virol Methods* 208:92-96. doi:10.1016/j.jviromet.2014.08.001
- Stevenson GW et al. (2013) Emergence of porcine epidemic diarrhea virus
  in the United States.  *J Vet Diagn Invest* 25(5):649-654.
  doi:10.1177/1040638713501675
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63


## Recipe 28 -- One-shot BioID cartridge build (`lamp-forge cartridge`)

Spinning up a new BioVind BioID cartridge SKU means running seven LAMP-Forge
commands in sequence -- `panel`, `pool`, `lod`, `rt-check`, `preorder`,
`panel-export` -- and then reconciling the outputs by hand into a single
GO / WARNING / NO_GO call that decides whether the panel goes to vendor.

Recipe 28 collapses that workflow into one command driven by a single
**cartridge manifest YAML** that names every channel, the shared extraction
chain, the device-window budget, the pooling parameters, and the vendor
format.  The output is a signed, content-addressed cartridge build that
your CI can gate on.

### Cartridge manifest

`config/cartridges/cartridge_respiratory_18ch.yaml` describes an 18-channel respiratory
panel (17 detect channels + a 16S rRNA sample-adequacy control):

```yaml
cartridge_id: resp_18ch_v1
chemistry: rt-lamp
max_channels: 18

extraction:
  sample_volume_ul: 1000.0
  extraction_efficiency: 0.50
  eluate_volume_ul: 100.0
  reaction_input_ul: 5.0

device_window_min: 60.0

pool:
  stock_conc_um: 500.0
  total_pool_volume_ul: 500.0

vendor: idt

targets:
  - label: SARS2_N
    sets_json: ../results/sars_cov_2_N/primer_sets.json
    target_na: rna
    role: detect
  - label: IAV_M
    sets_json: ../results/influenza_a_M/primer_sets.json
    target_na: rna
    role: detect
  # ...
  - label: CTRL_16S
    sets_json: ../results/general_16s/primer_sets.json
    target_na: dna
    role: control
```

Each `sets_json` path points at a `primer_sets.json` produced by a prior
per-target `lamp-forge run`.  Paths in the YAML are resolved relative to the
manifest file's directory, so a manifest in `config/` referring to
`../results/...` works regardless of the user's CWD.

### One-shot build

```bash
lamp-forge cartridge config/cartridges/cartridge_respiratory_18ch.yaml \
  --out-dir results/resp_18ch_v1
```

The command writes eight artifacts into `--out-dir`:

| Artifact | Producer | What it is |
|---|---|---|
| `cartridge_manifest.json` | this module | Signed top-level build record (see below). |
| `panel.json`, `panel_primers.csv`, `panel_report.html` | `panel.run_panel` | Cross-dimer-clean per-channel set selection. |
| `pool_plan.csv` | `pool.build_pool_plan` | Pipetting sheet for the multiplex working stock. |
| `lod_panel.csv` | local | One LoD row per channel + a `WORST(...)` summary row. |
| `rt_check.csv` | `rt_lamp.write_rt_check_csv` | Stacked RT-LAMP Tm check across every channel. |
| `preorder_<label>.csv` | `preorder.write_preorder_csv` | One per channel. |
| `vendor_order.csv` | `vendor_export.write_vendor_csv` | Single IDT/Twist order CSV across every channel. |
| `build_report.html` | local | Human-readable cartridge summary + LoD bar chart. |

### Cartridge-level status rollup

`cartridge_manifest.json` carries a single `build_status` field:

```jsonc
{
  "schema_version": "1.0",
  "cartridge_id": "resp_18ch_v1",
  "chemistry": "rt-lamp",
  "max_channels": 18,
  "n_targets": 18,
  "build_status": "GO",
  "reasons": [
    "All 18 channel(s) pass preorder, panel is clean, and every RNA channel is RT-LAMP optimised."
  ],
  "per_target": [
    {
      "label": "SARS2_N",
      "set_id": "region_03_set_017",
      "composite_score": 0.873,
      "target_na": "rna",
      "preorder_status": "GO",
      "rt_check": { "n_sets": 5, "n_sets_rt_ok": 5 },
      "lod": { "lod_copies_per_rxn": 2.996, "lod_copies_per_ml": 119.84 }
    }
    // ... 17 more
  ],
  "panel_summary": {
    "worst_inter_assay_dg": -3.42,
    "is_clean": true,
    "flagged_cross_dimers": []
  },
  "panel_lod_worst": { "label": "IBR_gB", "lod_copies_per_ml": 119.84 },
  "manifest_hash": "8f2c... (SHA-256 of manifest + per-target outputs)"
}
```

Rollup rules:

- Any channel with `preorder_status == "NO_GO"` -> cartridge **NO_GO**.
- Panel `is_clean == false` (flagged inter-assay cross-dimers) -> cartridge **NO_GO**.
- Any RNA channel with zero RT-optimised sets -> cartridge **NO_GO**.
- Any channel `WARNING` (TTP near limit, some RT sets sub-optimal) -> cartridge **WARNING**.
- Otherwise -> cartridge **GO**.

`manifest_hash` is the SHA-256 of the canonical JSON over the inputs
(cartridge id, chemistry, extraction, pool, vendor, every target's
`primer_sets.json` SHA-256) and the per-target outputs.  Two back-to-back
runs over the same inputs produce the same hash; changing a single primer
base in any `primer_sets.json` changes the hash.

### Wiring into a pre-order CI gate

`lamp-forge cartridge` returns a non-zero exit code when `build_status ==
"NO_GO"`.  Drop the command into the same CI step that uploads to your
vendor portal:

```yaml
# .github/workflows/preorder.yml (excerpt)
- name: Build cartridge
  run: |
    lamp-forge cartridge config/cartridges/cartridge_respiratory_18ch.yaml \
      --out-dir results/resp_18ch_v1
- name: Upload vendor CSV
  if: success()
  run: ./scripts/upload_to_idt.sh results/resp_18ch_v1/vendor_order.csv
- name: Archive signed manifest
  if: always()
  uses: actions/upload-artifact@v4
  with:
    name: cartridge-manifest
    path: results/resp_18ch_v1/cartridge_manifest.json
```

A dirty panel or a failed RT check now physically blocks the order
upload -- and the archived `cartridge_manifest.json` is the audit trail
for every cartridge SKU that does ship.

**Field note.** The cartridge command never re-designs primers; it composes
the existing pipeline's outputs.  When a build comes back **NO_GO** for
panel uncleanness, the fix is upstream: re-run `lamp-forge run` for the
offending channel(s) in a different conserved window, or tighten the Tm
window for sub-optimal RT-LAMP sets, then re-build.  The signed
`manifest_hash` makes it obvious whether a re-build is over the same
underlying primer sets or a different set.

**References.**

- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63
- Tanner NA et al. (2012) Visual detection of isothermal nucleic acid
  amplification using pH-sensitive dyes.  *BioTechniques* 53(2):81-89.
  doi:10.2144/0000113902
- New England Biolabs.  WarmStart(R) LAMP Kit (DNA & RNA) instruction
  manual, E1700.  2024.

---

## Recipe 29 -- Classical Swine Fever Virus (CSFV) via RT-LAMP (on-farm / checkpoint biosecurity)

**Goal.** Detect Classical Swine Fever Virus (CSFV; also called hog cholera or
European swine fever) broadly across all three genotype groups (1, 2, 3) from
whole blood, serum, or tonsil swab on the farm or at a border checkpoint -- fast
enough to support a quarantine decision before pigs move.

**Why it's interesting.** CSFV is a WOAH-listed Tier 1 notifiable disease with
near-100% case fatality in naive herds.  The disease is eradicated from North
America, the EU, and most of Oceania, but it remains endemic in parts of Asia,
Central America, and sub-Saharan Africa.  Critically, CSFV produces a clinical
picture that is **nearly identical to African Swine Fever Virus (ASFV)**: acute
fever, hemorrhagic diathesis, and rapid spread.  A BioVind portable device that
runs CSFV (Recipe 29) and ASFV (Recipe 7) in the same panel resolves this
differential on-farm in under 60 minutes -- the fastest route to a quarantine
or trade-restriction decision when both diseases are differential diagnoses.

CSFV has three genotype groups (1, 2, 3) with nine subtypes.  Subtype 2.1 accounts
for the majority of contemporary outbreaks in East and Southeast Asia.  A pan-CSFV
assay must capture all three genotype groups, so the primer design needs a diverse,
genotype-spanning input set.

CSFV is a **positive-sense ssRNA** virus (*Flaviviridae*, Genus *Pestivirus A*),
so detection requires **RT-LAMP**.  One-step RT-LAMP (Bst 2.0 WarmStart + NEB RTx)
at 63-65 degC delivers a yes/no result in under 30 minutes, compatible with the
BioVind-style portable platform.

**Target.** `NS5B` (RNA-directed RNA polymerase, ~1.8 kb) -- the most conserved
ORF across all CSFV genotypes and the target used in published CSFV RT-PCR and
RT-LAMP assays (Wang et al. 2018; Rodriguez-Prieto et al. 2017).  The alternative
E2 glycoprotein gene is more subtype-variable and better suited to
subtype-discriminating assays outside the scope of a first-pass pan-CSFV screen.

Use the ready-made config at `config/csfv_NS5B.yaml`, or paste the block below:

```yaml
target:
  name: csfv_NS5B
  taxon_id: 11096          # Classical swine fever virus (Pestivirus A)
  gene: NS5B               # RNA polymerase; alt annotation: "NS5b" / "p58" / "RdRp"
  max_sequences: 35        # span genotype groups 1/2/3 + major subtypes
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- BVDV shares pestivirus NS5B family
  min_coverage_threshold: 0.85

conservation:
  window_size: 25
  entropy_threshold: 0.30        # RdRp catalytic constraint; looser than housekeeping genes
  min_region_length: 200

primer:
  tm_min: 63.0                   # RT-LAMP one-step: co-activity floor for NEB RTx / Bst 2.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 43.0                   # CSFV NS5B ~51-55% GC; margins for subtype diversity
  gc_max: 62.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/csfv_NS5B
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config note.** `primer.tm_min: 63.0` -- one-step RT-LAMP at 63-65 degC;
same rationale as Recipes 11 (PRRSV), 12 (AIV), 14 (FMDV), and 18 (NDV).
`conservation.entropy_threshold: 0.30` is slightly looser than for bacterial
housekeeping genes (0.20 bits) -- RNA-virus synonymous drift is higher.

**SPECIFICITY CHALLENGE: BVDV.** CSFV is a pestivirus closely related to Bovine
Viral Diarrhea Virus (BVDV, Pestivirus B), which is endemic in cattle worldwide
and can be present on mixed cattle/pig farms.  The CSFV and BVDV NS5B genes share
~50-60% nucleotide identity -- high enough that a loosely designed assay will
cross-react.  `min_identity_threshold: 0.85` and inclusion of BVDV in the off-target
panel are **essential** for any CSFV assay intended for field deployment on
mixed-species farms.

**Off-target panel.** The key differentials are other pestiviruses and co-circulating
swine pathogens:

| File | Source | Why |
|---|---|---|
| `bvdv.fasta` | NC_001461.1 (BVDV-1) + NC_002539.1 (BVDV-2) | **Most important off-target** -- closest pestivirus relative; endemic in cattle on mixed farms |
| `bdv.fasta` | NC_003679.1 | Border disease virus (Pestivirus D); sheep pestivirus in the same genus |
| `prrsv.fasta` | PRRS virus *(representative)* | Co-circulates in pigs; RNA arterivirus |
| `asfv_p72.fasta` | NC_001659.1 region | ASFV (DNA Asfarviridae); primary clinical differential -- a DNA virus; must not cross-flag |
| `sus_scrofa_fragment.fasta` | Subset of *Sus scrofa* genome | Host DNA from pig blood or tissue |

**Expected output.** NS5B yields several broadly-conserved RdRp-motif windows.
Top sets should have zero flagged hits across the panel, including BVDV -- if BVDV
hits appear, raise `min_identity_threshold` to 0.90 and restrict to windows in the
most conserved catalytic core (GDD motif region).  Because ASFV is a large DNA virus
genetically unrelated to CSFV, there should be no cross-reactivity with ASFV.

**Verify RT-LAMP readiness.**

```bash
lamp-forge rt-check \
  --input results/csfv_NS5B/primer_sets.json \
  --na-type rna \
  --out-csv results/csfv_NS5B/rt_check.csv
```

Sets marked **NOT OPTIMIZED** have at least one core primer below 63.0 degC.
Because `primer.tm_min: 63.0` is already set in this config, suboptimal sets
are unusual -- but verify before ordering.

**Estimate LOD before ordering** (whole blood 200 uL, 50% RNA extraction,
50 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 200 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 200 x 0.5 x (5/50) = 10 uL per reaction -> LOD_95 approx
300 copies/mL of whole blood.  CSFV viral loads in blood during the acute febrile
phase can reach 10^6-10^8 TCID50/mL, providing orders-of-magnitude headroom for
detection before clinical signs peak.

**Multiplex CSFV with other swine biosecurity targets** to resolve the
ASFV/CSFV clinical overlap and complete a swine differential panel:

```bash
lamp-forge panel \
  --set CSFV=results/csfv_NS5B/primer_sets.json \
  --set ASFV=results/asfv_p72/primer_sets.json \
  --set PRRSV=results/prrsv_ORF7/primer_sets.json \
  --top-per-target 5 \
  --out results/swine_biosecurity_3plex
```

Then generate the pooling sheet and export for ordering:

```bash
lamp-forge pool \
  --panel results/swine_biosecurity_3plex/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/swine_biosecurity_3plex/pool_sheet.csv

lamp-forge export \
  --input results/csfv_NS5B/primer_sets.json \
  --format idt \
  --target-label CSFV_NS5B \
  --out orders/csfv_NS5B_idt_order.csv
```

**Field note.** CSFV is a WOAH-listed Tier 1 notifiable disease.  A positive
result must trigger **immediate notification to the national veterinary authority**
and confirmatory testing at an approved reference laboratory before any depopulation
or trade-restriction decisions.  Position this assay as a **triage/screening** tool.
Wet-lab validation must be performed at an approved facility using inactivated
positive controls or synthetic RNA standards.

**References.**

- Wang A, Cai X, Chen H et al. (2018) Rapid and sensitive detection of classical
  swine fever virus by reverse transcription loop-mediated isothermal amplification.
  *Vet Microbiol* 216:66-71. doi:10.1016/j.vetmic.2018.02.003
- Rodriguez-Prieto V, Vicente-Rubiano M, Sanchez-Matamoros A et al. (2017)
  Real-time and conventional RT-LAMP for detection of classical swine fever virus.
  *J Virol Methods* 250:7-12. doi:10.1016/j.jviromet.2017.09.014
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 30 -- Poultry biosecurity four-target panel (`lamp-forge poultry-risk`)

**Goal.**  A single-cartridge, on-farm LAMP panel that simultaneously screens
commercial broiler, layer, and breeder flocks for the four most economically
important and / or WOAH-notifiable poultry diseases: Avian Influenza A (AIV),
Newcastle Disease (NDV), Infectious Bursal Disease (IBDV / Gumboro), and
Infectious Bronchitis (IBV).

**Why it complements the existing panels.** The farm biosecurity panel
(Recipe 8 / `lamp-forge farm-risk`) covers a mixed-livestock view (ASFV, FMDV,
AIV, NDV, PRRSV) optimised for integrated swine-and-poultry operations.
Dedicated poultry farms need a tighter panel that includes IBDV and IBV --
the two non-WOAH-listed but commercially devastating diseases -- and omits
swine-specific targets.  This recipe builds that dedicated panel.

**Pathogens and target genes.**

| Target | Gene   | Assay type | WOAH status | Weight |
|--------|--------|-----------|-------------|--------|
| AIV    | M gene | RT-LAMP   | Tier 1 (HPAI) | 60 |
| NDV    | M gene | RT-LAMP   | Listed        | 35 |
| IBDV   | VP2    | LAMP      | Listed (some regions) | 20 |
| IBV    | N gene | RT-LAMP   | Not globally listed | 10 |

Both AIV and NDV are RNA viruses requiring RT-LAMP (RTx + Bst 2.0).  IBDV is
a double-stranded RNA virus; the VP2 assay also uses a one-step RT-LAMP
protocol.  IBV (gamma coronavirus) is single-stranded positive-sense RNA and
likewise runs as RT-LAMP.  All four can be designed with
`--na-type rna` / `--vertical farm` in the scaffold command.

**Designing the four assays (one command per target).**

```bash
# AIV -- all subtypes via the M gene (pan-IAV)
lamp-forge scaffold \
  --target-name influenza_a_M \
  --vertical farm \
  --na-type rna \
  --taxon-id 11520 \
  --gene M \
  --out config/influenza_a_M.yaml

# NDV -- pan-NDV via the M gene
lamp-forge scaffold \
  --target-name ndv_M_gene \
  --vertical farm \
  --na-type rna \
  --taxon-id 11234 \
  --gene M \
  --out config/ndv_M_gene.yaml

# IBDV -- VP2 (most conserved; covers classic and vvIBDV)
lamp-forge scaffold \
  --target-name ibdv_VP2 \
  --vertical farm \
  --na-type rna \
  --taxon-id 12227 \
  --gene VP2 \
  --out config/ibdv_VP2.yaml

# IBV -- N gene (pan-IBV; serotype agnostic)
lamp-forge scaffold \
  --target-name ibv_N_gene \
  --vertical farm \
  --na-type rna \
  --taxon-id 11120 \
  --gene N \
  --out config/ibv_N_gene.yaml
```

Run each design:

```bash
for cfg in influenza_a_M ndv_M_gene ibdv_VP2 ibv_N_gene; do
  lamp-forge run --config config/${cfg}.yaml
done
```

**Off-target panel.** A poultry panel must not cross-react with:

| File | Source | Why |
|---|---|---|
| `gallus_gallus_fragment.fasta` | GRCg7b chr1 subset | Host (*Gallus gallus*) DNA in all swab samples |
| `ibv_Ma5.fasta` | NC_001451.1 | Vaccine-strain IBV (Massachusetts 5); must not be misidentified as field strain |
| `ndv_lentogenic.fasta` | Representative La Sota strain | Vaccine NDV; distinguish from velogenic field strains |
| `ibdv_classic.fasta` | Classic IBDV reference | Verify vvIBDV assay also captures classic strains |
| `mycoplasma_gallisepticum.fasta` | Representative MG strain | Common co-infection; must not produce false positive |

**RT-LAMP readiness check.** Verify all four assay sets meet the 63 degC
minimum incubation-temperature requirement for one-step RT-LAMP:

```bash
for target in influenza_a_M ndv_M_gene ibdv_VP2 ibv_N_gene; do
  echo "=== ${target} ==="
  lamp-forge rt-check \
    --input results/${target}/primer_sets.json \
    --na-type rna \
    --out-csv results/${target}/rt_check.csv
done
```

Sets flagged **NOT OPTIMIZED** have at least one core primer below 63 degC;
tighten `primer.tm_min: 63.0` in the corresponding config and re-run design.

**Multiplex cross-dimer screen.** Once all four assays have clean primer
sets, screen for inter-assay cross-dimerisation (ΔG < -5 kcal/mol is
flagged as incompatible):

```bash
lamp-forge panel \
  --set AIV=results/influenza_a_M/primer_sets.json \
  --set NDV=results/ndv_M_gene/primer_sets.json \
  --set IBDV=results/ibdv_VP2/primer_sets.json \
  --set IBV=results/ibv_N_gene/primer_sets.json \
  --top-per-target 5 \
  --out results/poultry_4plex
```

The panel command emits `panel.json` (compatible set) and `panel_report.html`
(cross-dimer heatmap).  The BioVind BioID cartridge has up to 18 channels,
so a four-target poultry panel fits comfortably alongside a 16S extraction
control, leaving headroom for future targets such as Marek's Disease Virus
(MDV) or Mycoplasma gallisepticum (MG).

**Equimolar pooling sheet.**

```bash
lamp-forge pool \
  --panel results/poultry_4plex/panel.json \
  --stock-conc 100 \
  --total-volume 500 \
  --out results/poultry_4plex/pool_sheet.csv
```

For four targets the combined pool concentration is 4 x 44 uM = 176 uM, which
exceeds a standard 100 uM stock.  Request 200 uM resuspension from IDT or
Twist, or order pre-pooled plates with custom concentrations:

```bash
lamp-forge export \
  --input results/influenza_a_M/primer_sets.json \
  --format idt \
  --target-label AIV_M \
  --out orders/AIV_M_idt.csv
# ... repeat for each target, then concatenate the CSV files
```

**Interpreting a result with `lamp-forge poultry-risk`.**

After the BioID device reads the cartridge, enter the positivity flags
directly to get a structured alert:

```bash
# All four negative (clean flock)
lamp-forge poultry-risk

# IBDV detected in broiler chicks (Gumboro outbreak)
lamp-forge poultry-risk --ibdv

# AIV + IBV co-detected (presumptive HPAI; secondary respiratory co-infection)
lamp-forge poultry-risk --aiv --ibv

# NDV + IBDV co-detected (velogenic NDV in immunosuppressed flock)
lamp-forge poultry-risk --ndv --ibdv

# Read flags from a JSON file produced by the BioID device integration
lamp-forge poultry-risk \
  --input-json device_output/run_20260601_flock42.json \
  --out-json results/poultry/flock42_assessment.json \
  --out-csv  results/poultry/flock42_assessment.csv
```

The command prints a colour-coded alert table:

```
Poultry biosecurity LAMP panel results:
  AIV   (M gene)  [-]
  NDV   (M gene)  [-]
  IBDV  (VP2  )   [+]
  IBV   (N gene)  [-]

Alert level : MODERATE  (score 20/100)

Interpretation:
  Infectious bursal disease virus (IBDV / Gumboro) detected via the VP2 assay...
Recommended action:
  Notify the flock health veterinarian. Assess birds aged 2-8 weeks for...
```

When AIV is detected, two additional warnings appear:

```
  REPORT REQUIRED: notify the national veterinary authority before any birds are moved or culled.
  ZOONOTIC RISK: all personnel in the affected house must wear N95 respirators, eye protection, and gloves.
  DEPOPULATION RISK: consult the national authority before any birds are moved; emergency depopulation may be required.
```

**Alert level reference.**

| Score | Level    | Trigger                         | Immediate action |
|-------|----------|---------------------------------|-----------------|
| 60+   | CRITICAL | AIV detected (pan-IAV positive) | PPE, report to national authority, no bird movement |
| 35-59 | HIGH     | NDV, or NDV + IBDV / IBV        | Quarantine, confirmatory PCR, vet call |
| 20-34 | MODERATE | IBDV alone or IBDV + IBV        | Vet assessment, vaccination review |
| 10-19 | LOW      | IBV alone                       | Respiratory management, serotype S1 sequencing |
| 0     | NEGATIVE | All negative                    | Confirm extraction control positive |

**LOD estimate for a cloacal-swab sample** (100 uL swab eluate, 60% RNA
extraction efficiency, 50 uL eluate, 5 uL to RT-LAMP reaction):

```bash
lamp-forge lod \
  --sample-volume 100 \
  --efficiency 0.60 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample volume = 100 x 0.6 x (5/50) = 6 uL.
LOD_95 ~ 3.0 / 6 ~ 500 copies/mL of swab eluate.
AIV viral loads in cloacal swabs from HPAI-infected birds can reach
10^6-10^9 TCID50/mL, providing ample headroom for detection within the
BioID 30-60 min run window.

**Field notes.**

- **AIV positives must be treated as presumptive HPAI** until the HA/NA
  subtype is confirmed by a WOAH Reference Laboratory.  Under national HPAI
  response plans, premises are quarantined and birds must not be moved until
  the authority gives clearance.  Do not depopulate based solely on an on-farm
  LAMP result; premature depopulation can void indemnity claims.
- **NDV positives require F-gene pathotyping** at a reference laboratory to
  distinguish velogenic (notifiable) from lentogenic (vaccine-origin) strains.
  A field-positive NDV LAMP result is a presumptive positive pending this
  confirmation.
- **IBDV immunosuppression amplifies all other infections.** A flock with a
  positive IBDV result is at elevated risk from secondary respiratory viruses
  (including IBV), E. coli, and vaccine failures.  Review the entire health
  programme, not just the IBDV result.
- **IBV serotype diversity.** The N-gene assay detects all IBV serotypes but
  cannot distinguish them.  If the assay is positive, send tracheal swab
  material for S1 gene sequencing to match the circulating serotype against
  current vaccine strains; serotype mismatch is the most common cause of
  vaccine failure.

**References.**

- Swayne DE (2012) Impact of vaccines and vaccination on global control of
  avian influenza. *Avian Dis* 56(4 suppl):818-828.
  doi:10.1637/10138-091511-Review.1
- Alexander DJ (2000) Newcastle disease and other avian paramyxoviruses.
  *Rev Sci Tech OIE* 19(2):443-462. doi:10.20506/rst.19.2.1231
- van den Berg TP (2000) Acute infectious bursal disease in poultry:
  a review. *Avian Pathol* 29(3):175-194. doi:10.1080/030794500412385
- Cavanagh D (2007) Coronavirus avian infectious bronchitis virus.
  *Vet Res* 38(2):281-297. doi:10.1051/vetres:2006055
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

## Recipe 31 -- Poultry biosecurity outbreak trend analysis (`lamp-forge poultry-trend`)

**Goal.** Convert a chronological series of poultry biosecurity LAMP panel
results (Recipe 30) into a single trajectory assessment -- EMERGING,
STABLE_CLEAR, STABLE_ENDEMIC, RESOLVING, or INSUFFICIENT_DATA -- that tells
farm managers and veterinarians whether the flock situation is getting better,
holding steady, or deteriorating.

**Why it matters for BioVind.** The BioID cartridge does not replace an
outbreak-response protocol; it generates a data stream.  A single CRITICAL
result triggers an immediate response; a trend across monthly or weekly
monitoring intervals enables smarter decisions about when to escalate, when
to relax quarantine, and whether a consistently IBV-positive flock is in a
managed endemic state (normal) or genuinely deteriorating.

**Workflow.**

Step 1 -- generate a per-interval JSON file for each monitoring round:

```bash
# Month 1 -- IBV only
lamp-forge poultry-risk --ibv \
  --out-json results/flock-A/2026-01.json

# Month 2 -- IBDV appears alongside IBV
lamp-forge poultry-risk --ibdv --ibv \
  --out-json results/flock-A/2026-02.json

# Month 3 -- AIV detected (presumptive HPAI event)
lamp-forge poultry-risk --aiv \
  --out-json results/flock-A/2026-03.json
```

Step 2 -- compute the trend across all three intervals:

```bash
lamp-forge poultry-trend \
  --poultry-result results/flock-A/2026-01.json \
  --poultry-result results/flock-A/2026-02.json \
  --poultry-result results/flock-A/2026-03.json \
  --out-json results/flock-A/trend_Q1.json \
  --out-csv results/flock-A/trend_Q1.csv
```

Expected output:

```
Poultry biosecurity trajectory: 3 sample(s)

Sample ID    Date          AIV   NDV  IBDV   IBV  Alert
--------------------------------------------------------------
2026-01      2026-01-15    [-]   [-]   [-]   [+]  LOW
2026-02      2026-02-15    [-]   [-]   [+]   [+]  MODERATE
2026-03      2026-03-15    [+]   [-]   [-]   [-]  CRITICAL

Trend direction: EMERGING
  AIV NEWLY DETECTED: notify the national veterinary authority immediately;
  quarantine the affected house before any birds are moved.
  Worst alert level: CRITICAL

Interpretation: WOAH-notifiable pathogen(s) newly detected in the most
  recent sample: AIV. Biosecurity situation is escalating -- immediate
  response required.

Recommended action: Notify the national veterinary authority immediately
  regarding AIV detection. Enforce flock movement restrictions...
```

**Using a monitoring-spreadsheet CSV instead of per-run JSON files.**

If the farm already keeps a monitoring spreadsheet, export it to CSV and
supply it directly:

```
sample_id,date,aiv,ndv,ibdv,ibv
FlockA-2026-01,2026-01-15,0,0,0,1
FlockA-2026-02,2026-02-15,0,0,1,1
FlockA-2026-03,2026-03-15,1,0,0,0
```

```bash
lamp-forge poultry-trend \
  --csv monitoring/flock_A_2026.csv \
  --out-json results/flock-A/trend_2026.json
```

**IBV endemic flock scenario.**

A flock in which IBV has been consistently detected across all monitoring
intervals is flagged as `ibv_endemic: true` in the JSON output.  The
`poultry-trend` command also prints:

```
  IBV endemic flock: consistently positive -- review IBV vaccination
  programme and serotype match with the flock health veterinarian.
```

This is the expected background state for many commercial flocks; the
important signal is whether AIV or NDV breaks through on top of the
endemic IBV baseline.

**Trend direction reference.**

| Direction          | Meaning                                                | Priority action |
|--------------------|--------------------------------------------------------|----------------|
| EMERGING           | New WOAH-notifiable pathogen detected, or burden rising | Immediate notification; escalate biosecurity |
| STABLE_CLEAR       | All negative across all intervals                       | Routine surveillance |
| STABLE_ENDEMIC     | Pathogens consistently present, no trend                | Review vaccination; consult vet |
| RESOLVING          | Burden decreasing, or notifiable pathogen cleared       | Confirm with reference lab before relaxing controls |
| INSUFFICIENT_DATA  | Only one sample available                               | Collect next sample and re-run |

**References.**

- WOAH (2024) List of OIE-notifiable diseases, infections and infestations.
- Swayne DE (2012) Impact of vaccines and vaccination on global control of
  avian influenza. *Avian Dis* 56(4 suppl):818-828.
  doi:10.1637/10138-091511-Review.1
- Alexander DJ (2000) Newcastle disease and other avian paramyxoviruses.
  *Rev Sci Tech OIE* 19(2):443-462. doi:10.20506/rst.19.2.1231
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 32 -- Brucella IS711 LAMP: farm biosecurity and human point-of-care (`lamp-forge brucella-risk` / `lamp-forge brucella-trend`)

**Goal.** Design a pan-*Brucella* LAMP assay targeting the IS711 insertion
sequence, then use `lamp-forge brucella-risk` to convert a portable IS711
LAMP result into a structured alert — species-level interpretation, WOAH
notification flag, and recommended clinical or veterinary action — and
`lamp-forge brucella-trend` to track a herd or patient surveillance series
over time.

**Why it is important.** Brucellosis is one of the world's most prevalent
bacterial zoonoses, with ~500,000 new human infections estimated annually
(WHO).  In livestock it is a WOAH-listed Tier 1 notifiable disease with
major trade and welfare consequences:

- **Bovine brucellosis** (*B. abortus*, cattle): causes late-term abortion
  storms, reduced milk yield, and long-lasting infertility.  Mandatory
  national eradication programmes are in force in the US, EU, and many
  other jurisdictions; a positive pen-side result triggers immediate herd
  quarantine and reporting to the national veterinary authority.
- **Small-ruminant brucellosis** (*B. melitensis*, goats and sheep): the
  dominant cause of human brucellosis globally.  *B. melitensis* has been
  eliminated from the US but remains endemic in North Africa, the Middle
  East, Central Asia, and Southern Europe; imported small ruminants are
  the main re-introduction route.
- **Porcine brucellosis** (*B. suis* biovars 1 and 3): significant zoonotic
  risk in pig workers; biovars 2 and 4 are less zoonotic.
- **Human brucellosis** (undulant fever): febrile illness caused by all
  three species above; requires doxycycline + rifampin combination therapy
  for 6 weeks; missed or late diagnosis leads to relapse and chronic
  illness.

A portable IS711 LAMP assay running on a BioVind BioID device is the
fastest route to a quarantine or public-health notification decision in
all three scenarios — farm biosecurity (cattle, small ruminant, swine)
and human point-of-care (rural clinic, urgent care, farm-worker
occupational health).

**Why IS711.**
IS711 (also annotated IS6501 / ISBm2) is a *Brucella*-specific insertion
sequence present in 5–40 copies per genome across all nine *Brucella*
species.  Its multi-copy nature gives 10–100x lower analytical LOD than
single-copy chromosomal targets (Abd El Wahed et al. 2015 *PLoS NTD*),
making it the gold-standard molecular target for brucellosis diagnostics
in field-deployable settings.  IS711 is not present in *Ochrobactrum* spp.
or other alpha-proteobacterial near-neighbours, so a correctly designed
assay is genus-specific.

The internal conserved region used in published LAMP assays shows <2%
nucleotide divergence across *B. abortus*, *B. melitensis*, *B. suis*,
*B. ovis*, *B. canis*, and newer marine *Brucella* species — exceptional
conservation for an isothermal assay target.

**Note on species discrimination.** IS711 LAMP is a **genus-level** assay.
It cannot distinguish *B. abortus* from *B. melitensis* from *B. suis*.
`lamp-forge brucella-risk` infers the most probable species from the
*sample context* (cattle, small ruminant, swine, human, or environmental)
and adjusts the alert level and recommended action accordingly.  For
definitive species identification, follow up with culture or a species-
discriminating PCR at a reference laboratory.

**Target sequences.** NCBI taxon 234 (*Brucella* genus), gene `IS711`,
`max_sequences: 30` retrieves sequences spanning all major species and
biovars: *B. abortus* biovars 1–7/9, *B. melitensis* biovars 1–3,
*B. suis* biovars 1–5, *B. ovis*, and *B. canis*.

**Off-target panel.** The most important off-targets for IS711 are other
alpha-proteobacteria that share genomic synteny with *Brucella*:

| File | Source | Why |
|---|---|---|
| `ochrobactrum_anthropi.fasta` | ATCC 49188 RefSeq | Closest IS711-free relative; reclassified *Brucella anthropi* in 2020; farm and clinical environment |
| `rhizobium_leguminosarum.fasta` | RL3841 RefSeq (NC_008380.1) | Soil alpha-proteobacterium present in farm environments |
| `agrobacterium_tumefaciens.fasta` | C58 RefSeq (NC_003062.2) | Plant pathogen; found in soil and irrigation water on farms |
| `bartonella_henselae.fasta` | Houston-1 RefSeq (NC_005956.1) | Human-infecting alpha-proteobacterium; clinical differential for human POC samples |
| `bovine_genome_fragment.fasta` | Subset of ARS-UCD1.3 chr1 | Host DNA from bovine blood or tissue samples; dominant nucleic acid in animal samples |
| `human_chr_fragment.fasta` | Subset of GRCh38 chr1 | Host DNA from human blood or serum samples |

**Config.** Use the ready-made config at `config/brucella_is711.yaml`, or
paste the block below:

```yaml
target:
  name: brucella_is711
  taxon_id: 234              # Brucella genus -- pan-species (abortus, melitensis, suis, ovis, canis)
  gene: IS711                # insertion sequence IS711 (also annotated IS6501 / ISBm2)
  max_sequences: 30          # span major Brucella spp. and biovars
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85
  min_coverage_threshold: 0.85

conservation:
  window_size: 30
  entropy_threshold: 0.20
  min_region_length: 220

primer:
  tm_min: 60.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 45.0               # Brucella genome ~57% GC; IS711 is ~55-60% GC
  gc_max: 68.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/brucella_is711
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `gc_min: 45.0` -- *Brucella* has one of the highest GC contents among
  animal pathogens (~57% genomic GC).  IS711 itself is ~55–60% GC.  The
  standard LAMP-Forge default of `gc_min: 40` admits primers that are
  too AT-rich for this target; raising to 45 keeps candidates realistic.
  Unlike the AT-rich targets in Recipe 22 (*C. difficile*, `gc_min: 25`)
  and Recipe 5 (KPC, `gc_min: 30`), *Brucella* is at the high-GC end of
  the primer design space.

- `gc_max: 68.0` -- Caps primers below the GC ceiling of the IS711
  transposase ORF (~68% GC in the most GC-dense sub-regions) while
  leaving room for all biologically reasonable candidates.

- `tm_min: 60.0` -- IS711 is a chromosomal **DNA** target; standard LAMP
  (not RT-LAMP) is appropriate.  Do **not** raise to the 63 degC one-step
  RT-LAMP co-activity floor used for RNA-target assays (PRRSV, FMDV,
  AIV, CSFV, NDV, PEDV, PDCoV in Recipes 11-14, 18, 27-29).

- `entropy_threshold: 0.20` -- IS711 is exceptionally conserved within
  the genus; the internal region used in published LAMP assays shows <2%
  inter-species divergence, so the standard 0.20-bit entropy threshold
  (same as well-conserved housekeeping genes rpoB and speB) is appropriate.

- `min_identity_threshold: 0.85` / `min_coverage_threshold: 0.85` -- The
  most important off-target is *Ochrobactrum anthropi* (reclassified as
  *Brucella anthropi* in 2020 by some authorities), which shares ~70-80%
  average nucleotide identity with *Brucella* but lacks IS711.  The
  0.85/0.85 threshold is the minimum needed to surface any residual
  cross-reactive candidate before wet-lab validation.

**Run the primer design pipeline.**

```bash
docker compose run --rm lamp-forge run --config /work/config/brucella_is711.yaml
```

**Estimate LOD before ordering** (whole blood 200 uL, 60% DNA extraction
efficiency, 50 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 200 \
  --efficiency 0.60 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 200 x 0.60 x (5/50) = 12 uL per reaction -> LOD_95
approx 200 copies/mL of whole blood.  *Brucella* bacteraemia during the
acute febrile phase can reach 10^3-10^5 CFU/mL; the multi-copy IS711
(5-40 copies per cell) further multiplies effective target molecules ~20-fold,
putting the functional LOD closer to 10-40 *Brucella* CFU/mL -- well within
the clinically detectable range.

**Check TTP for the BioVind device window.**

```bash
lamp-forge ttp \
  --copies 200 \
  --preset lamp
```

At 200 copies/reaction (with IS711 multi-copy amplification equivalent to
~10-40 cells) the predicted TTP is ~43 min -- inside the 60-min BioVind
BioID window with headroom for sample-matrix inhibition.

**Verify config before running.**

```bash
lamp-forge validate \
  --config config/brucella_is711.yaml \
  --no-check-dirs
```

**Export primers for ordering.**

```bash
lamp-forge export \
  --input results/brucella_is711/primer_sets.json \
  --format idt \
  --target-label Brucella_IS711 \
  --out orders/brucella_is711_idt_order.csv
```

**Interpreting a result with `lamp-forge brucella-risk`.**

IS711 detects *Brucella* genus but the clinical or veterinary action depends
on sample context.  `lamp-forge brucella-risk` accepts the IS711 LAMP result
(`--is711` / `--no-is711`) and the sample context (`--context`) and returns a
structured alert:

```bash
# Positive result from a cattle pen-side blood sample
lamp-forge brucella-risk --is711 --context cattle

# Positive result from a small-ruminant (goat/sheep) serum sample
lamp-forge brucella-risk --is711 --context small_ruminant

# Positive result from a swine herd blood sample
lamp-forge brucella-risk --is711 --context swine

# Positive result from a human patient (rural clinic / farm-worker POC)
lamp-forge brucella-risk --is711 --context human

# Negative result -- all contexts give NEGATIVE alert
lamp-forge brucella-risk --no-is711 --context cattle

# Environmental sample (water trough, fomite wipe, soil)
lamp-forge brucella-risk --is711 --context environmental

# Save structured output for downstream workflows or trend analysis
lamp-forge brucella-risk \
  --is711 \
  --context cattle \
  --out-json results/herd-A/2026-01-is711.json \
  --out-csv  results/herd-A/2026-01-is711.csv
```

**Alert level reference.**

| Context | IS711 | Alert | Score | Notifiable | Zoonotic | Immediate report |
|---|---|---|---|---|---|---|
| Human | Positive | CRITICAL | 90 | Yes | Yes | Yes |
| Cattle | Positive | HIGH | 60 | Yes | Yes | Yes |
| Small ruminant | Positive | HIGH | 60 | Yes | Yes | Yes |
| Swine | Positive | HIGH | 60 | Yes | Yes | Yes |
| Environmental | Positive | LOW | 25 | No | No | No |
| Any | Negative | NEGATIVE | 0 | No | No | No |

A CRITICAL result (human context) triggers the public-health notification
pathway: notify the local health department immediately and start empiric
doxycycline + rifampin while confirmatory serology is sent to the reference
laboratory.  HIGH results (livestock) trigger the national veterinary
authority notification pathway and immediate herd quarantine.  LOW results
(environmental) require follow-up animal testing but not immediate reporting.

**Tracking a herd surveillance series with `lamp-forge brucella-trend`.**

Use `lamp-forge brucella-trend` to convert a chronological sequence of IS711
results into a trend assessment: NEWLY_DETECTED, PERSISTENT, RESOLVING,
STABLE_CLEAR, or INSUFFICIENT_DATA.  This drives the key operational decisions:
when to notify the veterinary authority, when to conduct whole-herd serology,
and when it is safe to apply for quarantine release.

*From individual JSON files (one file per monitoring interval):*

```bash
# Monthly cattle surveillance -- four consecutive months
lamp-forge brucella-risk --is711   --context cattle \
  --out-json results/herd-A/2026-01.json
lamp-forge brucella-risk --is711   --context cattle \
  --out-json results/herd-A/2026-02.json
lamp-forge brucella-risk --no-is711 --context cattle \
  --out-json results/herd-A/2026-03.json
lamp-forge brucella-risk --no-is711 --context cattle \
  --out-json results/herd-A/2026-04.json

lamp-forge brucella-trend \
  --brucella-result results/herd-A/2026-01.json \
  --brucella-result results/herd-A/2026-02.json \
  --brucella-result results/herd-A/2026-03.json \
  --brucella-result results/herd-A/2026-04.json \
  --out-json results/herd-A/trend_Q1.json
```

Expected trend output for the series above (positive, positive, negative,
negative -- RESOLVING):

```
Brucella IS711 surveillance trajectory: 4 sample(s)

Sample ID      Date          IS711  Context         Alert
-----------------------------------------------------------
                              [+]    cattle          HIGH
                              [+]    cattle          HIGH
                              [-]    cattle          NEGATIVE
                              [-]    cattle          NEGATIVE

Trend direction: RESOLVING
  IS711 cleared: confirm with serology before lifting quarantine.
  Worst alert level: HIGH
```

*From a monitoring-spreadsheet CSV (header required):*

```
sample_id,date,is711_positive,context
HerdA-2026-01,2026-01-15,1,cattle
HerdA-2026-02,2026-02-15,1,cattle
HerdA-2026-03,2026-03-15,0,cattle
HerdA-2026-04,2026-04-15,0,cattle
```

```bash
lamp-forge brucella-trend \
  --csv monitoring/herd_A_2026.csv \
  --out-json results/herd-A/trend_2026.json
```

**Trend direction reference.**

| Direction | Meaning | Priority action |
|---|---|---|
| NEWLY_DETECTED | IS711 positive for the first time in this series | Notify national veterinary authority; initiate herd quarantine |
| PERSISTENT | IS711 positive in >= 75% of >= 3 intervals | Whole-herd serology to identify and remove source animal |
| RESOLVING | IS711 was positive, now negative | Confirm clearance with serology before applying for quarantine release |
| STABLE_CLEAR | All intervals negative | Routine quarterly surveillance; no action required |
| INSUFFICIENT_DATA | Only one sample available | Collect next sample and re-run |

**Multiplex IS711 with farm biosecurity targets** for a single-run differential
panel on the BioVind device:

```bash
lamp-forge panel \
  --set Brucella_IS711=results/brucella_is711/primer_sets.json \
  --set ASFV=results/asfv_p72/primer_sets.json \
  --set FMDV=results/fmdv_3dpol/primer_sets.json \
  --top-per-target 5 \
  --out results/farm_4plex_brucella
```

Brucella IS711 is a DNA LAMP assay running at 60-65 degC, compatible with
co-incubation alongside one-step RT-LAMP assays for RNA viruses (ASFV is a
large DNA virus, also compatible).  Use `lamp-forge cartridge` to assign
channel slots:

```bash
lamp-forge cartridge \
  --panel results/farm_4plex_brucella/panel.json \
  --out results/farm_4plex_brucella/cartridge_build.json
```

**Field and safety notes.**

- *Brucella* spp. are BSL-3 pathogens (CDC) and USDA Select Agents
  (*B. abortus*, *B. melitensis*, *B. suis*).  All live-culture work must
  be conducted in a BSL-3 facility.  LAMP on extracted DNA or
  heat-inactivated samples is appropriate for a BSL-2 or field setting.
- Sample types for livestock: EDTA whole blood, serum, milk (ring test),
  vaginal swab post-abortion, fetal stomach contents.  For human POC:
  EDTA whole blood or serum from febrile patients with livestock exposure.
- A positive IS711 LAMP result is a presumptive positive for *Brucella*
  genus.  **Do not use this result alone for regulatory reporting** --
  confirmatory testing at an approved reference laboratory (culture,
  species-specific PCR, complement fixation test or ELISA) is required
  before quarantine orders or trade restrictions.
- Run a no-template control (NTC) and a positive extraction control
  (purified *Brucella* reference-strain DNA or heat-killed *B. abortus*
  standard, e.g. ATCC 23448) in every batch.
- Include the universal 16S rRNA internal extraction control (Recipe 19)
  to confirm that inhibition is not causing false-negative results in
  whole-blood or milk samples.
- Environmental samples (water troughs, fomite wipes, soil) may contain
  PCR-detectable IS711 DNA from persistence of dead *Brucella* cells or
  free DNA.  A positive environmental result does **not** confirm active
  infection in animals; follow up with animal pen-side testing before
  regulatory reporting.

**References.**

- Abd El Wahed A, El-Deeb A, El-Tholoth M et al. (2015) Loop-mediated
  isothermal amplification on a chip for detection of *Brucella* spp. in
  milk samples. *PLoS Negl Trop Dis* 9(10):e0004124.
  doi:10.1371/journal.pntd.0004124
- Al-Dahouk S, Nockler K, Scholz HC et al. (2010) Evaluation of
  *Brucella* IS711-based loop-mediated isothermal amplification for
  human brucellosis. *Ann Clin Microbiol Antimicrob* 9:15.
  doi:10.1186/1476-0711-9-15
- Godfroid J, Cloeckaert A, Liautard JP et al. (2005) From the discovery
  of the Malta fever's agent to the discovery of a marine mammal reservoir,
  brucellosis has continuously been a re-emerging zoonosis.
  *Vet Res* 36(3):313-326. doi:10.1051/vetres:2005003
- Godfroid J, Nielsen K, Saegerman C (2010) Diagnosis of brucellosis in
  livestock and wildlife. *Croat Med J* 51(4):296-305.
  doi:10.3325/cmj.2010.51.296
- WOAH (2021). Bovine brucellosis. *Terrestrial Animal Health Code*,
  Chapter 8.4; *Terrestrial Manual*, Chapter 3.1.2.
- WHO (2006). Brucellosis in humans and animals. WHO/CDS/EPR/2006.7.
  World Health Organization, Geneva.
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 33 -- Neonatal calf diarrhea outbreak trajectory (`lamp-forge bov-enteric-trend`)

**Goal.** Convert a weekly LAMP monitoring series for a calf group into a
trend assessment -- EMERGING, RESOLVING, STABLE_CLEAR, or STABLE_ENDEMIC --
with automatic detection of ETEC K99 new-onset events (IV-antibiotic emergency)
and Salmonella persistence (zoonotic herd-investigation trigger).

**Why it matters for BioVind.**
Neonatal calf diarrhea (scours) is the leading cause of calf mortality
globally.  A single LAMP result answers "what pathogen is present now?"; the
trend answers "is this getting better or worse?" and surfaces two
time-critical events that change the veterinary response:

1. **ETEC K99 newly detected** -- near-certain mortality in calves <4 days
   old without IV antibiotics within 2 h.  No other pathogen in the panel
   has this time-pressure, and a first detection after clean intervals is
   the highest-urgency escalation in neonatal medicine.
2. **Salmonella persistent** -- positive in every monitored interval suggests
   a herd-level contamination source or a chronic shedder; culture and herd
   investigation are required.

The trend command also reports:

- **Zoonotic risk count** -- intervals where Salmonella or Cryptosporidium
  is positive (both are notifiable food-safety pathogens).
- **Antibiotic required count** -- intervals where ETEC K99 or Salmonella
  is detected (the only two targets in the panel that require antimicrobial
  therapy; prevents overuse in viral/parasitic scours).
- **Worst alert level** across the monitoring series.

**Panel targets (same as `lamp-forge bov-enteric-risk`).**

| Target     | Gene      | Pathogen                        | Alert if sole positive |
|---|---|---|---|
| ETEC_K99   | faeG      | ETEC K99/F5 enterotoxigenic E. coli | CRITICAL |
| Salmonella | invA      | Salmonella sp. (any serovar)    | HIGH |
| Crypto     | 18S/hsp70 | Cryptosporidium parvum          | MODERATE |
| BCoV       | N gene    | Bovine coronavirus (enteric)    | MODERATE |
| RotaA      | VP6       | Bovine rotavirus A              | LOW |

**Typical workflow (JSON inputs).**

```bash
# Week 1 -- rotavirus only
lamp-forge bov-enteric-risk --rota-a \
  --out-json results/calf-grp-1/week-01.json

# Week 2 -- rotavirus + coronavirus
lamp-forge bov-enteric-risk --rota-a --bcov \
  --out-json results/calf-grp-1/week-02.json

# Week 3 -- ETEC K99 onset (outbreak event)
lamp-forge bov-enteric-risk --etec-k99 --rota-a \
  --out-json results/calf-grp-1/week-03.json

# Trend analysis across the three weeks
lamp-forge bov-enteric-trend \
  --enteric-result results/calf-grp-1/week-01.json \
  --enteric-result results/calf-grp-1/week-02.json \
  --enteric-result results/calf-grp-1/week-03.json \
  --out-json results/calf-grp-1/trend_wk1-3.json
```

Expected output:

```
Bovine neonatal enteric trajectory: 3 sample(s)

Sample ID    Date          ETEC  SALM  CRPT  BCOV  ROTA  Alert
-----------------------------------------------------------------
week-01                     [-]   [-]   [-]   [-]   [+]  LOW
week-02                     [-]   [-]   [-]   [+]   [+]  MODERATE
week-03                     [+]   [-]   [-]   [+]   [+]  CRITICAL

Trend direction: EMERGING
  ETEC K99 NEWLY DETECTED: IV antibiotic therapy required within 2 h; ...
  Worst alert level: CRITICAL

Interpretation: ETEC K99 (faeG) newly detected in the most recent sample
  after 2 consecutive negative interval(s). ETEC K99 causes near-certain
  mortality (>90%) in calves <4 days old without IV antibiotic therapy ...

Recommended action: CRITICAL: Start IV or oral electrolyte rehydration
  immediately. Administer IV antibiotics (ampicillin or amoxicillin-
  clavulanate) within 2 h of diagnosis -- do NOT wait for sensitivity
  results. Isolate affected calves. Sanitise pens ...
```

**Typical workflow (CSV input -- monitoring spreadsheet).**

The CSV route ingests a monitoring spreadsheet directly.
Column order is flexible; only the header names matter.

```
sample_id,date,etec_k99,salmonella,crypto,bcov,rota_a
CG1-week01,2026-01-07,0,0,0,0,1
CG1-week02,2026-01-14,0,0,0,1,1
CG1-week03,2026-01-21,1,0,0,1,1
CG1-week04,2026-01-28,0,0,0,0,1
CG1-week05,2026-02-04,0,0,0,0,0
```

```bash
lamp-forge bov-enteric-trend \
  --csv monitoring/calf_grp1_2026.csv \
  --out-json results/calf-grp-1/trend_2026.json \
  --out-csv  results/calf-grp-1/timeline_2026.csv
```

**Scenario: Salmonella persistent (zoonotic herd investigation).**

If Salmonella is positive at every monitoring interval, the trend command
raises a SALMONELLA PERSISTENT flag alongside the STABLE_ENDEMIC direction.
This pattern suggests either a herd-level environmental contamination source
(contaminated water trough, feed supply) or a persistently shedding animal
in the group.

```bash
# Build three Salmonella-positive JSON files programmatically
for w in 01 02 03; do
  lamp-forge bov-enteric-risk --salmonella \
    --out-json results/calf-grp-2/week-${w}.json
done

lamp-forge bov-enteric-trend \
  --enteric-result results/calf-grp-2/week-01.json \
  --enteric-result results/calf-grp-2/week-02.json \
  --enteric-result results/calf-grp-2/week-03.json
```

Expected key lines in the output:

```
Trend direction: STABLE_ENDEMIC
  SALMONELLA PERSISTENT: detected in every interval -- herd-level
    contamination or chronic shedder suspected; submit for culture and
    herd investigation.
  Zoonotic risk in 3 interval(s) (Salmonella and/or Cryptosporidium):
    enforce gloves and handwashing.
  Antibiotics indicated in 3 interval(s) (ETEC K99 and/or Salmonella).
```

**Output files.**

| File | Content |
|---|---|
| `--out-json` | Full trend assessment JSON (direction, event flags, interpretation, timeline) |
| `--out-csv` | Chronological timeline with per-pathogen 0/1 columns and per-sample alert level |

**Integration with the BioID multiplex panel.**

The BioID device's five-channel enteric panel maps directly to the five
targets above.  After each cartridge result, export the positivity flags to
a JSON file with `lamp-forge bov-enteric-risk --out-json` (flags can also be
set via `--input-json` from an upstream data integration), then feed the
JSON files into `lamp-forge bov-enteric-trend` to build the surveillance
timeline.

For whole-herd integration, see Recipe 28 (`lamp-forge cartridge`) which
packages all design artefacts into a single cartridge manifest ready for
the BioID instrument build.

**References.**

- Cho Y & Yoon K (2014) An overview of calf diarrhea -- infectious
  etiology, diagnosis, and intervention.  *J Vet Sci* 15(1):1-17.
  doi:10.4142/jvs.2014.15.1.1
- Gulliksen SM et al. (2009) Enteropathogens and risk factors for
  diarrhea in Norwegian dairy calves.
  *J Dairy Sci* 92(10):5057-5066. doi:10.3168/jds.2009-2080
- Smith DR (2012) Field diseases of beef calves in the first 45 days of
  life: enteric and respiratory diseases.
  *Vet Clin Food Anim* 28(1):57-82. doi:10.1016/j.cvfa.2012.01.001
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 34 -- Human Respiratory Syncytial Virus (hRSV) via RT-LAMP (human point-of-care)

**Goal.** Detect human RSV broadly across both RSV-A and RSV-B subtypes from
a nasal swab at a rural urgent care clinic, pharmacy, or telemedicine-adjacent
facility -- fast enough to guide antiviral therapy decisions before the patient
leaves the building.

**Why it's interesting.** hRSV is the **leading cause of acute lower
respiratory tract illness in infants** globally -- an estimated 33 million
RSV-associated ALRI episodes and 3.6 million hospitalisations per year in
children under 5 (Shi et al. Lancet 2017).  In older adults and
immunocompromised patients, RSV causes disease comparable to influenza A.
From 2023 onward, two RSV vaccines (Abrysvo, Arexvy) have been approved for
older adults, and a maternal vaccine (Abrysvo) for infant protection --
creating a new clinical-decision point: does this patient need antiviral
nirsevimab (monoclonal antibody, paediatric) or is a respiratory illness due to
influenza A (oseltamivir)?  A portable LAMP panel that simultaneously reports
hRSV + influenza A + SARS-CoV-2 resolves the differential in under 60 minutes
at a BioVind-style bedside device, directing the right antiviral to the right
patient within the treatment window.

**Why this is a BioVind fit.**  hRSV has historically been under-tested at
point-of-care because rapid antigen tests have only moderate sensitivity
(~80%) and PCR is not available outside of hospital labs.  A portable
RT-LAMP panel changes this equation for rural clinics, urgent cares, and
even home-based telehealth monitoring -- exactly BioVind's POC vertical.

**Two subtypes, one assay.**  RSV-A and RSV-B circulate simultaneously each
winter season, with inter-subtype nucleotide identity of ~77% in the N gene.
A pan-RSV design requires both subtypes in the input sequence set.  The N gene
(nucleocapsid, ~1203 nt CDS) is the most conserved coding region across both
subtypes -- more so than the fusion (F) or attachment (G) genes -- and is the
basis of published one-step RT-LAMP assays (Nie et al. 2017; Zhang et al. 2021).

**hRSV is a negative-sense ssRNA pneumovirus** (Pneumoviridae,
Orthopneumovirus) -- so detection requires **RT-LAMP**.  One-step RT-LAMP
(Bst 2.0 WarmStart + NEB RTx at 63-65 degC) delivers a yes/no result in under
30 minutes, fully compatible with the BioVind BioID portable platform.

**Target.** N gene (nucleocapsid), anchoring on RSV-A (NCBI taxon 11250) and
adding RSV-B accessions explicitly for cross-subtype conservation.  Use the
ready-made config at `config/hrsv_N_gene.yaml`, or paste the block below:

```yaml
target:
  name: hrsv_N_gene
  taxon_id: 11250          # Human orthopneumovirus (RSV-A anchor species)
  accessions:
    # Add RSV-B accessions for pan-RSV cross-subtype conservation.
    # Representative starting points:
    #   AY353550.1   RSV-B prototype strain B1 (N-gene region)
    #   KJ627327.1   RSV-B contemporary BA9 genotype
    # Pull additional accessions spanning RSV-A clades A1, A2, A2b and
    # RSV-B lineages BA7-BA10 for full pan-RSV coverage.
    []
  gene: N                  # nucleocapsid protein; most conserved across RSV-A and RSV-B
  max_sequences: 35        # span RSV-A clades A1/A2/A2b + RSV-B BA lineages
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets
  min_identity_threshold: 0.85   # tight -- hMPV shares Pneumoviridae family
  min_coverage_threshold: 0.80

conservation:
  window_size: 25
  entropy_threshold: 0.30        # RSV N gene varies ~23% between subtypes; RNA-virus drift
  min_region_length: 200

primer:
  tm_min: 63.0                   # RT-LAMP one-step: co-activity floor for NEB RTx / Bst 2.0
  tm_max: 65.0
  tm_match_tolerance: 2.0
  gc_min: 35.0                   # RSV N gene ~37-42% GC; AT-rich pneumovirus
  gc_max: 55.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 150

output:
  dir: results/hrsv_N_gene
  top_n: 10
  generate_html: true
  generate_csv: true
```

**Key config notes.**

- `primer.tm_min: 63.0` -- one-step RT-LAMP at 63-65 degC; same rationale as
  PRRSV (Recipe 11), AIV (Recipe 12), FMDV (Recipe 14), NDV (Recipe 18), PEDV
  (Recipe 27), and CSFV (Recipe 29).  hRSV is negative-sense ssRNA; the RT
  enzyme must remain co-active with Bst polymerase at 63-65 degC.
- `conservation.entropy_threshold: 0.30` -- RSV-A and RSV-B share only ~77%
  N-gene nt identity.  A threshold below 0.20 bits would be too strict for
  cross-subtype conservation analysis; 0.30 bits finds windows conserved across
  both subtypes while excluding the highly variable central portion of the N gene.
- `gc_min: 35.0, gc_max: 55.0` -- the hRSV N gene is approximately 37-42% GC
  (RSV is notably AT-rich relative to most human respiratory RNA viruses).  The
  20-unit window accommodates inter-subtype GC variation without admitting primers
  too far outside the primer3 design space.
- `off_targets.min_identity_threshold: 0.85` -- human metapneumovirus (hMPV,
  Pneumoviridae) shares the paramyxovirus-like nucleocapsid fold with RSV.  At
  0.85, the threshold is strict enough to flag hMPV cross-reactivity risk while
  permitting the RSV-A/B diversity needed for pan-RSV coverage.

**Off-target panel.** Co-circulating respiratory pathogens and human host DNA
from nasal swab:

| File | Source | Why |
|---|---|---|
| `hmpv.fasta` | Human metapneumovirus (NC_004148.2) | Closest relative; Pneumoviridae family; shares N-protein structural motifs |
| `hpiv1.fasta` | Human parainfluenza virus 1 (NC_003461.1) | Co-circulating Paramyxoviridae; causes croup in same age group as RSV |
| `hpiv3.fasta` | Human parainfluenza virus 3 (NC_001796.2) | Most common PIV; bronchiolitis in infants overlaps with RSV presentation |
| `influenza_a_h1n1.fasta` | Influenza A H1N1 *(representative segments)* | Primary co-circulating respiratory pathogen in RSV season; separate panel target |
| `sars_cov_2.fasta` | SARS-CoV-2 N gene *(representative)* | Co-circulating respiratory pathogen; separate panel target (Recipe 2) |
| `human_chr_fragment.fasta` | Subset of GRCh38 chr1 | Host DNA from nasal swab; dominant background matrix |

**Expected output.** The N gene (~1203 nt) typically yields 2-4 conserved
windows across RSV-A and RSV-B.  With entropy_threshold 0.30, expect one major
window in the N-terminal RNA-binding domain (~nt 50-400 of the N CDS) and one
in the C-terminal region (~nt 800-1100).  Published RT-LAMP primer sets (Nie
et al. 2017) land in the N-terminal domain -- a top-scoring LAMP-Forge set
should co-localise there.  Zero flagged off-target hits are expected against
hMPV at 0.85 identity.

**Verify RT-LAMP readiness.**

```bash
lamp-forge rt-check \
  --input results/hrsv_N_gene/primer_sets.json \
  --na-type rna \
  --out-csv results/hrsv_N_gene/rt_check.csv
```

Sets marked **NOT OPTIMIZED** have at least one core primer below 63.0 degC.
Because `primer.tm_min: 63.0` is already set, all sets should pass -- verify
before ordering.

**Estimate LOD before ordering** (nasal swab in 400 uL VTM, 50% RNA
extraction, 50 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 400 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 400 x 0.5 x (5/50) = 20 uL per reaction.  LOD_95 approx
150 copies/mL of VTM.  hRSV in symptomatic patients typically ranges from
10^3 to 10^7 copies/mL in nasal washes, providing orders-of-magnitude
headroom above the achievable LOD.

**Multiplex with other POC respiratory targets.** hRSV, influenza A, and
SARS-CoV-2 are the three major viral causes of lower respiratory tract illness
in winter.  Design each target independently, then combine:

```bash
lamp-forge panel \
  --set RSV=results/hrsv_N_gene/primer_sets.json \
  --set IAV=results/influenza_a_M/primer_sets.json \
  --set SARS2=results/sars_cov_2_N/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/poc_respiratory_3plex
```

Then generate the pooling sheet (3 RNA targets = 132 uM minimum; request
200 uM resuspension from IDT or Twist) and export to vendor:

```bash
lamp-forge pool \
  --panel results/poc_respiratory_3plex/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/poc_respiratory_3plex/pool_sheet.csv

lamp-forge export \
  --input results/hrsv_N_gene/primer_sets.json \
  --format idt \
  --target-label RSV_N \
  --out orders/hrsv_N_idt_order.csv
```

**Confirm TTP within the BioVind device window.**

```bash
lamp-forge ttp --preset rt-lamp --device-window 60
```

At LOD_95 (~3 copies/reaction), TTP is expected around 50-60 min for
one-step RT-LAMP; clinical specimens at 10^3+ copies/reaction read out in
30-45 min, well within the 60-min BioVind BioID window.

**Wet-lab note.** The nasal swab matrix (mucus, human genomic DNA) is
moderately inhibitory.  Validate with a spiked VTM extraction protocol
(e.g. MagMAX CORE Nucleic Acid Purification Kit) and confirm LOD using
contrived specimens at 10^2-10^4 copies/mL before clinical deployment.
Note that hRSV samples may require a lysis step at 95 degC for 5 min to
denature the negative-sense genome before the RT step when using column-free
extraction.

**References.**

- Shi T, McAllister DA, O'Brien KL et al. (2017) Global, regional, and
  national disease burden estimates of acute lower respiratory infections due
  to respiratory syncytial virus in young children in 2015: a systematic review
  and modelling study.  *Lancet* 390(10098):946-958.
  doi:10.1016/S0140-6736(17)30938-8
- Nie K, Zhang W, Han J et al. (2017) Rapid and sensitive detection of human
  respiratory syncytial virus using one-step reverse transcription loop-mediated
  isothermal amplification.  *Arch Virol* 162(4):1077-1083.
  doi:10.1007/s00705-016-3182-3
- Zhang Y, Chen Q, Zhang X et al. (2021) Rapid detection of human respiratory
  syncytial virus by reverse transcription loop-mediated isothermal amplification.
  *J Clin Lab Anal* 35(3):e23687. doi:10.1002/jcla.23687
- Zlateva KT, Lemey P, Moes E et al. (2004) Genetic variability and molecular
  evolution of the human respiratory syncytial virus subgroup B attachment G
  protein.  *J Virol* 78(9):4675-4683. doi:10.1128/JVI.78.9.4675-4683.2004
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

## Recipe 35 -- Swine Respiratory Disease Complex (PRDC) three-target panel (`lamp-forge swine-respiratory-risk`)

**Goal.** A portable BioVind-style three-channel LAMP panel that distinguishes
the three core drivers of the Porcine Respiratory Disease Complex (PRDC):
PRRSV (RNA), PCV2 (DNA), and *Mycoplasma hyopneumoniae* (DNA).  Together,
these three pathogens account for the majority of swine respiratory mortality
and production loss in commercial grow-finish pigs worldwide.

**Why it matters.**  PRDC is the most economically significant disease complex
in commercial pig production -- USD 664 M/yr from PRRSV alone in the US
(Holtkamp et al. 2013), EUR 385 M/yr from PCV2 in the EU (Segalés et al.
2005), and over USD 1 B/yr globally from *M. hyopneumoniae* (Maes et al.
2018).  The three pathogens are not independent: PRRSV destroys alveolar
macrophages, PCV2 further dysregulates innate immunity, and Mhp eliminates
the mucociliary escalator -- together producing mortality and production
losses far exceeding the sum of their individual effects.  A portable pen-side
result distinguishes "PRRSV only (restrict movement)" from "full PRDC
(restrict + treat + vaccinate + call vet)" at the point of care, before the
situation escalates.

### Targets

| Config | Pathogen | NA type | Gene | Assay |
|---|---|---|---|---|
| `prrsv_orf7.yaml` | PRRSV (pan-genotype) | ssRNA | ORF7 | RT-LAMP |
| `pcv2_cap.yaml` | Porcine circovirus 2 (pan-genotype) | ssDNA | cap (ORF2) | LAMP |
| `mycoplasma_hyop_p97.yaml` | *Mycoplasma hyopneumoniae* | DNA | p97 adhesin | LAMP |

### Designing primers for each target

Design PRRSV, PCV2, and Mhp primers independently (one config per run):

```bash
# PRRSV (RT-LAMP -- RNA virus)
lamp-forge run --config config/prrsv_orf7.yaml
lamp-forge rt-check --input results/prrsv_ORF7/primer_sets.json --na-type rna

# PCV2 (standard LAMP -- ssDNA virus, dsDNA replication intermediate)
lamp-forge run --config config/pcv2_cap.yaml

# Mycoplasma hyopneumoniae (standard LAMP -- DNA bacterium)
lamp-forge run --config config/mycoplasma_hyop_p97.yaml
```

### Multiplex panel cross-dimer check

Before pooling, confirm the three primer sets are compatible:

```bash
lamp-forge panel \
  --set PRRSV=results/prrsv_ORF7/primer_sets.json \
  --set PCV2=results/pcv2_cap/primer_sets.json \
  --set Mhp=results/mycoplasma_hyop_p97/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/prdc_panel
```

### Pooling the working stock

PRRSV is RNA (RT-LAMP); PCV2 and Mhp are DNA (LAMP).  For a one-step
reaction that runs all three simultaneously, use a WarmStart LAMP Master Mix
that contains RTx (e.g. NEB WarmStart LAMP Kit E1700 at 65 degC).  The pool
uses standard LAMP stoichiometry per target:

```bash
lamp-forge pool \
  --panel results/prdc_panel/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/prdc_panel/prdc_pool_sheet.csv
```

At 3 targets, minimum stock concentration is 3 x 44 = 132 uM; request
200 uM resuspension from IDT or Twist.

### Interpreting a pen-side result

Feed the positivity flags from the device into the risk engine:

```bash
# Example -- PRRSV and PCV2 co-detected (classic PRDC)
lamp-forge swine-respiratory-risk --prrsv --pcv2

# Example -- Mhp alone (endemic enzootic pneumonia)
lamp-forge swine-respiratory-risk --mhp

# Save for trend tracking
lamp-forge swine-respiratory-risk --prrsv --pcv2 \
  --out-json results/barn-A/2026-01.json
```

**Scoring model and alert levels:**

| Pattern | Score | Alert |
|---|---|---|
| PRRSV + PCV2 + Mhp | 90 | CRITICAL |
| PRRSV + PCV2 | 70 | CRITICAL |
| PRRSV + Mhp | 60 | CRITICAL |
| PCV2 + Mhp | 50 | HIGH |
| PRRSV alone | 40 | HIGH |
| PCV2 alone | 30 | MODERATE |
| Mhp alone | 20 | LOW |
| All negative | 0 | NEGATIVE |

PRRSV triggers `movement_restriction_required`.  Two or more positives
triggers `complex_prdc` -- the multi-pathogen synergy warning.

### Tracking trend across monitoring intervals

```bash
lamp-forge swine-respiratory-trend \
  --resp-result results/barn-A/2026-01.json \
  --resp-result results/barn-A/2026-02.json \
  --resp-result results/barn-A/2026-03.json \
  --out-json results/barn-A/prdc_trend_Q1.json

# Or from a monitoring CSV:
lamp-forge swine-respiratory-trend \
  --csv results/barn-A/monitoring.csv \
  --out-json results/barn-A/prdc_trend.json \
  --out-csv results/barn-A/prdc_timeline.csv
```

CSV format:

```
sample_id,date,prrsv,pcv2,mhp
BarnA-2026-01,2026-01-15,1,0,0
BarnA-2026-02,2026-02-15,1,1,0
BarnA-2026-03,2026-03-15,0,1,1
```

Trend directions: EMERGING (PRRSV newly detected or burden rising),
STABLE_CLEAR, STABLE_ENDEMIC (Mhp/PCV2 endemic), RESOLVING.

### Veterinary decision matrix

| Alert | Key flag | Immediate action |
|---|---|---|
| CRITICAL | `complex_prdc=True` | Full PRDC protocol; movement restriction; vet call |
| HIGH | `movement_restriction_required=True` | Restrict movement; assess vaccination; vet call |
| MODERATE | -- | Verify PCV2 vaccination; monitor; vet consultation |
| LOW | -- | Review Mhp metaphylaxis schedule; monitor |
| NEGATIVE | -- | Routine surveillance; confirm extraction control positive |

**References.**

- Holtkamp DJ et al. (2013) Assessment of the economic impact of PRRSV on
  United States pork producers.  *J Swine Health Prod* 21(2):72-84.
- Segalés J, Allan GM, Domingo M (2005) Porcine circovirus diseases.
  *Anim Health Res Rev* 6(2):119-142.  doi:10.1079/AHR2005106
- Maes D et al. (2018) Swine enzootic pneumonia: a review.  *Vet Rec Open*
  5(1):e000228.  doi:10.1136/vetreco-2017-000228
- Harms PA et al. (2002) Experimental reproduction of PRDC with PRRSV, Mhp,
  and PCV2.  *J Swine Health Prod* 10(3):105-111.
- Zhang X et al. (2019) Loop-mediated isothermal amplification for Mhp
  detection.  *J Vet Diagn Invest* 31(4):571-576.
  doi:10.1177/1040638719843372
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63.  doi:10.1093/nar/28.12.e63

---

## Recipe 36 -- Swine neonatal enteric three-target panel (`lamp-forge swine-enteric-risk`)

**Goal.** A portable BioVind-style three-channel LAMP panel that simultaneously
detects the three dominant pathogens of neonatal piglet diarrhea -- Porcine
Epidemic Diarrhea Virus (PEDV), Porcine Deltacoronavirus (PDCoV), and Porcine
Rotavirus A -- and maps the result to a quarantine / biosecurity decision within
60 minutes, pen-side in a farrowing barn.

**Why it matters.** Neonatal diarrhea is the leading cause of pre-weaning
piglet mortality in commercial swine production worldwide.  PEDV alone caused
the death of over 10% of the US sow herd in the 2013-14 outbreak.  The three
pathogens are clinically indistinguishable (watery yellow diarrhea, rapid
dehydration, prostration), yet they require fundamentally different biosecurity
responses:

| Pathogen | Neonatal CFR | Key action |
|---|---|---|
| PEDV | ~100% (< 1 wk) | Pen quarantine, barn lockdown, stop all movement |
| PDCoV | 40-80% (< 1 wk) | Quarantine, intensive supportive care |
| Rotavirus A | 10-40% | Oral rehydration therapy, monitor for coronavirus co-infection |

A PCR lab turnaround of 24-48 h is too slow to prevent barn-to-barn fomite
spread.  A pen-side LAMP result in under 60 min drives the right intervention
before the attending herdsperson leaves the farrowing unit.

**Why this is a BioVind fit.** BioVind's farm-biosecurity vertical targets
exactly this scenario: a portable device, a multiplex LAMP panel, and an
actionable on-device clinical report.  All three pathogens are RNA-virus targets,
so each channel is an RT-LAMP assay fully compatible with BioVind BioID at
63-65 degC.

**Primer design -- one target per run, then combine.**

Recipe 27 already covers PEDV standalone design.  This recipe adds PDCoV and
Porcine Rotavirus A and combines all three into the neonatal enteric panel.

```bash
# Design each target independently (network required)

lamp-forge run --config config/pedv_N_gene.yaml
lamp-forge run --config config/pdcov_N_gene.yaml
lamp-forge run --config config/porcine_rota_a.yaml
```

**Cross-dimer compatibility screen.**

```bash
lamp-forge panel \
  --set PEDV=results/pedv_N_gene/primer_sets.json \
  --set PDCoV=results/pdcov_N_gene/primer_sets.json \
  --set RotaA=results/porcine_rota_a/primer_sets.json \
  --top-per-target 5 \
  --dimer-dg-threshold -5.0 \
  --out results/swine_enteric_3plex
```

**Verify RT-LAMP readiness for each channel.** All three are RNA targets:

```bash
lamp-forge rt-check \
  --input results/pedv_N_gene/primer_sets.json --na-type rna \
  --out-csv results/pedv_N_gene/rt_check.csv

lamp-forge rt-check \
  --input results/pdcov_N_gene/primer_sets.json --na-type rna \
  --out-csv results/pdcov_N_gene/rt_check.csv

lamp-forge rt-check \
  --input results/porcine_rota_a/primer_sets.json --na-type rna \
  --out-csv results/porcine_rota_a/rt_check.csv
```

All core primers (F3, B3, FIP, BIP) must be >= 63.0 degC.  Rotavirus A is
dsRNA; chaotropic lysis or 95 degC heat denaturation releases ssRNA template
that is efficiently primed by NEB RTx at 63-65 degC.

**Estimate LOD before ordering** (rectal swab in 500 uL PBS, 50% RNA
extraction, 50 uL eluate, 5 uL to reaction):

```bash
lamp-forge lod \
  --sample-volume 500 \
  --efficiency 0.50 \
  --eluate-volume 50 \
  --reaction-input 5
```

Effective sample = 500 x 0.5 x (5/50) = 25 uL per reaction.  LOD_95 approx
120 copies/mL of PBS swab eluate.  Faecal PEDV titres in acutely infected
neonatal piglets typically range from 10^5 to 10^9 copies/g (He et al. 2019),
providing orders-of-magnitude headroom above the assay LOD.

**Pool the working stock.** Three RNA targets, minimum stock concentration
3 x 44 = 132 uM; request 200 uM resuspension from IDT or Twist:

```bash
lamp-forge pool \
  --panel results/swine_enteric_3plex/panel.json \
  --stock-conc 200 \
  --total-volume 500 \
  --out results/swine_enteric_3plex/enteric_pool_sheet.csv
```

**Export to IDT for synthesis:**

```bash
lamp-forge export \
  --input results/pedv_N_gene/primer_sets.json \
  --format idt \
  --target-label PEDV_N \
  --out orders/pedv_idt_order.csv

lamp-forge export \
  --input results/pdcov_N_gene/primer_sets.json \
  --format idt \
  --target-label PDCoV_N \
  --out orders/pdcov_idt_order.csv

lamp-forge export \
  --input results/porcine_rota_a/primer_sets.json \
  --format idt \
  --target-label RotaA_VP6 \
  --out orders/rota_a_idt_order.csv
```

**Interpreting a pen-side result.** Feed the three positivity flags from the
device into the swine enteric risk engine:

```bash
# PEDV detected alone (critical)
lamp-forge swine-enteric-risk --pedv

# PDCoV and Rotavirus A co-detected (no PEDV)
lamp-forge swine-enteric-risk --pdcov --rota-a

# All three co-detected (worst-case neonatal enteric outbreak)
lamp-forge swine-enteric-risk --pedv --pdcov --rota-a

# Save for trend tracking
lamp-forge swine-enteric-risk --pedv \
  --out-json results/barn-A/farrowing/2026-01.json
```

**Scoring model and alert levels:**

| Detection pattern | Score | Alert | Key action |
|---|---|---|---|
| PEDV (any co-infection) | 50+ | CRITICAL | Pen quarantine; barn lockdown; stop all pig movement |
| PDCoV alone or + RotaA | 30-45 | HIGH | Quarantine unit; intensive supportive care; vet call |
| Rotavirus A alone | 15 | MODERATE | Oral rehydration therapy; monitor 24-48 h |
| Any single positive | 1-14 | LOW | Review with site veterinarian |
| All negative | 0 | NEGATIVE | Confirm extraction control positive |

`immediate_quarantine_required` is set True whenever PEDV is detected,
regardless of co-infection status.  `biosecurity_lockdown` is set True
for any porcine enteric coronavirus (PEDV or PDCoV).

### Tracking trend across monitoring intervals

Time-series trend analysis converts a sequence of panel results into a
trajectory (EMERGING / STABLE_CLEAR / STABLE_ENDEMIC / RESOLVING /
INSUFFICIENT_DATA) and flags the critical PEDV first-detection and
clearance events:

```bash
lamp-forge swine-enteric-trend \
  --enteric-result results/barn-A/farrowing/2026-01.json \
  --enteric-result results/barn-A/farrowing/2026-02.json \
  --enteric-result results/barn-A/farrowing/2026-03.json \
  --out-json results/barn-A/farrowing/trend_Q1.json \
  --out-csv results/barn-A/farrowing/trend_Q1.csv
```

Or from a monitoring CSV (header required):

```bash
lamp-forge swine-enteric-trend \
  --csv results/barn-A/farrowing/monitoring.csv \
  --out-json results/barn-A/farrowing/trend.json
```

CSV format:

```
sample_id,date,pedv,pdcov,rota_a
BarnA-2026-01,2026-01-15,1,0,0
BarnA-2026-02,2026-02-15,1,0,1
BarnA-2026-03,2026-03-15,0,0,1
```

Boolean columns accept `0`/`1` or `true`/`false` (case-insensitive).
The `date` column is optional.

Trend directions and their meaning:

| Direction | Meaning | Immediate action |
|---|---|---|
| EMERGING | PEDV newly detected or pathogen burden rising | Full outbreak protocol; movement restriction; vet call |
| STABLE_CLEAR | All targets negative across intervals | Routine monitoring; continue surveillance schedule |
| STABLE_ENDEMIC | Rotavirus A or PDCoV endemic but no PEDV | Standard management; review vaccination; no movement restriction |
| RESOLVING | Pathogen burden declining after a detection event | Maintain restrictions until 3 consecutive negatives |
| INSUFFICIENT_DATA | Fewer than 2 intervals | Collect additional samples before trend classification |

PEDV first-detection (`pedv_newly_detected`) and clearance (`pedv_cleared`)
events are flagged explicitly in the JSON output.

### Wet-lab notes

- **Sample type.** Rectal swab in 500 uL PBS, or intestinal content from
  acutely affected piglets.  Intestinal content from the small intestine
  (jejunum / ileum) yields the highest viral titres for PEDV and PDCoV.
  Swabs from multiple litters should be pooled before extraction to increase
  surveillance sensitivity.
- **Lysis.** Use a chaotropic lysis buffer (guanidinium thiocyanate) for
  one-step extraction; this simultaneously lyses enveloped coronaviruses
  and denatures Rotavirus A dsRNA.  A 5-minute 95 degC incubation step
  is optional for Rotavirus A with column-free extraction protocols but is
  not required with NEB WarmStart RTx at 63-65 degC when combined with
  guanidinium lysis.
- **Internal control.** Include a 16S rRNA LAMP channel (Recipe 19) or
  an RNA spike-in (MS2 phage) to confirm extraction success.  A negative
  internal control invalidates the run -- do not report a clean result from a
  panel where the extraction control is also negative.
- **Co-infection.** All three pathogens can co-circulate in a single piglet.
  The scoring model is additive (PEDV + RotaA = score 65, CRITICAL).  Confirm
  PCR subtype information at an accredited laboratory if regulatory notification
  or vaccine efficacy assessment is required.

**References.**

- Stevenson GW, Hoang H, Schwartz KJ et al. (2013) Emergence of porcine
  epidemic diarrhea virus in the United States: clinical signs, lesions, and
  viral genomic sequences.  *J Vet Diagn Invest* 25(5):649-654.
  doi:10.1177/1040638713501675
- Song D, Zhou X, Peng Q et al. (2015) Newly emerged porcine deltacoronavirus
  associated with diarrhoea in swine in China: identification, prevalence and
  full-length genome sequence analysis.  *Transbound Emerg Dis* 62(6):575-580.
  doi:10.1111/tbed.12399
- Saif LJ, Pensaert MB, Sestak K, Yeo S-G, Jung K (2012) Coronaviruses.
  In: Zimmerman JJ et al. (eds) *Diseases of Swine*, 10th ed.
  Wiley-Blackwell, pp 501-524.
- Khamrin P, Maneekarn N, Hidaka S et al. (2010) Loop-mediated isothermal
  amplification for rapid detection of rotavirus in clinical stool specimens.
  *J Virol Methods* 163:384-389.  doi:10.1016/j.jviromet.2009.11.013
- He X, Guo S, Wang Y et al. (2019) Development and evaluation of a rapid and
  sensitive reverse-transcription loop-mediated isothermal amplification method
  for the detection of porcine epidemic diarrhea virus.
  *Transbound Emerg Dis* 66(1):431-439.  doi:10.1111/tbed.13044
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63.  doi:10.1093/nar/28.12.e63

## Recipe 37 -- POC urinary tract infection four-target panel (`lamp-forge uti-risk`)

**Goal.** A portable BioVind-style four-channel LAMP panel that simultaneously
detects the four most common community-acquired UTI uropathogens -- *Escherichia
coli*, *Klebsiella pneumoniae*, *Enterococcus faecalis*, and *Staphylococcus
saprophyticus* -- and maps the result to an organism-directed prescribing
recommendation within 60 minutes, at the point of care.

**Why it matters.** Urinary tract infection is the most common bacterial
infection in primary-care and urgent-care settings (~7 million outpatient
visits per year in the United States alone).  Empiric antibiotic prescribing
is the standard of care because urine cultures take 24-48 hours.  Yet resistance
to first-line agents (fluoroquinolones, trimethoprim-sulfamethoxazole) now
exceeds 20% for *E. coli* in many locales, and ESBL-producing *K. pneumoniae*
is intrinsically resistant to cephalosporins.  The right antibiotic depends on
the organism, and a wrong-class prescription accelerates resistance while failing
the patient.

A same-visit LAMP result resolves the two questions that drive prescribing:
**what is the organism?** and **is there a resistance concern?**

| Target | Gene | Community UTI prevalence | Key antibiotic note |
|---|---|---|---|
| *E. coli* | uidA | 75-80% | Nitrofurantoin / TMP-SMX first-line; ESBL strains require fosfomycin or carbapenem |
| *K. pneumoniae* | phoE | 5-10% | ESBL/KPC risk; send culture before prescribing; avoid cephalosporins |
| *E. faecalis* | ddl | 5-10% | Intrinsic cephalosporin resistance; use ampicillin or nitrofurantoin |
| *S. saprophyticus* | sdrF | 5-15% young women | Nitrofurantoin / trimethoprim; low-count significant |

**Why this is a BioVind fit.** BioVind's human point-of-care vertical targets
rural clinics, urgent-care centres, and telemedicine-adjacent collection sites
where same-visit results change prescribing behaviour.  All four targets are
DNA-LAMP (no RT step required), making the panel technically simpler than the
respiratory viral panels and compatible with any isothermal LAMP device at
63-65 degC.

**Primer design -- one target per run, then combine.**

```bash
# E. coli uidA (beta-D-glucuronidase -- E. coli-specific, absent in Salmonella/Shigella)
lamp-forge run --config config/ecoli_uidA.yaml

# K. pneumoniae phoE (phosphate outer-membrane porin, chromosomally encoded)
lamp-forge run --config config/kpneu_phoE.yaml

# E. faecalis ddl (D-Ala:D-Ala ligase, species-specific variant)
lamp-forge run --config config/efaec_ddl.yaml

# S. saprophyticus sdrF (serine-aspartate repeat protein F)
lamp-forge run --config config/ssaph_sdrF.yaml
```

**Config template for *E. coli* uidA.**

```yaml
target:
  name: ecoli_uidA_UTI
  taxon_id: 562                # Escherichia coli
  gene: uidA
  max_sequences: 30
  email: you@your-inst.edu

off_targets:
  fasta_dir: input/off_targets_uti
  min_identity_threshold: 0.80
  min_coverage_threshold: 0.80

conservation:
  window_size: 30
  entropy_threshold: 0.20
  min_region_length: 200

primer:
  tm_min: 60.0
  tm_max: 65.0
  gc_min: 40.0
  gc_max: 65.0
  hairpin_dg_threshold: -2.0
  dimer_dg_threshold: -5.0
  amplicon_size:
    f2_b2_min: 120
    f2_b2_max: 160

output:
  dir: results/ecoli_uidA
  top_n: 5
  generate_html: true
  generate_csv: true
```

**Off-target panel for UTI (`input/off_targets_uti/`).**

| File | NCBI accession | Why |
|---|---|---|
| `kpneu_ref.fasta` | NC_009648.1 | Cross-species screen for uidA assay |
| `ecoli_ref.fasta` | NC_000913.3 | Cross-species screen for phoE assay |
| `efaecium_ref.fasta` | NC_017960.1 | Distinguish *E. faecalis* from *E. faecium* |
| `efaecalis_ref.fasta` | NC_004668.1 | Cross-species screen for sdrF assay |
| `saureus_ref.fasta` | NC_002951.2 | Confirm sdrF is *S. saprophyticus*-specific |
| `human_urothelium_fragment.fasta` | GRCh38 chr17 region | Host DNA from urine epithelial cells |

**Cross-dimer compatibility screen.**

```bash
lamp-forge panel \
  --set Ecoli=results/ecoli_uidA/primer_sets.json \
  --set Kpneu=results/kpneu_phoE/primer_sets.json \
  --set Efaec=results/efaec_ddl/primer_sets.json \
  --set Ssaph=results/ssaph_sdrF/primer_sets.json \
  --dimer-threshold -5.0 \
  --out-json results/uti_panel/panel.json \
  --out-csv  results/uti_panel/panel_summary.csv
```

**Equimolar pool for the working-stock mix.**

All four targets are DNA-LAMP; use standard pooling stoichiometry
(FIP/BIP at 16 uM, LF/LB at 4 uM, F3/B3 at 2 uM per target in the 10x stock):

```bash
lamp-forge pool \
  --panel-json results/uti_panel/panel.json \
  --stock-conc 200 \
  --total-volume 200 \
  --out-csv results/uti_panel/pool_sheet.csv
```

**Interpreting a panel result with `lamp-forge uti-risk`.**

Once the device returns positivity flags from a patient urine specimen, feed
the flags into the risk assessment engine to get a structured prescribing alert:

```bash
# E. coli-positive uncomplicated UTI
lamp-forge uti-risk --ecoli

# K. pneumoniae detected -- culture urgently needed
lamp-forge uti-risk --kpneu --out-json results/uti/kpneu_alert.json

# Polymicrobial E. coli + E. faecalis
lamp-forge uti-risk --ecoli --efaec --out-csv results/uti/mixed_uti.csv

# S. saprophyticus -- young woman, uncomplicated
lamp-forge uti-risk --ssaph

# All negative -- send culture
lamp-forge uti-risk
```

**Example output for K. pneumoniae-positive UTI.**

```
POC UTI LAMP Panel Assessment
========================================
Alert level : HIGH  (score 35/100)
Detected    : K.pneu

  Gram-negative uropathogen detected -- urine culture recommended.
  K. pneumoniae: ESBL/KPC resistance risk -- avoid cephalosporins empirically.
  Antibiotic guidance: organism identity changes first-line selection.

Interpretation: Klebsiella pneumoniae detected via the phoE assay. K. pneumoniae
accounts for 5-10% of community UTI and carries a significant risk of ESBL or,
in healthcare-adjacent patients, KPC-type carbapenemase production.

Recommended action: HIGH: Send urine culture with susceptibility testing before
committing to empiric therapy. If treatment cannot be deferred, nitrofurantoin
or fosfomycin are reasonable empiric choices for uncomplicated cystitis while
awaiting culture results ...
```

**Alert level reference.**

| Alert level | Score | Meaning | First-line guidance |
|---|---|---|---|
| HIGH | >= 35 | Gram-negative uropathogen (*E. coli* or *K. pneumoniae*) | Nitrofurantoin or TMP-SMX (if susceptible); send culture for *K. pneumoniae* |
| MODERATE | >= 20 | *E. faecalis* alone | Ampicillin or nitrofurantoin; rule out VRE by culture |
| LOW | >= 10 | *S. saprophyticus* alone | Nitrofurantoin 5-7 days or trimethoprim 7 days |
| NEGATIVE | 0 | All targets negative | Send conventional urine culture if clinical suspicion persists |

**LOD estimation before ordering primers.**

```bash
# Typical urine: 200 uL input, 85% extraction, 50 uL eluate, 5 uL reaction input
lamp-forge lod \
  --sample-volume 200 \
  --extraction-efficiency 0.85 \
  --eluate-volume 50 \
  --reaction-input 5 \
  --confidence 0.95
```

Expect a reportable LOD around 3-6 CFU/mL -- well below the 10^3 CFU/mL
threshold for *S. saprophyticus* significance and comfortably within the
10^5 CFU/mL standard for gram-negative UTI.

**Validation plan for regulatory filing.**

BioVind acceptance criteria: Sn >= 95%, Sp >= 99%.

```bash
lamp-forge validation-plan \
  --sn-min 0.95 --sp-min 0.99 \
  --p-expected-pos 0.97 --p-expected-neg 1.00 \
  --power 0.80 \
  --out-csv results/uti_panel/validation_plan.csv
```

**References.**

- Foxman B (2002) Epidemiology of urinary tract infections: incidence,
  morbidity and economic costs. *Am J Med* 113(Suppl 1A):5S-13S.
  doi:10.1016/S0002-9343(02)01054-9
- Gupta K, Hooton TM, Naber KG et al. (2011) International clinical practice
  guidelines for the treatment of acute uncomplicated cystitis and pyelonephritis
  in women. *Clin Infect Dis* 52(5):e103-e120. doi:10.1093/cid/ciq257
- Terlizzi ME, Gribaudo G, Maffei ME (2017) UroPathogenic *Escherichia coli*
  (UPEC) infections: virulence factors, bladder responses, antibiotic, and
  non-antibiotic antimicrobial strategies. *Front Microbiol* 8:1566.
  doi:10.3389/fmicb.2017.01566
- Peirano G, Pitout JD (2019) Extended-spectrum beta-lactamase-producing
  *Enterobacteriaceae*: update on molecular epidemiology and treatment options.
  *Drugs* 79(14):1529-1541. doi:10.1007/s40265-019-01180-3
- Hedman P, Ringertz O, Lindstrom M, Olsson K (1991) The origin of
  *Staphylococcus saprophyticus* from cattle and pigs. *Scand J Infect Dis*
  23(6):769-773. doi:10.3109/00365549109024318
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63

---

## Recipe 38 -- Bovine mastitis herd surveillance trajectory (`lamp-forge mastitis-trend`)

**Goal.** Track the bovine mastitis situation in a dairy herd across consecutive
LAMP panel results to determine whether contagious pathogens (*Staphylococcus aureus*,
*Streptococcus agalactiae*) are newly emerging, persisting, or being eradicated --
and whether environmental mastitis (*Streptococcus uberis*, *E. coli*) is a recurring
herd-level problem requiring husbandry correction.

**Why it matters.** A single LAMP result tells you what pathogens are present in
one milk sample right now.  Trend analysis across monthly or quarterly surveillance
intervals tells you whether your management programme is working.  The decision to
cull a high-SCC *S. aureus* cow, declare a herd free of *S. agalactiae*, or commission
a cubicle hygiene review all require a time series -- one result is not enough.

`lamp-forge mastitis-trend` converts a sequence of `mastitis-risk` results into one
of six trend directions:

| Direction | Meaning | Immediate action |
|---|---|---|
| NEWLY_DETECTED | Contagious pathogen appears for first time | Segregate cow; whole-herd LAMP screen within 48 h |
| PERSISTENT | Contagious pathogen in current AND prior sample(s) | Whole-herd screen to find carrier cow; vet protocol |
| RESOLVING | Contagious pathogen was present, now absent | Maintain surveillance; confirm 3 consecutive negatives |
| STABLE_CLEAR | All samples negative | Routine monitoring; continue teat disinfection |
| STABLE_ENDEMIC | No contagious pathogen; environmental recurring | Review cubicle hygiene, bedding, dry-cow protocol |
| INSUFFICIENT_DATA | Fewer than 2 samples | Collect next monitoring sample first |

**Typical workflow -- quarterly contagious mastitis surveillance.**

```bash
# Month 1 -- S. aureus detected at milking
lamp-forge mastitis-risk --saur \
  --out-json results/herd-A/2026-01.json

# Month 2 -- S. aureus + S. uberis co-detected
lamp-forge mastitis-risk --saur --sube \
  --out-json results/herd-A/2026-02.json

# Month 3 -- all negative (eradication programme completed)
lamp-forge mastitis-risk \
  --out-json results/herd-A/2026-03.json

# Analyse the trajectory
lamp-forge mastitis-trend \
  --mastitis-result results/herd-A/2026-01.json \
  --mastitis-result results/herd-A/2026-02.json \
  --mastitis-result results/herd-A/2026-03.json \
  --out-json results/herd-A/trend_Q1.json \
  --out-csv results/herd-A/trend_Q1.csv
```

Expected output (month 3 negative after two positives):

```
Bovine mastitis surveillance trajectory: 3 sample(s)

Sample ID     Date          SAUR  SAGA  SUBE  ECOLI  Alert
---------------------------------------------------------------
2026-01       2026-01-15     [+]   [-]   [-]   [-]   HIGH
2026-02       2026-02-15     [+]   [-]   [+]   [-]   HIGH
2026-03       2026-03-15     [-]   [-]   [-]   [-]   NEGATIVE

Trend direction: RESOLVING
  Contagious pathogen cleared: confirm with 3 consecutive negatives
    before ending enhanced surveillance.
  Worst alert level: HIGH
```

**CSV monitoring-sheet input.** For herds where results are recorded in a
spreadsheet rather than individual JSON files:

```bash
lamp-forge mastitis-trend \
  --csv monitoring/herd_A_2026.csv \
  --out-json results/herd-A/trend_2026.json
```

CSV format (header required):

```
sample_id,date,saur,saga,sube,ecoli
HerdA-2026-01,2026-01-15,1,0,0,0
HerdA-2026-02,2026-02-15,1,0,1,0
HerdA-2026-03,2026-03-15,0,0,0,0
```

Boolean columns accept `0`/`1` or `true`/`false`.  The `date`, `sube`,
and `ecoli` columns are optional.

**Persistent contagious mastitis -- the carrier-cow problem.**  When the same
contagious pathogen is detected in >= 75 % of >= 3 monitoring intervals,
`persistent_contagious_likely` is flagged and a highlighted warning appears in
the CLI output.  This is the signal to escalate to whole-herd individual-cow
LAMP screening:

```bash
lamp-forge mastitis-risk --saur --out-json results/herd-B/2026-01.json
lamp-forge mastitis-risk --saur --out-json results/herd-B/2026-02.json
lamp-forge mastitis-risk --saur --sube --out-json results/herd-B/2026-03.json

lamp-forge mastitis-trend \
  --mastitis-result results/herd-B/2026-01.json \
  --mastitis-result results/herd-B/2026-02.json \
  --mastitis-result results/herd-B/2026-03.json
# Output: PERSISTENT -- persistent_contagious_likely=True
```

**Environmental mastitis monitoring -- STABLE_ENDEMIC.**  When there is no
contagious pathogen across all samples but *S. uberis* or *E. coli* appear in
two or more intervals, the direction is `STABLE_ENDEMIC`:

```bash
lamp-forge mastitis-risk --sube --out-json results/herd-C/2026-01.json
lamp-forge mastitis-risk --sube --ecoli --out-json results/herd-C/2026-02.json
lamp-forge mastitis-risk --sube --out-json results/herd-C/2026-03.json

lamp-forge mastitis-trend \
  --mastitis-result results/herd-C/2026-01.json \
  --mastitis-result results/herd-C/2026-02.json \
  --mastitis-result results/herd-C/2026-03.json
# Output: STABLE_ENDEMIC
# Recommended action: Review cubicle and bedding management...
```

**References.**

- Hogeveen H, Huijps K, Lam TJGM (2011) Economic aspects of mastitis: new
  developments. *N Z Vet J* 59(1):16-23. doi:10.1080/00480169.2011.552541
- Bradley AJ, Green MJ (2004) The importance and control of mastitis in
  dairy cattle. *Vet Clin N Am Food Anim Pract* 20(3):469-491.
  doi:10.1016/j.cvfa.2004.06.010
- Keefe G (2012) Update on control of *Staphylococcus aureus* and
  *Streptococcus agalactiae* for management of mastitis. *Vet Clin N Am
  Food Anim Pract* 28(2):203-216. doi:10.1016/j.cvfa.2012.03.010
- Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Res* 28(12):e63. doi:10.1093/nar/28.12.e63
