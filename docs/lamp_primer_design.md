# LAMP primer geometry — a designer's reference

This is a concise reference for the geometric and thermodynamic constraints
that LAMP-Forge enforces. The canonical source is **Notomi et al. (2000)
"Loop-mediated isothermal amplification of DNA," *Nucleic Acids Research*
28(12):E63**, supplemented by Eiken Chemical's *A Guide to LAMP Primer
Design* (the company's design recommendations are open and widely cited).

## The eight binding sites

A complete LAMP target is divided into eight regions across both strands:

```
                  F3c               F2c              F1c
   ──────────────|─────|─────|─────|─────|─────|─────|─────|──────────────
                                                                            template (5'→3')
   ──────────────|─────|─────|─────|─────|─────|─────|─────|──────────────
                                                                            template (3'→5')
                  B1                B2               B3
                          F1c arm                  B1c arm
                          (in FIP)                 (in BIP)
```

(Equivalent diagrams in Notomi 2000 Fig. 1; the convention is that *c*
suffixes mean "complement of," so F1c is the reverse-complement of F1.)

The six designed primers and what they hit:

| Primer | Anneals to | Strand of template | Role |
|---|---|---|---|
| **F3** | F3c | sense | Outer displacement (forward) |
| **B3** | B3c | sense | Outer displacement (backward) |
| **FIP** | F2c (with F1c arm dangling) | sense | Inner forward, chimeric |
| **BIP** | B2c (with B1c arm dangling) | sense | Inner backward, chimeric |
| **LF** | F2c–F1c gap region (loop) | sense | Loop-accelerator, forward |
| **LB** | B1c–B2c gap region (loop) | sense | Loop-accelerator, backward |

## Geometric constants (the numbers we enforce)

```
5' ─F3─┬─F2─┬───────F1───────[loop]───────B1c──────┬─B2─┬─B3─ 3'
       │    │                                       │    │
       │ ←0-20→                                  ←0-20→  │
       │      ←F2-F1 = 40-60 bp→                       │
       │                       ←F2-B2 = 120-160 bp→    │
       └────────────────── F3-B3 amplicon ─────────────┘
```

| Span | Default range | Notes |
|---|---|---|
| F3 length | 18-22 nt | Tm-driven within band |
| F2 length | 18-22 nt | Tm-driven; identical band to F3 |
| F1c length | 20-22 nt | Slightly longer — anchors the chimera |
| B1c length | 20-22 nt | Same as F1c |
| B2 length | 18-22 nt | Same as F2 |
| B3 length | 18-22 nt | Same as F3 |
| Loop primer length | 18-22 nt | Sits in unpaired loop |
| F3 → F2 gap | 0-20 bp | "Just clear" displacement spacing |
| B2 → B3 gap | 0-20 bp | Mirror of above |
| F2 → F1 gap | 40-60 bp | Loop size driver — tighter is worse |
| B1 → B2 gap | 40-60 bp | Same |
| F2 5' → B2 5' | 120-160 bp | Inner amplicon; affects reaction kinetics |
| F3 5' → B3 5' | <250 bp | Full amplicon; longer slows the reaction |

## Thermodynamic constraints

| Parameter | Default | Reasoning |
|---|---|---|
| Tm | 60-65 °C | Bst polymerase optimum |
| Tm spread across set | ≤ 2 °C | Mismatched Tms desync the reaction |
| GC content | 40-65% | Lower → weak binding; higher → secondary structure |
| Hairpin ΔG | ≥ -2 kcal/mol | Below this, intramolecular folding competes |
| Homodimer ΔG | ≥ -5 kcal/mol | Below this, dimers reduce active primer concentration |
| Heterodimer ΔG (across pairs) | ≥ -5 kcal/mol | Especially watch FIP/BIP × loop primers |

Tm is calculated with primer3's salt-corrected nearest-neighbour
thermodynamic tables (Allawi & SantaLucia 1997 parameters). The salt
defaults we use match typical LAMP buffer conditions:

- Monovalent (Na+ equivalents): 50 mM
- Divalent (Mg2+): 8 mM
- dNTPs: 0.8 mM
- Primer concentration: 0.2 µM (effective concentration after dNTP/Mg
  corrections)

If you're working with a non-standard buffer (e.g. some isothermal master
mixes use 6 mM Mg2+ for stringency), adjust `lamp_forge.primer_design.tm()`
arguments before running.

## Common pitfalls

**1. The F2 5' end matters more than the F2 3' end.** During the reaction,
F2 primes from its 3' end but the loop architecture is templated by F2 + F1c
spacing measured from F2's 5' end. LAMP-Forge places F2 coordinates as
5'-on-template; if you read primers off a synthesis-house dropdown that
reports the 3'-on-template position, translate carefully.

**2. FIP is written 5' → 3' as F1c-F2, not F2-F1c.** F1c sits first (it's the
"looped-back" arm) and F2 sits second (it's the actual annealing arm). The
ordering matters for synthesis. LAMP-Forge's CSV output uses the standard
F1c-F2 ordering.

**3. Loop primers occasionally interfere with FIP/BIP.** The loop primers sit
in the F2-F1 (or B1-B2) gap. If they overlap with the FIP/BIP annealing
region at all, expect dimer or stalling. LAMP-Forge enforces non-overlap by
design — loop primer scan windows are bounded by the inner-side coordinates
of F2 and F1c (mirror for B side).

**4. AT-rich regions amplify weakly.** Bst polymerase has slightly worse
fidelity and slightly lower activity in low-GC contexts. If your conserved
region is <40% GC, raise `tm_min` to 62 °C and accept that you may get fewer
candidate sets.

## Recommended downstream validation

Even after a clean LAMP-Forge run, **always**:

1. **Run a temperature gradient (60-66 °C, 1 °C steps) on synthetic
   double-stranded target DNA** at 10⁵ copies/reaction. The optimum will
   often not be exactly 65 °C.
2. **Test specificity on your real off-target panel** at 10⁶ copies, with
   the optimum incubation time + 50%. False positives appear in the tail.
3. **Sensitivity ladder:** 10⁶, 10⁵, ..., 10⁰ copies / reaction, in
   triplicate. LAMP usually limits at 10-100 copies for a well-designed
   set. If you bottom out at 10³, the set is poorly performing.
4. **Inclusivity panel:** confirm amplification on clinical isolates spanning
   your target's known diversity. The conservation step picks regions that
   should work on your training set; clinical reality always has more variation.

For diagnostic deployment, also add: sample-matrix testing (sputum, stool,
buffer), an internal amplification control (a co-amplified positive control
sequence), and a contrived clinical panel before patient samples.

## Further reading

- Notomi T et al. (2000). Loop-mediated isothermal amplification of DNA.
  *Nucleic Acids Research* 28:E63.
- Nagamine K, Hase T, Notomi T (2002). Accelerated reaction by loop-mediated
  isothermal amplification using loop primers. *Molecular and Cellular
  Probes* 16:223-229. **(The LF/LB paper.)**
- Tanner NA, Evans TC Jr (2014). Loop-Mediated Isothermal Amplification for
  Detection of Nucleic Acids. *Current Protocols in Molecular Biology*
  105:15.14.1-15.14.14.
- Becherer L et al. (2020). Loop-mediated isothermal amplification (LAMP) —
  review and classification of methods for sequence-specific detection.
  *Analytical Methods* 12:717-746.
