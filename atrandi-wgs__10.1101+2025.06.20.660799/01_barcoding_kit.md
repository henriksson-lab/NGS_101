# Atrandi Single-Microbe DNA Barcoding Kit — protocol quick reference

Source: `DGPM02323198001_Single_Microbe_DNA_Barcoding_Kit_V3.pdf`
Doc. No. **DGPM02323198001**, Revision **V3**, 03 November 2023 (V1 2 Oct 2023; V2 13 Oct 2023).
Atrandi Biosciences / Droplet Genomics. Research Use Only.
Kit P/N **CKP-BARK1**. 3 experiments, ~10,000 barcoded cells sequenced per experiment.

> **Evidence marking.** 🟢 = verbatim from the source document · 🟡 = derived or inferred
> (reasoning given; arithmetic verified where possible) · 🔴 = not published anywhere.
> Any claim not marked 🟢 must reach the final page carrying its uncertainty.

> This document ends at barcoding. DNA purification and Illumina library prep are in the
> companion guide DGPM02323206001 — see [`02_library_prep.md`](02_library_prep.md).

---

## 1. Principle

Semi-permeable capsules (SPCs): picolitre aqueous core inside a hydrogel shell. The shell is a
size-selective membrane — small reagents (enzymes, detergents, salts) diffuse in and out freely,
while genomic DNA and cells are retained. This is what allows a harsh alkaline lysis and multiple
buffer exchanges without losing single-cell compartmentalisation, and it is what makes split-pool
barcoding possible without any single-cell isolation.

Workflow: encapsulation → fixation → lysis → whole-genome amplification → debranching →
end-prep → 4 rounds of split-pool barcode ligation → (companion doc) release, library prep, sequencing.

Total procedure time ≈ 8 h 25 min (3 h 35 min hands-on, 4 h 50 min incubation).

## 2. Barcode architecture

- 4 barcode libraries: **A, B, C, D**, each **24 variants**, each **8 nt**.
- Diversity: 24 × 24 × 24 × 24 = **331,776** combinations.
- One 96-well plate, split by column — all 8 rows (A–H) used in each block:

| Barcode | Plate columns | n |
|---|---|---|
| A | 1–3 | 24 |
| B | 4–6 | 24 |
| C | 7–9 | 24 |
| D | 10–12 | 24 |

- Ligation order is **A → B → C → D**. 🟡 D is ligated last, so D ends up outermost / proximal to
  the sequencing adapter, and is read **first** in Read 2. (Confirmed independently by the read
  layout in `03_our_protocol.md` §8.)
- Final construct (Figure 3):

```
5'- P7 - Index i7 - Adapter - D - C - B - A - {Fragment} - Adapter - Index i5* - P5 -3'
                             \_______________/
                                cell barcode              (*i5 optional)
```

- Spacing in the read is **8 + 4 + 8 + 4 + 8 + 4 + 8 (+4)** — 8 nt barcode parts separated by
  **4 nt linkers**. See `02_library_prep.md` for read offsets.
- Barcode sequences (A–D, 24 each) are tabulated in the **companion** document, Tables 2 and 3.
  **Linker sequences are not given in either document.**

## 3. Kit contents

| Group | Component | P/N | Amount |
|---|---|---|---|
| Cell lysis | 2X Lysis Buffer (LB)* | CRP-LB1 | 3 × 0.4 mL |
| | Neutralization Buffer (NB) | CRP-NB1 | 3 × 5 mL |
| WGA | dNTP Mix (10 mM each) | CRP-DNTP1 | 140 µL |
| | Primer Mix | CRP-PRM1 | 70 µL |
| | WGA Polymerase | CRP-WGP1 | 70 µL |
| | WGA Enhancer | CRP-WGE1 | 25 µL |
| | 10X WGA Reaction Buffer | CRP-WGB1 | 140 µL |
| Debranching / end-prep | Debranching Enzyme | CRP-DBE1 | 50 µL |
| | 10X Debranching Buffer | CRP-DBB1 | 100 µL |
| | End Prep Enzyme Mix | CRP-EPE1 | 30 µL |
| | End Prep Reaction Buffer | CRP-EPB1 | 70 µL |
| Barcoding | Ligation Enzyme | CRP-LGE1 | 250 µL |
| | 10X Ligation Buffer | CRP-LGB1 | 750 µL |
| | Barcodes plate | CRP-PLT1 | 3 plates |
| | Ligation Adapter | CRP-LGA1 | 10 µL |
| Buffers | 50X Wash Buffer (WB)** | CRP-WB1 | 3 × 1 mL |
| | 10X STOP Buffer (SB)** | CRP-SB1 | 3 × 1 mL |
| | Water, nuclease-free (NFW) | CRP-WTR1 | 3 × 1 mL |

Storage −20 ± 5 °C.
\* add 100 µL of 1 M DTT to the tube before use. \*\* dilute to 1X in NFW before use.

**Not provided:** 100% methanol (fixation); 1 M DTT (lysis); 10% Pluronic F-68 / Poloxamer 188
(Thermo 24040032, for WGA); 0.1 M DTT (WGA); 1.5 mL non-stick tubes (e.g. Eppendorf LoBind);
SYTO™9 (Thermo S34854, SPC imaging); adhesive plate seals.

Note: the **Ligation Adapter (CRP-LGA1)** ships in *this* kit but is only used in the companion
library-prep protocol.

## 4. Protocol

### Starting material
Encapsulate on the **FLUX SPC Generator** or **ONYX** system. After emulsion breaking and washing,
obtain ≥150 µL packed SPCs; working range 150–500 µL.
**Start with ≥30,000 encapsulated cells to recover 10,000 sequenced genomes.**
Target occupancy ~0.1 (over-loading causes barcode cross-contamination).

### Cell fixation (steps 1–3)
1. Chill 100% methanol on ice or −20 °C.
2. To 100 µL SPCs add 900 µL cold methanol **drop by drop** while agitating.
3. Store at −20 °C / −80 °C for long-term storage.

### Preparation (steps 4–5)
Dilute 50X WB and 10X SB to 1X in NFW.

### Cell lysis (steps 6–13)
*Skip 6–7 if already fixed after encapsulation.*
6–7. Methanol fixation as above; incubate 30 min (or overnight) at ≤ −20 °C.
8. Wash **3×** with 1 mL 1X WB (vortex, 1 min @ 1000 × g, aspirate).
9. Bring volume to 500 µL with 1X WB.
10. Prepare fresh alkaline 2X LB by adding 100 µL of 1 M DTT; add 500 µL onto the SPCs, mix.
11. **Rotate 15 min at RT (20–25 °C).** Proceed immediately.
12. Wash **5×** with 1 mL Neutralization Buffer (NB).
13. Wash **5×** with 1 mL 1X WB.

> Alkaline lysis + DTT. **Do not use mechanical lysis with SPCs.**

### Whole genome amplification (steps 14–20)

| Reagent | µL |
|---|---|
| SPCs | 150 |
| NFW | 97.5 |
| 10X WGA Reaction Buffer | 37.5 |
| dNTP Mix | 37.5 |
| Primer mix | 18.75 |
| 0.1 M DTT | 3.75 |
| 10% Pluronic F-68 | 3.75 |
| WGA Enhancer | 7.5 |
| WGA Polymerase | 18.75 |
| **Total** | **375** |

Aliquot into 5 PCR tubes, 75 µL each.

| Temp | Time |
|---|---|
| 45 °C | 15 min |
| 65 °C | 10 min |
| 4 °C | hold |

17. Pool into 1.5 mL tube; rinse tubes with 100 µL 1X WB; pellet 1 min @ 1000 × g.
18. Wash **3×** with 1 mL WB.
19. *Recommended:* image with 10× SYTO9 to check occupancy of amplified genomes (Figure 9).
20. **Stop point:** 4 °C in Stop Buffer up to 7 days; wash 3× with WB before continuing.

> 🟡 **INFERRED.** The 45 °C / 15 min + 65 °C profile is consistent with isothermal
> strand-displacement (MDA-like) amplification. 🔴 Enzyme identity is proprietary — Atrandi state
> nothing beyond "WGA Polymerase".

### DNA debranching (steps 21–24)

| Reagent | µL |
|---|---|
| SPCs | 150 |
| 10X Debranching Buffer | 30 |
| Debranching Enzyme | 15 |
| NFW | 105 |
| **Total** | **300** |

22. **37 °C, 1 h, 1000 rpm** on a thermomixer.
23. Wash 3× with 1 mL 1X WB.
24. **Stop point:** 4 °C in WB up to 2 days.

> 🟡 **INFERRED.** Presumed to resolve the hyperbranched / cruciform structures produced by
> strand-displacement WGA into ligatable linear duplex ends. 🔴 The enzyme is not named; a
> structure-specific resolvase such as T7 endonuclease I would fit, but Atrandi do not say.

### DNA end preparation (steps 25–30)

**Barcode ligation must follow immediately.**

| Reagent | µL |
|---|---|
| SPCs | 150 |
| End Prep Reaction Buffer | 21 |
| End Prep Enzyme Mix | 9 |
| **Total** | **180** |

Aliquot into 2 PCR tubes, 90 µL each.

| Temp | Time |
|---|---|
| 20 °C | 30 min |
| 65 °C | 30 min |
| 4 °C | hold |

28. Pool, rinse with 100 µL 1X WB, pellet 1 min @ 1000 × g.
29. Wash 3× with 1 mL 1X WB.
30. Proceed immediately to step 31.

> 🟡 **INFERRED.** 20 °C = end repair/blunting, 65 °C = dA-tailing + inactivation — this is the
> standard NEB-style end-prep profile. Implies **3'-dA ends**, hence a T-overhang on the round-A
> barcode cassette. This inference carries real weight: see `05_barcode_cassette_model.md` §2.

### Split-and-pool barcoding (steps 31–44)

Thaw the barcode plate at RT, spin 1000 × g for 30 s. Cut the seal only over the wells in use.

**Barcode A ligation (steps 31–40)** — master mix:

| Reagent | µL |
|---|---|
| SPCs | 150 |
| 10X Ligation Buffer | 60 |
| Ligation Enzyme | 20 |
| NFW | 70 |
| **Total** | **300** |

32. Distribute to 8-strip tubes / 8 wells, 36 µL each (or directly to 24 wells).
33. 8-channel pipette **10 µL master mix into each well containing 10 µL barcodes A** (columns 1–3); mix.
34. Seal with a fresh adhesive seal.
35. **22 °C (or RT) for 15 min.**
36. Spin briefly.
37. Add **40 µL 1X Stop Buffer** per well.
38. 5 min at RT.
39. Pool into a 1.5 mL tube, pellet 1000 × g 1 min, remove supernatant.
    *Optional:* rinse wells with 40 µL SB to maximise recovery.
40. Wash **5×** with 1 mL 1X WB (1 min @ 1000 × g).

41. **Barcode B** — repeat 31–40 with columns 4–6.
42. **Barcode C** — repeat 31–40 with columns 7–9.
43. **Barcode D** — repeat 31–40 with columns 10–12.
44. **Stop point:** barcoded sample keeps at 4 °C in 1X Stop Buffer for a few weeks.

> Each round: 24-way split, 15 min ligation, stop, pool, wash. 🟡 Reaction volume per well is
> 20 µL (10 µL mix + 10 µL barcode), quenched with 40 µL SB — arithmetic from the stated volumes.

## 5. Troubleshooting (Table 1)

| Issue | Cause | Solution |
|---|---|---|
| No MDA-positive capsules | Cell conc. too low | Adjust before encapsulation |
| Low yield | Inefficient lysis | Adjust lysis conditions. **Do not use mechanical lysis with SPCs** |
| Unexpected microorganisms in data | Contamination | Laminar flow hood, decontaminate surfaces, filter tips; UV the SPC chip and Core Reagent before encapsulation. **Do not expose Shell Reagent to light** |
| High barcode cross-contamination | Cell conc. too high | Dilute to occupancy ~0.1 |
| Few SPCs recovered | SPC stickiness | Non-stick (LoBind) tubes, pipette carefully |

## 6. Notes and gaps

- Enzyme identities are proprietary throughout (WGA polymerase, primer mix, debranching enzyme,
  end-prep mix, ligase). Thermal profiles are the only mechanistic clue.
- **The 4 nt linker sequences are not disclosed** in either Atrandi document.
- The barcode plate supplies barcodes as pre-formed ligation-ready adapters; their duplex
  structure (overhang, phosphorylation, blocking) is **not shown** in this guide.
- Relevant for our protocol: we replace this kit's WGA step with **PTA**, and run **NEBNext
  library prep separately** — see the project notes.
