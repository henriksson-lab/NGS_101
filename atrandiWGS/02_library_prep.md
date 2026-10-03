# Atrandi Single-Microbe DNA Barcoding — Library Prep for Sequencing — quick reference

Source: `DGPM02323206001_Single_Microbe_DNA_Barcoding_Kit_Library_Prep_For_Seq_V3.pdf`
Doc. No. **DGPM02323206001**, Revision **V3**, 03 November 2023.
(V1 2 Oct 2023; V2 13 Oct 2023 — *library amplification enzyme changed to NEBNext Ultra II Q5
Master Mix*; V3 — *updated demultiplexing by multiple barcode combinations*.)
Atrandi Biosciences / Droplet Genomics. Research Use Only.

> Input: barcoded SPCs from the Single-Microbe DNA Barcoding Kit (CKP-BARK1) —
> see [`01_barcoding_kit.md`](01_barcoding_kit.md).

> **Evidence marking.** 🟢 = verbatim from the source document · 🟡 = derived or inferred
> (reasoning given; arithmetic verified where possible) · 🔴 = not published anywhere.
> Any claim not marked 🟢 must reach the final page carrying its uncertainty.


---

## 1. Library structure (Figure 2)

```
5'- P7 - Index i7 - Adapter - D - C - B - A - {Fragment} - Adapter - Index i5* - P5 -3'
                             \_______________/
                                cell barcode                         (* i5 optional)
```

So on the i7 side, reading inward: adapter, then the four barcode blocks in the order
**D, C, B, A** (reverse of the ligation order A→B→C→D), then the genomic fragment.

### Ligation Adapter sequence

| Name | Sequence (5'→3') | Purification |
|---|---|---|
| Top adapter | `/5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG*T` | HPLC |
| Bottom adapter | `/5AmMC6/GCTCTTCCGATCT` | HPLC |

`*` = phosphorothioate linkage. `/5Phos/` = 5' phosphate. `/5AmMC6/` = 5' amino-C6 (blocks ligation
at that end). The **bottom** adapter's 3'-terminal T is the single-base overhang that pairs with
the dA-tailed insert.

🟡 **It is not a forked / Y adapter** — derived, but verified base-by-base in code (Atrandi give
only the two oligo sequences, not the structure). The bottom oligo is 13 nt and *all 12*
of its 5'-proximal bases pair with the top strand's first 12, leaving only **one** single-stranded
arm — a 12 bp stem plus a 20 nt 3' tail on the top strand:

```
                 <------ 12-bp stem ------>  <--------- 20-nt single-stranded 3' arm --------->
 TOP  5'-/5Phos/ G A T C G G A A G A G C      G T C G T G T A G G G A A A G A G T G T -3'
                 | | | | | | | | | | | |
 BOT        3'-T C T A G C C T T C T C G /5AmMC6/-5'
              ^
              3'-T overhang = the ligation end
```

Verified: `revcomp(GCTCTTCCGATC) = GATCGGAAGAGC` (12/12 Watson-Crick);
`revcomp(GTCGTGTAGGGAAAGAGTGT) = ACACTCTTTCCCTACACGAC` = the **first 20 nt of the TruSeq Read 1
primer**; and `GCTCTTCCGATCT` = the **last 13 nt of the TruSeq Read 1 primer**.

**Consequence — this adapter installs only the TruSeq Read 1 / P5 side, at both ends of every
fragment. It carries no Read 2 / P7 arm at all.** After ligation the top strand read outward from
the insert is `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` (leading A = the dA tail) — character-for-character
Illumina's canonical Read 2 adapter-trimming sequence.

So the **P7 / Read 2 side must be installed by the split-pool barcode cassette**, not here. See
`04_nebnext_illumina.md` §4.

🟡 **INFERRED (Atrandi give no rationale).** The `/5AmMC6/` on the bottom strand sits at the fork
point and makes that 5' end permanently ligation-incompetent: unlike a bare 5'-OH it cannot be rescued by the polynucleotide-kinase
activity carried over from the FS end-repair mix, so adapter-adapter dimers cannot form.

### Indexing primers (index underlined in the PDF)

| Name | Sequence (5'→3') | Purification |
|---|---|---|
| Index i1 P7 | `CAAGCAGAAGACGGCATACGAGAT` **`AACCTG`** `GTGACTGGAGTTCAGACGTGTGCTCT*T` | HPLC |
| P5 | `AATGATACGGCGACCACCGAGATCTACAC` `TCTTTCCCTACACGA*C` | HPLC |

- `CAAGCAGAAGACGGCATACGAGAT` = Illumina P7; `AATGATACGGCGACCACCGAGATCTACAC` = Illumina P5.
- The i7 index here is **6 nt** (`AACCTG`); the P5 primer carries **no index** as written.
- Note: these are given as *examples*. The protocol actually specifies NEBNext Multiplex Oligos
  Dual Index Primers Set 2 (NEB #E7780S) for the PCR step.

> Atrandi's note: "Atrandi Biosciences' products may work with a number of library preparation
> technologies. Library preparation using the kit sold by NEB for use with Illumina sequencing is
> presented here as a validated example."

## 2. Reagents

| Step | Reagent | Vendor / P/N |
|---|---|---|
| Barcoded DNA purification | 1X STOP Buffer, 1X Wash Buffer, NFW | Atrandi (in CKP-BARK1) |
| | Release Reagent | Atrandi (in SPC Innovator Kit CKN-G11 / SPC Generation Kit CKP-G34) |
| | AMPure XP Reagent | Beckman Coulter A63880 |
| Fragmentation | NEBNext Ultra II FS DNA Library Prep Kit for Illumina | NEB **E7805S** |
| Adapter ligation | Ligation Adapter | Atrandi (in CKP-BARK1) |
| Purification | 80% ethanol | MLS |
| PCR | NEBNext Multiplex Oligos for Illumina (Dual Index Primers Set 2) | NEB **E7780S** |
| PCR | NEBNext Ultra II Q5 Master Mix | NEB (part of E7805S) |

Also needed: magnetic rack; 0.2 mL non-stick PCR tubes; LoBind 1.5 mL tubes; Qubit + dsDNA HS
(Thermo Q33238 / Q32851); Agilent 2100 Bioanalyzer + HS DNA Kit (G2939BA / 5067-4626).

## 3. Protocol

### I. Barcoded DNA purification
1. Pellet SPCs stored in 1X SB, 1 min @ 1000 × g; discard supernatant.
2. Wash **3×** with 1 mL 1X WB.
3. Bring to 200 µL with 1X WB.
4. **Split into 4 × 0.2 mL PCR tubes, 50 µL each** — all subsequent volumes are *per tube*.
6. Add **2 µL Release Reagent**, 5 min at RT. (2 µL per ≤50 µL SPCs; scale proportionally.)
7. Add 48 µL NFW → 100 µL total.
8–9. Add **80 µL AMPure XP** (0.8X), mix.
10. 5 min RT.
11. Magnet; remove supernatant.
12–13. **2× wash** with 180 µL 80% ethanol, ~30 s each.
14. Spin, magnet, remove residual, air-dry 1 min.
15. Elute in **27 µL NFW**, 5 min RT.
16–17. Magnet until clear; transfer eluate.
18. Quantify by Qubit dsDNA HS.

**Stop point:** −20 °C up to a week.

> 🟡 **INFERRED.** The Release Reagent presumably dissolves the hydrogel shell, freeing barcoded
> DNA into solution. 🔴 Composition undisclosed.

### II. Fragmentation (NEBNext Ultra II FS)

| Reagent | µL |
|---|---|
| Purified DNA | 26 |
| NEBNext Ultra II FS Reaction Buffer | 7 |
| NEBNext Ultra II FS Enzyme Mix | 2 |
| **Total** | **35** |

| Temp | Time |
|---|---|
| 37 °C | **10 min** |
| 65 °C | 30 min |
| 4 °C | hold |

- **CRITICAL:** pre-heat the block to 37 °C before moving samples from ice.
- **CRITICAL:** exactly 10 min at 37 °C (this sets fragment size).
- **CRITICAL:** vortex the FS Enzyme Mix 5–8 s before use, keep on ice.

> One-tube enzymatic fragmentation + end repair + dA-tailing.

### III. Ligation

| Reagent | µL |
|---|---|
| Fragmented DNA (step II) | 35 |
| NEBNext Ultra II Ligation Master Mix | 30 |
| NEBNext Ligation Enhancer | 1 |
| **Ligation Adapter (1.5 µM)** | 2.5 |
| **Total** | **68.5** |

**20 °C, 15 min**, then 4 °C hold.

> The Atrandi single-tailed adapter is used here instead of the NEBNext hairpin adaptor — so there
> is **no USER enzyme step** (there is no hairpin loop to open).

### IV. 0.8X AMPure XP purification
1. Add 31.5 µL NFW to the 68.5 µL ligation → 100 µL.
2–3. Add **80 µL AMPure XP**, mix.
4. 5 min RT. 5. Magnet, discard supernatant.
6–7. 2× 180 µL 80% ethanol, ~30 s.
8. Dry 1 min.
9. Elute in **40 µL NFW**, 5 min RT. 10–11. Magnet, transfer.

### V. PCR

| Reagent | µL |
|---|---|
| Purified DNA | 40 |
| Indexing primer mix (5 µM each) | 10 |
| NEBNext Ultra II Q5 Master Mix | 50 |
| **Total** | **100** |

| Step | Temp | Time | Cycles |
|---|---|---|---|
| Initial denaturation | 98 °C | 45 s | 1 |
| Denaturation | 98 °C | 20 s | |
| Annealing | **54 °C** | 30 s | **8–14** |
| Extension | 72 °C | 20 s | |
| Final extension | 72 °C | 2 min | 1 |
| Hold | 4 °C | — | — |

\*14 cycles for <100 ng input; 8 cycles for 500 ng input.

### VI. 0.6–0.8X double-sided AMPure XP
1–4. Add **60 µL beads** to the 100 µL PCR (0.6X), 5 min RT, magnet — **keep the supernatant**.
5–6. Add **20 µL fresh beads** to that supernatant (→0.8X), 5 min RT.
7. Magnet, discard supernatant.
8–10. 2× 180 µL 80% ethanol, ~30 s.
11–12. Remove residual, dry 1 min.
13. Elute in **30 µL NFW**, 5 min RT. 14. Magnet.
15. Qubit dsDNA HS.
16. Bioanalyzer HS DNA for QC (Figure 4 shows a typical trace).

## 4. Sequencing

| Platform | Loading conc. | PhiX |
|---|---|---|
| MiSeq | 10 pM | 1% |
| NextSeq 550 | 1.8 pM | 1% |
| NovaSeq (SP flow cell) | 300 pM | 1% |

**Paired-end run configuration:**

| Read | Length | Content |
|---|---|---|
| Read 1 | **128 bp** | genomic insert |
| i7 index | 6 bp (optional) | sample index |
| i5 index | 6 bp (optional) | sample index |
| Read 2 | **172 bp** | **cell barcode + genomic insert** |

Index lengths depend on the indexing primers used.

## 5. Barcode structure in Read 2 (Figure 6)

Four 8 nt barcode parts with 4 nt linkers between them: **8 + 4 + 8 + 4 + 8 + 4 + 8**
(a 4th linker follows barcode A, before the insert).

Read-2 offsets, **0-based, half-open**, as given in the PDF pseudocode:

| Part | Slice | 1-based positions |
|---|---|---|
| Barcode **D** | `read[0:8]` | 1–8 |
| linker | — | 9–12 |
| Barcode **C** | `read[12:20]` | 13–20 |
| linker | — | 21–24 |
| Barcode **B** | `read[24:32]` | 25–32 |
| linker | — | 33–36 |
| Barcode **A** | `read[36:42]` ⚠ | 37–44 |
| linker | — | 45–48 |

> **Slip in the PDF.** It states barcode A is `read[36:42]` — only 6 nt, contradicting "each barcode
> part comprises 8 nucleotides". Read it as **`read[36:44]`**. The manual is not exact throughout;
> treat the Bascet source as authoritative for offsets.

> ⚠ The example sequence rendered in Figure 6
> (`GAACTGGAGAGTATGCGCCCATTCCACCCGTATTCGGGCCGTTGT…`) is **illustrative only** — none of its
> 8-mers at the barcode offsets match the A/B/C/D whitelists in Tables 2–3. Do **not** read the
> linker sequences off that figure. **The real linker sequences are not disclosed in either
> Atrandi document.**

## 6. Barcode whitelists (Tables 2 and 3)

All sequences 5'→3', 8 nt. Well position refers to the barcode plate.

### Barcode A (columns 1–3)
| Well | Seq | Well | Seq | Well | Seq |
|---|---|---|---|---|---|
| A1 | GTAACCGA | A2 | TCCTCAAC | A3 | TGGTCTCA |
| B1 | GACAGCAT | B2 | GATGGTCT | B3 | CATACCGT |
| C1 | CGGTAAGT | C2 | GTGACTCT | C3 | CTGTGAAC |
| D1 | CACTGACA | D2 | AGATACCG | D3 | TTGGATGC |
| E1 | AGACTCAC | E2 | CAGCAACT | E3 | CACTGTTG |
| F1 | TCAGAGGA | F2 | GCAGACTA | F3 | TGACGAAG |
| G1 | GTTGTCTG | G2 | GTATTGCC | G3 | AGTCACTG |
| H1 | TGATTCGG | H2 | ATGTCTGC | H3 | ACAATCCG |

### Barcode B (columns 4–6)
| Well | Seq | Well | Seq | Well | Seq |
|---|---|---|---|---|---|
| A4 | TACAACCG | A5 | GCTGGATA | A6 | CATCGTTG |
| B4 | TCTGGAAC | B5 | CGAACTTG | B6 | AGAATGCC |
| C4 | GCAATGAC | C5 | TCCAGTCT | C6 | CACCATCT |
| D4 | CGCCAATA | D5 | AGGTCGTA | D6 | TAACGGAG |
| E4 | TCTGCTTG | E5 | AGCATTGG | E6 | AGGACATC |
| F4 | AGACGAAC | F5 | ACAGTAGC | F6 | CAGGTGAA |
| G4 | GTGTGATG | G5 | TGCTATGG | G6 | AGACACCA |
| H4 | GGAACTGT | H5 | TATGCGAC | H6 | TGTTGGAC |

### Barcode C (columns 7–9)
| Well | Seq | Well | Seq | Well | Seq |
|---|---|---|---|---|---|
| A7 | TACAGCAG | A8 | TTCGGTAG | A9 | GATACCGA |
| B7 | GGAACCAA | B8 | CTTCAGCA | B9 | ACCTTCAG |
| C7 | GGACGAAT | C8 | AACTGCTC | C9 | CACCTTGA |
| D7 | CTGCTCAA | D8 | AGACCGAT | D9 | TGGTCACT |
| E7 | CGACATTC | E8 | AACGATGG | E9 | AGATGGTC |
| F7 | TTGGTAGC | F8 | GTCTGCTA | F9 | ATGGTGTG |
| G7 | CGCTACTA | G8 | AAGGTGAC | G9 | CCATACGA |
| H7 | CGACAAGA | H8 | AGGCTAAC | H9 | AGGCATTG |

### Barcode D (columns 10–12)
| Well | Seq | Well | Seq | Well | Seq |
|---|---|---|---|---|---|
| A10 | GTAATGCC | A11 | AGGTCCAA | A12 | ATCGACAG |
| B10 | GTGGTGAT | B11 | GGTCTCAT | B12 | GGTTGATG |
| C10 | CGACAGTA | C11 | TACCACAG | C12 | ATATGCGG |
| D10 | GTGTAAGC | D11 | ATACGACC | D12 | GGTCATTG |
| E10 | TGTGGTTG | E11 | CGGTAGAA | E12 | TATCCGGT |
| F10 | ACCTCGTA | F11 | GTGAGACT | F12 | CAACTTGC |
| G10 | CGCTAGTT | G11 | TATCAGCC | G12 | ACTAACCG |
| H10 | AATGGCAG | H11 | CAACCAAC | H12 | AAGGCTTG |

## 7. Data analysis

- R2 carries the combinatorial barcode; R1 is the mate.
- Recommended workflow: **sample a subset of reads → demultiplex to discover observed barcode
  combinations → whitelist → split the full FastQ**. Rationale: far fewer barcodes to process,
  separates demultiplexing logic from I/O, and allows custom filtering without a full pass.
- Recommended strategy: **cascading** — match D, then C, then B, then A at their fixed offsets,
  rather than matching the whole 44-mer at once. Allow error (e.g. Hamming distance 1).
- Suggested tools: Pheniqs, Seqkit, Samtools, or awk + GNU Parallel. FastQC for QC.
  Prefer BAM output for the split step.
- Suggested hardware: ≥72 vCPU, 144 GB RAM, 1 TB storage.

## 8. Troubleshooting (Table 4)

| Issue | Cause | Solution |
|---|---|---|
| Underclustering on flow cell | Library has no correct-structure barcodes | Check all barcoding + adapter ligation steps were done; more PCR cycles in step V |
| Low library yield | Undercycling | More PCR cycles in step V |
| | Too little barcoded input DNA | Measure DNA after step I; repeat barcoding with higher cell conc. |
| Library size wrong | Incorrect fragmentation | Adjust step II; see the NEBNext Ultra II FS manual |

## 9. Notes and gaps

- The 4 nt **linker sequences are not disclosed** — needed to draw the construct base-by-base.
- The **barcode oligo duplex structure** (overhangs, phosphorylation, blocking groups) is not shown;
  only the final barcode 8-mers are given.
- `read[36:42]` in the PDF pseudocode is wrong; should be `read[36:44]`.
- The Figure 6 example sequence is not a real barcoded read.
- Relevant for our protocol: we use **PTA** for genome amplification and run **NEBNext library prep
  separately** — the NEBNext hairpin adaptor + USER workflow and the dual-index primer sequences
  differ from the Atrandi example adapter above.

---

## 10. Discrepancy with the Bascet whitelist (added after cross-check)

> **Scope:** this affects *demultiplexing only*. The chemistry page draws barcodes as placeholders
> (`AAAAAAAA`/`BBBBBBBB`/`CCCCCCCC`/`DDDDDDDD`), so the real sequences are not needed there.

The barcode whitelist shipped in the authors' own demultiplexer
([Bascet](https://github.com/henriksson-lab/bascet),
`crates/bascet-cli/src/barcode/atrandi_barcodes.tsv`, archived here as
`ref/bascet_atrandi_barcodes.tsv`) does **not** fully agree with Tables 2–3 above.
Verified by direct set comparison:

| Round | Bascet `pos` | vs. PDF table | Result |
|---|---|---|---|
| B | 2 | Barcode B | **identical set** (24/24) |
| C | 3 | Barcode C | **identical set** (24/24) |
| D | 4 | Barcode D | **identical set** (24/24) |
| A | 1 | Barcode A | **only 12/24 shared** |

**Round 1 / Barcode A — the 12 that differ:**

| Only in Bascet | Only in the Atrandi PDF |
|---|---|
| GTCCGATT, TTGACCAC, TCCAGGAT | AGACTCAC, CAGCAACT, CACTGTTG |
| CGGTTGAT, CGACCTAT, TTCCACTC | TCAGAGGA, GCAGACTA, TGACGAAG |
| TGACAGTG, GGCATCAA, ACATCGTC | GTTGTCTG, GTATTGCC, AGTCACTG |
| GTCGGTAA, TTGATGGC, GTTACGGT | TGATTCGG, ATGTCTGC, ACAATCCG |

**Well mapping also differs.** Even where the sets match, the two sources disagree on which
well holds which sequence: Bascet enumerates each 3-column block **row-major** over the order
the PDF lists **column-major**. E.g. for barcode B the PDF gives A5 = `GCTGGATA` while Bascet
gives A5 = `TCTGGAAC` (= the PDF's B4); Bascet's A6 = `GCAATGAC` = the PDF's C4.

**Consequences.** For demultiplexing only the *set* matters, so rounds B/C/D are unaffected and
round A is not: 12 of its sequences would fail to match. For tracing a cell back to a physical
well, the mapping matters everywhere.

**Unresolved** — most likely a kit version/lot difference in the round-1 plate, or an error in one
of the two sources. Do not assume either is authoritative for our page until checked against the
actual plate used. Resolve by checking observed round-1 8-mers in our own R2 data against both lists.
