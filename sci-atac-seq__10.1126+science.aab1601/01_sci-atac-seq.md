# sci-ATAC-seq family: chromatin accessibility by combinatorial indexing

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or not available to us. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note. They are not yet asserted in a self-test, because
> this protocol has no `tools/` module yet (status `notes`).

This note covers three variants. They share one idea and differ in how the barcodes
get in:

| Variant | Paper | Round 1 barcode | Round 2 | Round 3 |
|---|---|---|---|---|
| **sci-ATAC-seq** (2-level) | Cusanovich 2015, *Science* | barcoded Tn5 (both ends) | indexed PCR (both ends) | – |
| **sci-ATAC-seq, 2018 protocol** | Cusanovich 2018, *Nature* | same Tn5 oligos | PCR, 10-nt indices, NextSeq | – |
| **sci-ATAC-seq3** (3-level) | Domcke 2020, *Science* | plain Nextera Tn5, then **N5 ligation** | **N7 ligation** | indexed PCR |

**sci-ATAC-seq**: Cusanovich DA, Daza R, Adey A, Pliner HA, Christiansen L, Gunderson KL,
Steemers FJ, Trapnell C, Shendure J. "Multiplex single-cell profiling of chromatin
accessibility by combinatorial cellular indexing." *Science* 2015;348(6237):910–914.
doi:[10.1126/science.aab1601](https://doi.org/10.1126/science.aab1601) · PMC4836442 ·
GEO GSE67446.

Other papers in the family:

- **Chemistry origin (CPT-seq):** Amini S, …, Shendure J, Gunderson KL, Steemers FJ.
  "Haplotype-resolved whole-genome sequencing by contiguity-preserving transposition and
  combinatorial indexing." *Nat Genet* 2014;46(12):1343–1349.
  doi:[10.1038/ng.3119](https://doi.org/10.1038/ng.3119). Here the indexed transposons,
  the two-level scheme and the custom sequencing primers were first published.
- Cusanovich DA, Reddington JP, Garfield DA, …, Shendure J, Furlong EEM. "The
  cis-regulatory dynamics of embryonic development at single-cell resolution." *Nature*
  2018;555(7697):538–542. doi:[10.1038/nature25981](https://doi.org/10.1038/nature25981).
  This is the improved two-level protocol, and its Table S12 holds the full oligo list.
- **sci-ATAC-seq3:** Domcke S, Hill AJ, Daza RM, …, Cusanovich DA, Shendure J. "A human
  cell atlas of fetal chromatin accessibility." *Science* 2020;370(6518):eaba7612.
  doi:[10.1126/science.aba7612](https://doi.org/10.1126/science.aba7612). Protocol:
  protocols.io doi:10.17504/protocols.io.be8mjhu6.
- Cited by upstream for the "semi-suppressive PCR" argument: nanoCAGE/CAGEscan,
  *Nat Methods* 2010, doi:10.1038/nmeth.1470 (not fetched for this note).

## Sources read

Fetched by `tools/get_sources.py` into
`$CHEM_DATA/sources/sci-atac-seq-family__10.1126+science.aab1601/` (never committed):

| File | What | Used for |
|---|---|---|
| `upstream_sci-ATAC-seq_family.html.txt` | scg_lib_structs page (both variants) | second source, checked below |
| `science.aab1601_PMC4836442.html.txt` | 2015 *Science* main text (author manuscript) | design, nuclei numbers, collision rate |
| `ng.3119_41588_2014_BFng3119_MOESM2_ESM.xlsx.txt` | Amini 2014 Supplementary Table 4 (oligos) | **transposon oligos, sequencing primers** |
| `ng.3119_PMC4409979.html.txt` | Amini 2014 Online Methods | transposome assembly, PCR, sequencing |
| `ng.3119_41588_2014_BFng3119_MOESM1_ESM.pdf.txt` | Amini 2014 Supplementary Figs. | Supp. Fig. 4: order of index cycles |
| `x_moesm2/nature25981-s2/Table_S12.xlsx.txt` | 2018 *Nature* Table S12, unzipped from `…MOESM2_ESM.zip` | **PCR primers (96 + 96), Tn5 oligos, seq primers** |
| `nature25981_PMC5866720.html.txt` | 2018 *Nature* Online Methods | 2018 step-by-step, cycling, NextSeq |
| `science.aba7612_PMC7785298.html.txt` | 2020 *Science* (sci-ATAC-seq3) methods | sci-ATAC-seq3 steps and read recipe |
| `scg_aba7612_domcke_table-s7.xlsx.txt` | original Domcke 2020 Supplementary Table S7, mirrored by scg_lib_structs | every sci-ATAC-seq3 splint, N5/N7 ligation oligo, PCR oligo and custom sequencing primer |
| `nature25981_…MOESM1_ESM.pdf`, `…MOESM3_ESM.xlsx`, `x_moesm2/…/Table_S11.xlsx` | reporting summary, enrichment tables, enhancer-cloning primers | not chemistry |

Could not be fetched (PMC download gate). The URLs are in `MANIFEST.tsv`:

| Missing | Why it matters |
|---|---|
| 2015 *Science* Supplementary Material (`NIHMS776051-supplement-Supplementary_Material.pdf`) and `Table_S2.txt` | **the 2015 methods and its own oligo/primer list**. 🔴 The 2015 PCR primers (8 or 10 nt indices?) are therefore unknown to us |
| sci-ATAC-seq3 supplementary PDF (`NIHMS1652075-supplement-manuscript_supplementary.pdf`) | figures and prose remain gated; the separate original Table S7 workbook is available and read |
| Amini 2014 `Figures___Tables.doc` | duplicate of the Springer supplement, probably |
| 2018 `Supplementary_tables.zip`, `Supplementary_table_4.zip` (PMC copies) | the Springer copy `MOESM2_ESM.zip` was obtained instead |
| protocols.io be8mjhu6 | detailed sci-ATAC-seq3 volumes, oligo concentrations |

---

## 1. What it is, and what is new

There is no physical single-cell isolation. Nuclei are barcoded **in bulk, per well**.
They are then pooled and redistributed at limiting numbers, and given a second barcode
per well. Any nucleus is unlikely to share its **combination** of wells with another,
so the combination acts as the cell barcode. 🟢 (2015 main text, Fig. 1A). The collision
rate is tuned by the number of nuclei per second-round well: 15–25 nuclei per well and
~11 % estimated collisions in 2015, 25 nuclei and ~12 % in 2018. 🟢

| | New thing here | Builds on |
|---|---|---|
| 1 | **Combinatorial indexing of intact nuclei** (2015). Tn5 inserts barcoded adapters inside nuclei, and the nuclei stay intact for a FACS re-sort | CPT-seq (Amini 2014) did the same on naked high-molecular-weight DNA; [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) |
| 2 | **Barcoded transposons**: 8 "T5" and 12 "T7" oligos with an 8-nt index **between** an extended primer site and an extended mosaic end. Every fragment carries the Tn5 well barcode at both ends | Nextera's s5/s7 + ME layout, lengthened on both sides of the index |
| 3 | **Custom sequencing primers** (Read 1/2 = "extended ME", Index 1/2 = their reverse complements) and a custom recipe with **dark cycles** across the constant region between the Tn5 index and the PCR index | Amini 2014 |
| 4 | (sci-ATAC-seq3) **No custom Tn5**: plain Nextera Tn5, then the 5'-OH of the transferred strand is phosphorylated with T4 PNK, and barcodes are **splint-ligated** to the s5 end (N5) and then to the s7 end (N7). Three levels, "384^3" 🟢 | the split-pool ligation idea of sci-RNA-seq3 / SPLiT-seq; [small-RNA ligation](../ref/concepts/small-rna-ligation.md) for splint ligation in general |

## 2. Oligos: two-level sci-ATAC-seq (2014 / 2015 / 2018)

### 2.1 Transposon (Tn5) oligos

🟢 Verbatim from Amini 2014 Supplementary Table 4 (sheet "Oligonucleotide_Seuqneces"
[sic]). The 2018 Table S12 (sheet "Trasposon Oligos" [sic]) says it is identical, and the
20 + 1 rows are identical character for character 🟡 *(computed)*.
Transferred strands ("Index Transposon Oligo"):

```
P5_i5_1_Universal_Connector_A_C15_ME   TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
P5_i5_8_Universal_Connector_A_C15_ME   TCGTCGGCAGCGTCTCCACGCGTACTGACGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
P7_i7_1_Universal_Connector_B_D15_ME   GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
P7_i7_12_Universal_Connector_B_D15_ME  GTCTCGTGGGCTCGGCTGTCCCTGTCCCTATCGCTCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
```

Non-transferred strand ("Universal Transposon Oligo"), written exactly as the table
gives it, with no leading slash:

```
pMENTS   5Phos/CTGTCTCTTATACACATCT
```

All 20 follow two templates 🟡 *(computed: each row was split on the constant parts,
and every one matches)*:

```
T5 (8 oligos, 63 nt):  TCGTCGGCAGCGTCTCCACGC  <t5:8>  GCGATCGAGGACGGC  AGATGTGTATAAGAGACAG
T7 (12 oligos, 69 nt): GTCTCGTGGGCTCGGCTGTCCCTGTCC  <t7:8>  CACCGTCTCCGCCTC  AGATGTGTATAAGAGACAG
```

The barcodes, in table order (P5_i5_1..8 and P7_i7_1..12), as they stand in the oligo.
🟢 They are substrings of the table rows. The list was extracted by computation.

```
t5: TATAGCCT ATAGAGGC CCTATCCT GGCTCTGA AGGCGAAG TAATCTTA CAGGACGT GTACTGAC
t7: CGAGTAAT TCTCCGGA AATGAGCG GGAATCTC TTCTGAAT ACGAATTC AGCTTCAG GCGCATTA
    CATAGCCG TTCGCGGA GCGCGAGA CTATCGCT
```

The minimum pairwise Hamming distance is 3 within each set. 🟡 *(computed)*

### 2.2 Sequencing primers

🟢 Amini 2014 Table 4 gives Read 1, Read 2 and Index 1. The 2018 Table S12 (sheet
"Sequencing Primers", "based on" Amini) adds **Index 2**:

```
Read 1 Sequencing Primer (Extended Mosaic End A)   GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
Read 2 Sequencing Primer (Extended Mosaic End B)   CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
Index 1 Sequencing Primer                          CTGTCTCTTATACACATCTGAGGCGGAGACGGTG
Index 2 Sequencing Primer   (2018 only)            CTGTCTCTTATACACATCTGCCGTCCTCGATCGC
```

### 2.3 PCR (second-round) primers, 2018

🟢 Table S12, sheets "P5 Primers" / "P7 Primers". The layout lines are copied as the
table states them:

```
P5: 5'-AATGATACGGCGACCACCGAGATCTACAC [Barcode] TCGTCGGCAGCGTC-3'
P7: 5'-CAAGCAGAAGACGGCATACGAGAT [Barcode] GTCTCGTGGGCTCGG-3'

P5_1_PCR_Primer   AATGATACGGCGACCACCGAGATCTACACCTCCATCGAGTCGTCGGCAGCGTC    Barcode CTCCATCGAG
P7_1_PCR_Primer   CAAGCAGAAGACGGCATACGAGATCCGAATCCGAGTCTCGTGGGCTCGG        Barcode CCGAATCCGA
```

There are 96 P5 primers (`P5_1..96`) and 96 P7 primers (`P7_1..96`). The full list is in
the source file. It is not copied here.

🔴 The 2015 paper's own PCR primers are in its gated supplement. Amini 2014 used the
**Illumina Nextera index kit** (FC-121-1012, 8 i5 × 12 i7) 🟢. Whether 2015 used that kit
or the 10-nt set is not known to us.

### 2.4 How they interlock: 🟡 (computed against `lib/`)

- **T5 = `nextera.S5` + `TCCACGC` + t5 + `GCGATCGAGGACGGC` + `nextera.ME`.** The first
  21 nt are what upstream calls "s5", and Amini's oligo name calls them "Universal
  Connector A". Its first 14 nt are exactly Nextera's s5.
- **T7 = `nextera.S7` + `CTGTCCCTGTCC` + t7 + `CACCGTCTCCGCCTC` + `nextera.ME`.** The
  first 27 nt are "Connector B", and the first 15 of them are exactly Nextera's s7.
- **Read 1 primer = `GCGATCGAGGACGGC` + ME** (34 nt), the 15 nt after the T5 index plus
  the ME ("C15"). **Read 2 = `CACCGTCTCCGCCTC` + ME** (34 nt, the T7 counterpart).
  Neither is `nextera.READ1_PRIMER`/`READ2_PRIMER`, so a standard Nextera Read 1/2
  primer would land **upstream of the Tn5 index** and read into it. The extended primers
  start reading at the genomic insert.
- **Index 1 = revcomp(Read 2 primer)** and **Index 2 = revcomp(Read 1 primer)**. Each
  index primer sits on the ME and reads outward across the Tn5 index into the PCR index.
- **The 2018 PCR primers are standard Nextera-shaped primers:** `illumina.P5` + 10-nt
  index + `nextera.S5` (53 nt), and `illumina.P7` + 10-nt index + `nextera.S7` (49 nt).
  In every row the barcode column is the bases between the constant parts, written in the
  primer's own orientation. All 96 barcodes of each set are distinct (minimum Hamming
  distance 3). So the PCR primers anneal with only **14 / 15 nt**, the Nextera s5/s7
  stems at the 5' end of each transposon. The Tn5 oligos are compatible with any Nextera
  index primer. This is what lets Amini use the Illumina kit.
- **The constant gap the index reads must cross** is the rest of the connector: on the i7
  side 27 nt (`S7` + `CTGTCCCTGTCC`), on the i5 side 21 nt (`S5` + `TCCACGC`).

## 3. Step by step

### 3.1 Transposome assembly (Amini 2014 Online Methods 🟢)

1. Anneal each of the 20 indexed oligos to pMENTS 1:1 at 100 µM in 10 mM Tris-HCl, 1 mM
   EDTA, 25 mM NaCl, pH 8.0. Heat to 95 °C for 5 min, then ramp at 0.1 °C/s to 25 °C.
2. Mix each annealed transposon 1:1 with EZ-Tn5 transposase (Epicentre), 12.5 µM final,
   at 37 °C for 1 h. QC on an 8 % TBE gel.
3. Make the **96 transposome mixes**: the 8 i5 and 12 i7 complexes are aliquoted 1:1 into
   a 96-well plate, giving one T5 × T7 pair per well. The working stock is 2.5 µM. (The
   methods literally say i5 into "columns 1–12" and i7 into "rows A-H, respectively",
   which cannot hold for 8 vs 12 oligos; 🟡 the sensible reading is T5 on rows A–H and T7
   on columns 1–12.) 🟡 Each well therefore holds a mix that can put **t5 on one end and
   t7 on the other**.

In 2018 the complexes were "custom and uniquely indexed Tn5 Transposomes (Illumina,
2.5 µM)". 🟢

### 3.2 sci-ATAC-seq, 2018 protocol (2018 Online Methods 🟢)

1. **Nuclei**: thaw, pellet at 500 g for 5 min, and resuspend in cold lysis buffer (10 mM
   Tris-HCl pH 7.4, 10 mM NaCl, 3 mM MgCl₂, 0.1 % IGEPAL CA-630, protease inhibitors).
   Stain with DAPI at 3 µM.
2. **Sort 2,500 nuclei per well** into 96 wells, each holding 9 µL lysis buffer and
   10 µL TD buffer (Illumina).
3. **Tagment**: 1 µL of indexed transposome per well, 55 °C for 30 min. The genomic
   fragments now carry t5/t7 at their ends, and the well identity is written into both
   ends of every fragment. 🟡
4. **Stop**: 20 µL of 40 mM EDTA + 1 mM spermidine, 37 °C for 15 min.
5. **Pool**, re-stain with DAPI, and **sort 25 nuclei per well** into new 96-well plates
   (4 plates per time point). Each well holds 12 µL of reverse-crosslinking buffer:
   11 µL EB, 0.5 µL Proteinase K (20 mg/mL) and 0.5 µL 1 % SDS.
6. **Reverse crosslink** overnight at 65 °C. The embryos were formaldehyde-fixed.
7. **Indexed PCR** in the same well: primers at 0.5 µM each (one P5 + one P7 from Table
   S12), 7.5 µL NPM and BSA (2× final). Cycle 72 °C 3 min (gap fill-in), 98 °C 30 s, then
   **15–25 cycles** of 98 °C 10 s / 63 °C 30 s / 72 °C 1 min. The cycle number comes from a
   SYBR test on spare wells.
8. **Pool** all wells, clean up on 4 × Zymo DCC-5 columns and then AMPure, and run
   Bioanalyzer QC.

The 2015 main text gives only the outline 🟢: 96 barcoded-Tn5 wells; pool; dilute; 15–25
nuclei per second-round well by FACS; lyse; PCR with indexed primers; pool; sequence.
The 2015 tagmentation and PCR conditions are in the gated supplement 🔴.

Plate bookkeeping 🟡 *(from the "Note" column of Table S12)*: every P5 primer is used on
**one** plate (8 per plate, 12 plates). The P7 sets of 12 are reused between some plates,
for example "2-4hr_Plate1, 10-12_Plate1". The unique P5 keeps each plate × well
combination unique.

### 3.3 The fragments that amplify: 🟡 (inferred, the standard Tn5 argument)

Each Tn5 well has both T5 and T7 complexes, so fragments come out as T5–T5, T7–T7 or
T5–T7, with the 9-nt duplicated target gap at each end
([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)). Only T5–T7 fragments are
amplified exponentially by a P5 + P7 primer pair. Upstream says the same-end products
are "not amplifiable due to semi-suppressive PCR". That claim is upstream's. The
fetched papers do not state it.

## 4. Final library (two-level), as segment lists: 🟡 (assembled from the oligos)

Top strand, 5'→3'. Index lengths are from 2018 (10-nt PCR indices). With the Nextera kit
of Amini 2014, `i5`/`i7` are 8 nt.

```
5'- P5 (AATGATACGGCGACCACCGAGATCTACAC)
  · i5 <10>                       PCR well, P5 side
  · Connector A (TCGTCGGCAGCGTCTCCACGC)     = s5 + TCCACGC
  · t5 <8>                        Tn5 well, P5 side
  · C15-A (GCGATCGAGGACGGC) · ME (AGATGTGTATAAGAGACAG)   = Read 1 primer site
  · <genomic insert>
  · ME' (CTGTCTCTTATACACATCT) · C15-B' (GAGGCGGAGACGGTG) = Index 1 primer site
  · t7' <8>                       Tn5 well, P7 side
  · Connector B' (GGACAGGGACAGCCGAGCCCACGAGAC)  = revcomp(s7 + CTGTCCCTGTCC)
  · i7' <10>
  · P7' (ATCTCGTATGCCGTCTTCTGCTTG) -3'
```

Insert-flanking length (P5 → ME, without i5/t5) = 29 + 21 + 15 + 19 = 84 nt. On the other
side it is 19 + 15 + 27 + 24 = 85 nt. With 10 + 8 nt of index on each side, the
non-genomic total is **205 nt**. 🟡 *(computed)*

The **upstream page agrees** with this final library segment for segment, including
10-nt i5/i7, 8-nt Tn5 barcodes and the 21/27-nt connectors.

## 5. Sequencing (two-level)

- **Amini 2014 🟢**: HiSeq 2000, custom recipe, custom Read 1 (on cBot), Index 1 and
  Read 2 primers, 51-bp reads. Supp. Fig. 4 gives the cycle order: read 1 (1–51); index 1
  = transposon i7 (52–59) then PCR i7 (60–67); index 2 = **PCR i5 (68–75) then
  transposon i5 (76–83)**; read 2 (84–134). The dark cycles across the connector are not
  numbered.
- **2018 🟢**: NextSeq High Output 300-cycle kit, 1.5 pM, custom primers and recipe,
  50 bp each end "in addition to the barcodes". The cell barcode is 4 parts: Tn5 + PCR
  on the P5 side and Tn5 + PCR on the P7 side. Each part is matched within 3 edits.
- **Index read lengths 🟡 *(computed, matches upstream)***: Index 1 = 8 (t7) + **27
  dark** + 10 (i7) = 45 cycles. Index 2 = 8 (t5) + **21 dark** + 10 (i5) = 39 cycles.
  Upstream gives exactly 45 and 39. The papers do not give these numbers.
- **Why Amini reads i5 in the opposite order 🟡 (inferred)**: Amini lists **no Index 2
  primer**. On the HiSeq 2000 the i5 read is primed from the P5 end (the standard
  forward-strand workflow), so it meets the PCR i5 first and then the Tn5 index. The 2018
  NextSeq run uses a custom Index 2 primer = revcomp(Read 1 primer). It anneals on the ME
  side of the top strand and reads Tn5 index first, then i5, which is the order upstream
  draws.

## 6. sci-ATAC-seq3 (Domcke 2020)

### 6.1 Oligos: 🟢 Supplementary Table S7

As upstream writes them (spaces added around the placeholders):

```
3LV2_N5_Splint                         GCCGACGACTGATTA/3ddC/
3LV2_N7_Splint                         CACGAGACGACAAGT/3ddC/
3LV2_N5_Oligo_{A01..H12}_Plate{1..4}   CACCGCACGAGAGGT [10-bp N5 barcode] GTAATCAG     (384)
3LV2_N7_Oligo_{A01..H12}_Plate{1..4}   CAGCACGGCGAGACT [10-bp N7 barcode] GACTTGTC     (384)
3LV2_P5_PCR_{A01..H12}                 AATGATACGGCGACCACCGAGATCTACAC [10-bp i5] CACCGCACGAGAGGT   (96)
3LV2_P7_PCR_{A01..H12}                 CAAGCAGAAGACGGCATACGAGAT [10-bp i7] CAGCACGGCGAGACT        (96)
Read 1 sequencing primer               TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Read 2 sequencing primer               GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
Index 1 sequencing primer              CTCCGAGCCCACGAGACGACAAGTC
Index 2 sequencing primer              ACACATCTGACGCTGCCGACGACTGATTAC
```

### 6.2 How they interlock: 🟡 (computed against `lib/`)

- **N5 splint** (counting the 3' ddC as the 16th base, `…TTAC`): its reverse complement
  is `GTAATCAG` + `TCGTCGGC`. That is the **last 8 nt of the N5 oligo** followed by the
  **first 8 nt of `nextera.S5`**. The splint bridges the nick between the N5 oligo's
  3'-OH and the PNK-phosphorylated 5' end of the s5 transferred strand, with 8 bp on each
  side.
- **N7 splint** likewise: revcomp = `GACTTGTC` + `GTCTCGTG`, the last 8 nt of the N7
  oligo and `nextera.S7[:8]`. Each splint pairs with only one of the two Nextera ends.
  This is why round 1 tags s5 ends only and round 2 s7 ends only. The 3'-ddC stops the
  splint from being extended or ligated itself.
- **Read 1 / Read 2 primers = `nextera.READ1_PRIMER` / `READ2_PRIMER`** exactly: plain
  Nextera, because the Tn5 is plain Nextera.
- **Index 1 = `nextera.ME_RC[-2:]` (`CT`) + `nextera.S7_RC` + revcomp(N7 tail)** (25 nt).
  **Index 2 = `nextera.ME_RC[-8:]` + `nextera.S5_RC` + revcomp(N5 tail)** (30 nt). Each
  ends right at the ligated barcode.
- **PCR primers anneal on the 15-nt 5' heads** of the N5/N7 oligos. So 15 dark cycles
  separate the ligation barcode from the PCR index on each side. That matches the paper's
  recipe 🟢 "index 1: 10 cycles+15 dark cycles+10 cycles, index 2: 10 cycles+15 dark
  cycles+10 cycles".
- N5 and N7 oligos are 33 nt each.

### 6.3 Step by step (2020 methods 🟢, condensed)

1. Fixed, frozen nuclei: thaw, Omni lysis buffer, dilute in ATAC-RSB + 0.1 % Tween-20.
2. **Tagment in bulk per sample**: 200,000 nuclei in TD buffer, 1× DPBS, 0.01 %
   digitonin and 0.1 % Tween-20, split over 4 wells (50,000 nuclei/well). Add 2.5 µL
   Nextera v2 enzyme per well and incubate 55 °C for 30 min. Stop with 40 mM EDTA + 1 mM
   spermidine, 37 °C for 15 min. There is no barcode at this step.
3. **Pool, wash, phosphorylate**: PNK buffer, 1 mM rATP, T4 PNK; each sample is spread
   over a total of 16 wells across the four 96-well plates, 37 °C for 30 min.
4. **N5 ligation**: T7 ligase buffer + N5 splint + T7 DNA ligase, then one of 384 N5
   oligos per well, 25 °C for 1 h. Stop with EDTA/spermidine.
5. **Pool, wash, N7 ligation**: N7 splint + T7 ligase, redistribute over 4 × 96 wells,
   one of 384 N7 oligos per well, 25 °C for 1 h. Stop.
6. **Pool, redistribute 1,000–3,000 nuclei per well** (4 × 96 wells). Reverse
   crosslink with EB + Proteinase K + 1 % SDS at 65 °C for 16 h.
7. **Indexed PCR**: NPM, BSA, one indexed P5 and one indexed P7 oligo. The cycle number
   comes from a SYBR test PCR. Pool, clean up with Zymo DCC-5 and then 1× AMPure, and run
   Tapestation QC.
8. **Sequence**: NovaSeq 6000, custom primers and recipe. Read 1 51, read 2 51, index 1
   10+15 dark+10, index 2 10+15 dark+10.

What happens at the ends 🟡 (inferred, agrees with upstream): standard Nextera leaves
the 5' end of each transferred strand unphosphorylated, so PNK is needed before ligation.
The barcode oligos are ligated to the top-strand 5' ends only. The 9-nt gap and the
missing complement are filled by the 72 °C step before PCR. Only fragments with an s5
end (N5-ligated) and an s7 end (N7-ligated) carry both PCR handles.

### 6.4 Final library (sci-ATAC-seq3): 🟡 assembled from 🟢 Table S7 oligos

```
5'- P5 · i5 <10> · N5 head (CACCGCACGAGAGGT) · N5 <10> · GTAATCAG
  · s5 (TCGTCGGCAGCGTC) · ME (AGATGTGTATAAGAGACAG)          = Read 1 site
  · <genomic insert>
  · ME' (CTGTCTCTTATACACATCT) · s7' (CCGAGCCCACGAGAC) · GACAAGTC   (Index 1 primer = CT + s7' + GACAAGTC)
  · N7' <10> · N7 head' (AGTCTCGCCGTGCTG) · i7' <10> · P7' -3'
```

The cell barcode is i5 + N5 + N7 + i7. The paper calls the design "384^3" 🟢: 384 N5
wells, 384 N7 wells and 4 × 96 PCR wells, each PCR well with its own P5/P7 pair 🟡. The sample has no barcode of its own. Each of the 24 samples per batch is
spread over 16 phosphorylation/N5 wells (24 × 16 = 384), so the N5 barcode identifies the
sample. 🟡 *(computed from the 2020 methods)*
How samples map to plates is in Table S1 / S7 🔴.

## 7. Checking the upstream page against the papers

| Upstream claim | Paper / table | Verdict |
|---|---|---|
| Barcoded Tn5 s5/s7 templates, ME bottom `/Phos/CTGTCTCTTATACACATCT` | Amini Table 4 = 2018 Table S12 (`5Phos/…`) | **agrees** |
| Read 1, Read 2, Index 1, Index 2 primers (two-level) | 2018 Table S12 | **agrees** (Index 2 is not in Amini 2014) |
| P5 index primer ends `…TCGTCGGCAGCGTCTCCACGC` (21 nt); P7 ends `…GTCTCGTGGGCTCGGCTGTCCCTGTCC` (27 nt) | 2018 Table S12: P5 ends `TCGTCGGCAGCGTC` (14), P7 ends `GTCTCGTGGGCTCGG` (15); Amini used Nextera-kit primers (also 14/15) | **disagrees**: upstream extends the PCR primers over the whole connector. The final library is the same either way, because the rest of the connector comes from the transposon |
| i5/i7 are 10 nt | 2018 Table S12: 10 nt; Amini 2014: Nextera kit (8 nt) | agrees for 2018. 2015 unknown 🔴 |
| Index 1 = 45 cycles (8 + 27 dark + 10), Index 2 = 39 (8 + 21 dark + 10) | not stated in the papers | consistent with the computed gaps 🟡 |
| Same-end products not amplifiable ("semi-suppressive PCR") | not stated in the fetched papers | upstream interpretation 🟡 |
| sci-ATAC-seq3: PNK, N5 then N7 splint ligation, 384 × 4 plates, reverse crosslink, indexed PCR, 10+15+10 index reads, 50-cycle reads | 2020 methods: same steps. Reads are **51** cycles each; upstream says 50 | agrees, except the 50/51 read length |
| sci-ATAC-seq3 oligo sequences | original Table S7 workbook mirror | **agree** 🟢 |

## 8. Open questions

- 🔴 The 2015 *Science* supplement (methods and Table S2): its tagmentation conditions,
  PCR primers (Nextera kit, 8 nt, or the 10-nt set?), cycle numbers and sequencer.
- 🟡 Barcode orientation in the reads: whether i5 is read as written in the primer or
  as its reverse complement depends on the instrument's i5 workflow. This was not
  checked against demultiplexing code (github shendurelab/human-atac;
  atlas.gs.washington.edu).
- 🔴 Why Amini's connectors were lengthened (s5 + 7 nt, s7 + 12 nt) and why the
  "C15" extensions were added before the ME. Presumably the reason is primer Tm for the
  custom read and index primers. Not stated.
- 🟡 Upstream draws the Index 2 primer of the two-level method annealed as a 19-nt ME
  only. The real primer (34 nt) also covers C15-A.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the upstream page, PMC full text of all four papers and
the Springer supplements. All PMC-hosted supplements (2015 methods, sci-ATAC-seq3
Table S7) hit the download gate. The Europe PMC `…_supplementary.zip` for the 2018 paper
was cut off at exactly 50,000,000 bytes and is not a valid zip. The oligo table (Table
S12) sat inside the Springer `MOESM2_ESM.zip`, which neither `get_sources.py` nor
`doctext.py` unpacks. It was unzipped by hand into `x_moesm2/` and converted with
`doctext.py`. After that, `scrape_primers.py` found it and ranked it with Amini's Table 4
and the upstream page. The 2018 `.ppt` extended-data files were not converted. The
reaction order came from reading the methods.
