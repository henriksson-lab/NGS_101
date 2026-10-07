# PETRI-seq — split-pool single-cell RNA-seq for bacteria

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes`). Claims that rest only on the upstream
> scg_lib_structs page are 🟡 (secondary source).

**PETRI-seq** (Prokaryotic Expression-profiling by Tagging RNA In situ and sequencing) —
Blattman SB, Jiang W, Oikonomou P, Tavazoie S. "Prokaryotic single-cell RNA sequencing by
in situ combinatorial indexing." *Nature Microbiology* 5, 1192–1201 (2020).
doi:[10.1038/s41564-020-0729-6](https://doi.org/10.1038/s41564-020-0729-6) · PMID 32451472 ·
PMC8330242 (author manuscript NIHMS1587539). Data: GEO GSE141018. The catalogue lists no
other paper for this method. Upstream scg_lib_structs drew it from the authors' protocol
website (oligo table and protocol accessed 26 Jan 2024), not from the paper.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/petri-seq__10.1038+s41564-020-0729-6/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41564-020-0729-6_PMC8330242.html.txt` | full text (PMC author manuscript) incl. Methods and Extended Data legends | **all methods**, reaction order, conditions, sequencing |
| `s41564-020-0729-6_41564_2020_729_MOESM1_ESM.pdf(.txt)` | Supplementary Information: Supp. Figs 1–2, **Supp. Table 1** (costs), **Supp. Table 2** (single-tube oligos) | linkers, blockers, TSO, PCR/qPCR oligos |
| `s41564-020-0729-6_41564_2020_729_MOESM3_ESM.xlsx(.txt)` | **Supp. Tables 3–5**: 96-well barcode oligos, experiment overview, statistics | **every barcode oligo** (3 × 96), which experiment used which variant |
| `s41564-020-0729-6_41564_2020_729_MOESM2_ESM.pdf.txt` | Reporting summary | sequencer, software, GEO |
| `upstream_PETRI-seq.html.txt` | scg_lib_structs page | second source, checked below |
| `MOESM4_ESM.gz`, `MOESM5..16_ESM.xlsx` | count matrix (Supp. Table 6), source data | not chemistry, not read |

Not fetched / not usable:

| What | Why |
|---|---|
| `…NIHMS1587539-supplement-…_Supp_Fig1-2-Supp_Tab1-2.pdf`, `…_Supp_Tab3-5.xlsx`, `…_Supp_Tab6.gz` | saved, but each is the PMC JavaScript download-gate HTML page, not the file. Same content is in the Springer `MOESM1` / `MOESM3` / `MOESM4` files, so nothing is lost. |
| `…_SourceDataED1..10.xlsx`, `…_SourceDataFig2/3.xlsx` (PMC) | `(manual)`, download gate; source data, not chemistry |
| Publisher PDF of the main text | not needed (PMC full text read) |
| Authors' protocol website (upstream's source) | not fetched; upstream quotes it 🔴 |

---

## 1. What it is

PETRI-seq is **SPLiT-seq adapted to bacteria**: fixed, permeabilised cells are the
compartments, and three rounds of split-pool barcoding (one RT, two ligations, 96 wells
each, 96³ ≈ 0.9 M combinations) label every cDNA with a cell barcode combination (BC)
plus a 7-nt UMI. ~10,000 cells are then lysed per sub-library. 🟢 (main text, Methods)

| | New thing here | Builds on |
|---|---|---|
| 1 | **Bacterial cell prep**: pellet, 4 % formaldehyde 16 h, 50 % ethanol, lysozyme (E. coli) or lysostaphin (S. aureus), in situ **DNase I** then heat-inactivation at 50 °C (not 70 °C, to keep cells intact) 🟢 | SPLiT-seq (eukaryotic fixation/permeabilisation) |
| 2 | **Round-1 RT with random hexamers only** — bacterial mRNA has no poly(A); all 96 round-1 oligos end in `NNNNNN` 🟢 | SPLiT-seq used dT + random mix ([reverse transcription](../ref/concepts/reverse-transcription.md)) |
| 3 | **Short (16-nt) ligation linkers** instead of SPLiT-seq's 30-nt ones, so the barcode read fits in 58 cycles (75-cycle kit, ~50 % cheaper); barcodes overlap the linker by one constrained base (S or W) 🟢 | [SPLiT-seq](../split-seq__10.1126+science.aam8999/01_split-seq.md) splint ligation |
| 4 | **Two blocking oligos per round** (a hairpin + the linker complement) instead of one 🟢 | SPLiT-seq blocking strand |
| 5 | **AMPure clean-up of lysate** instead of biotin/streptavidin pull-down; then **second-strand synthesis** (NEBNext) rather than template switching 🟢 | SPLiT-seq (streptavidin + TSO) |
| 6 | Nextera XT tagmentation, then PCR with **NEB i50x** on the barcode end and **Nextera N7xx** on the Tn5 end 🟢 (that i50x carries TruSeq Read 1: 🟡) | [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) |

## 2. Oligos

### 2.1 Barcode plates — 🟢 Supplementary Table 3 (`MOESM3_ESM.xlsx`, sheet `Table_S3`)

96 of each, IDT notation as given (`/5Phos/` = 5′ phosphate). Template and the first and
last well of each plate, verbatim:

```
Round 1 RT          /5Phos/GCCAGA <bc1:7> NNNNNN                              19 nt
  BC 1              /5Phos/GCCAGACAGAGAANNNNNN
  BC 96             /5Phos/GCCAGAGTAGACTNNNNNN
Round 2 ligation    /5Phos/GCTTCGC <bc2:7> CCTCCTAC                           22 nt
  BC 1              /5Phos/GCTTCGCGACCTTACCTCCTAC
  BC 96             /5Phos/GCTTCGCGTGCTCTCCTCCTAC
Round 3 ligation    AGAATACACGACGCTCTTCCGATCT <UMI:NNNNNNN> <bc3:7> GGTCCTTG     47 nt
  BC 1              AGAATACACGACGCTCTTCCGATCTNNNNNNNGTGTGAAGGTCCTTG
  BC 96             AGAATACACGACGCTCTTCCGATCTNNNNNNNGCATACTGGTCCTTG
```

Round 3 is **not** phosphorylated (it is the last piece; its 5′ end becomes the library
end). A **biotinylated round-3 set** is priced as an alternative in Supp. Table 1c, but
the Methods use AMPure, not streptavidin. 🟢

### 2.2 Single-tube oligos — 🟢 Supplementary Table 2 (`MOESM1_ESM.pdf`)

The `.txt` twin scrambles this table's columns (see §8); name–sequence pairs below were
read from a layout-preserving extraction of the same PDF and agree with the upstream page
for every oligo upstream lists.

```
SB83   Round 2 Linker Oligo        STCTGGCGTAGGAGGW
SB84   Round 2 Blocking 1          GCCAGASACGTTAGGCAGGACCTAACGT
SB85   Round 2 Blocking 2          WCCTCCTACGCCAGAS
SB80   Round 3 Linker Oligo        GCGAAGCCAAGGACCW
SB81   Round 3 Blocking 1          GCTTCGCTGCAATCGGACCTCGATTGCA
SB82   Round 3 Blocking 2          WGGTCCTTGGCTTCGC

SB14   Template Switch Oligo                         AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
SB15   PCR Primer for After Template Switch (+SB86)  AAGCAGTGGTATCAACGCAG
SB86   PCR Primer for After Template Switch (+SB15 or SB13)  AGAATACACGACGCTCTTCC
SB94   Hexamer for RT Without Barcoding              AGAATACACGACGCTCTTCCGATCTNNNNNN

SB10   rpsB-specific RT primer             /5Biosg/ACAGACATGTGCTCTTCCAGCTGAGAACGGCCTTCAC
SB12   SB10-binding rpsB qPCR primer (+SB13)       ACAGACATGTGCTCTTCCAG
SB13   rpsB qPCR primer (+SB12, SB86 or SB115)     ATACCAACTCTGATCCGGAC
SB5    rpsB qPCR primer                            AAACCGTTCGTCAGTCCATC
SB6    rpsB qPCR primer                            ATGTCTTTGATACCGCCCAG

SB110  rpsB RT primer, ligation test       /5Phos/GCCAGAGGCCAGGAATGAGAACGGCCTTCAC
SB111  Linker for ligation test (30 nt)    TTCCTGGCCTCTGGCGTAGGAGGTGGAAGA
SB113  Ligated primer, ligation test       AGAATACACGACGCTCTTCCACCTCCTAC
SB114  ligation-test positive control      AGAATACACGACGCTCTTCCACCTCCTACGCCAGAGGCCAGGAATGAGAACGGCCTTCAC
SB115  qPCR primer, ligation test (+SB13)  AGAGGCCAGGAATGAGAAC
```

(The table's descriptions for SB5/SB6 read "with SB5" / "with SB6" — each naming
itself, evidently a typo; the Methods state that SB5 and SB6 were used together as the
qPCR pair for genomic DNA / random-hexamer cDNA. 🟢)

### 2.3 Commercial primers (not listed in the paper's tables)

The paper names the kits only: **Nextera Index Kit v2 Set A N70x** and **NEB E7600 i50x**
🟢. Sequences 🟡 (upstream and `lib/`):

```
NEB i50x      AATGATACGGCGACCACCGAGATCTACAC <i5:8> ACACTCTTTCCCTACACGACGCTCTTCCGATCT
Nextera N7xx  CAAGCAGAAGACGGCATACGAGAT <i7:8> GTCTCGTGGGCTCGG
```

## 3. How the oligos interlock — 🟡 (computed)

**Barcode constraints.** Round-1 barcode position 1 is `C` in wells 1–48 and `G` in
49–96 (= `S`); rounds 2 and 3 end in `A` (1–48) or `T` (49–96) (= `W`). All 96 barcodes
of each round are distinct, minimum pairwise Hamming distance 3. These constrained bases
are exactly the degenerate bases of the linkers — the one-base overlap of barcode and
linker the Extended Data Fig. 1e legend mentions. 🟢 (upstream states the S/W rule too)

**Linkers (splints).**

- `revcomp(SB83)` = `W` · `CCTCCTAC` · `GCCAGA` · `S`: it pairs with the last barcode base
  and 8-nt 3′ constant of the round-2 oligo **and** the 6-nt 5′ constant and first barcode
  base of the round-1 RT primer — a 16-nt bridge, 9 + 7.
- `revcomp(SB80)` = `W` · `GGTCCTTG` · `GCTTCGC`: pairs with the round-3 oligo's 3′ end
  (last barcode base + 8-nt constant) and the round-2 oligo's 5′ `GCTTCGC` (7 nt).
- Each ligation joins a 3′-OH (incoming barcode oligo) to the 5′-phosphate of the piece
  already on the cDNA (round-1 primer, then round-2 oligo) — hence `/5Phos/` on rounds 1
  and 2 only.

**Blockers.**

- **SB85 = revcomp(SB83)** and **SB82 = revcomp(SB80)** exactly: the second blocker is the
  full complement of the linker, sequestering free linker.
- **SB84** = `GCCAGAS` + hairpin (stem `ACGTTAGG`, loop `CAGGA`, stem `CCTAACGT`); its 7-nt
  single-stranded 5′ arm is the round-1 primer's `GCCAGA`+`S`, i.e. it pairs with the
  round-1-side arm of SB83. **SB81** = `GCTTCGC` + hairpin (stem `TGCAATCG`, loop `GACCT`,
  stem `CGATTGCA`); its arm is the round-2 oligo's 5′ end, pairing with SB80. The Methods
  anneal SB84 and SB81 "to form an intramolecular hairpin" 🟢. Why a hairpin decoy is
  needed in addition to the full complement is not explained 🔴.
- The round-3 blocking mix also carries **EDTA** (600 µL 0.5 M in 1.2 mL), which stops the
  ligase 🟢; the round-2 blocking mix has none.

**Library primer sites.**

- Round-3 5′ constant `AGAATACACGACGCTCTTCCGATCT` (25 nt) is **not** a substring of
  `illumina.TRUSEQ_READ1`: they share only the 3′ **21 nt** `TACACGACGCTCTTCCGATCT`; the
  first 4 nt (`AGAA`) differ from TruSeq's `TCCC`. The NEB i50x primer
  (`illumina.nebnext_i5_primer` = P5 + i5 + full `TRUSEQ_READ1`) anneals over those 21 nt
  and overwrites the first 4 in the product. The same `AGAA…` start is SPLiT-seq's round-3
  design.
- **SB86** = the round-3 5′ constant's first 20 nt; **SB94** = round-3 constant + N6 (an
  unbarcoded hexamer RT primer carrying the same handle).
- Template-switch arm: **SB14** = full `rt.SMART_HANDLE` (23 nt) + `GAAT` + `rGrGrG`;
  **SB15** = `SMART_HANDLE[:20]`.
- Round-1 constant `GCCAGA` occurs inside SPLiT-seq's round-1 handle
  `AGGCCAGAGCATTCG`; whether it was derived from it is not stated.

**Ligation test oligos (Extended Data Fig. 2l).** SB114 = SB113 + SB110 exactly (a
pre-ligated positive control). `revcomp(SB111)` = SB113[-15:] + SB110[:15] (30-nt
splint), while SB83 bridges the same junction with 16 nt (SB113 ends `…ACCTCCTAC` = W +
round-2 constant; SB110 starts `GCCAGAG` = round-1 constant + S). SB115 sits in SB110 at
offset 3. SB10 and SB110 share the rpsB-binding 3′ end `…GAGAACGGCCTTCAC`; SB12 = the
first 20 nt of SB10.

## 4. Step by step

Volumes for the optimised "4×" protocol (Experiment 2.01); the "1×" values are in the
paper. All 🟢 from Methods unless marked.

**Annealing (before day 2).** Round-2 plate: SB83 + each round-2 oligo; round-3 plate:
SB80 + each round-3 oligo; 95 °C 3 min, ramp −0.1 °C/s to 20 °C. SB84 and SB81 heated to
94 °C and cooled slowly to form hairpins.

**Cell preparation.**

1. Pellet culture (5,525×g, 2 min, 4 °C), resuspend in ice-cold **4 % formaldehyde/PBS**,
   rotate 16 h at 4 °C.
2. Wash in PBS-RI (SUPERase·In), **50 % ethanol** in PBS-RI, two washes.
3. **Permeabilise** 15 min RT: 100 µg/mL lysozyme (E. coli) or 40 µg/mL lysostaphin
   (S. aureus) in TEL-RI (100 mM Tris pH 8, 50 mM EDTA, RNase inhibitor).
4. **DNase I** in situ, 30 min RT; Stop Solution, **50 °C 10 min**. Wash, count.

**Split-pool barcoding.**

5. **Round 1 — RT in situ**: 2 µL of 10 µM round-1 primer per well + 8 µL mix (Maxima H
   Minus RT, dNTPs, RNase inhibitor; 3 × 10⁷ cells in 960 µL). 50 °C 10 min, then a
   stepped ramp 8 → 15 → 20 → 30 → 42 °C (6 min) → 50 °C 16 min. The initial 50 °C step
   before primer annealing is presumably to denature RNA structure 🟡. Each cDNA now
   begins `5′-P-GCCAGA · bc1 · NNNNNN · cDNA…` ([reverse transcription](../ref/concepts/reverse-transcription.md)).
6. Pool; (2.01 only) Tween-20 0.04 % then 0.01 %; spin 10,000×g 20 min.
7. **Round 2 — ligation**: cells in T4 ligase buffer + T4 DNA ligase + BSA + RNase
   inhibitor, 5.76 µL into 2.24 µL annealed round-2 adaptor per well, **37 °C 30 min**.
   Add 2 µL blocking mix (SB84 + SB85), 37 °C 30 min. Pool.
8. **Round 3 — ligation**: 8.51 µL cells/ligase into 3.49 µL annealed round-3 adaptor,
   37 °C 30 min. Add 10 µL round-3 blocking mix (SB81 + SB82 + **EDTA**). Pool.
9. Wash (TEL-RI ± 0.01 % Tween), count; **aliquots of ~10,000 cells** into lysis buffer
   (50 mM Tris pH 8, 25 mM EDTA, 200 mM NaCl) + proteinase K, **55 °C 1 h**. Store −80 °C.

After step 8 the first-strand cDNA reads, 5′→3′ 🟡 (computed from the oligos):

```
5'- AGAATACACGACGCTCTTCCGATCT · UMI(7) · bc3(7) · GGTCCTTG · GCTTCGC · bc2(7) · CCTCCTAC · GCCAGA · bc1(7) · NNNNNN · cDNA (antisense to RNA) -3'
```

**Library preparation.**

10. **AMPure XP 1.8×** on the lysate, elute 20 µL.
11. **Second strand**: NEBNext Second Strand Synthesis (E6111S), 16 °C 2.5 h; AMPure 1.8×.
    (Mechanism — RNase H nicking + DNA Pol I — is the kit's, not stated in the paper 🟡.)
12. **Tagment** with Nextera XT (25 µL TD, 20 µL dscDNA, 5 µL ATM), neutralise (NT).
13. **PCR**: N70x (Nextera v2 Set A) + i50x (NEB E7600) in NPM; 8 cycles, qPCR side
    reaction to find the exponential phase, then +8 (2.01) or +11 (1.06SaEc, 1.10) cycles
    = 16 or 19 total. AMPure 1×. Only fragments with **round-3 handle at one end and an
    s7 Tn5 end at the other** amplify; s5-tagmented ends have no primer, which the
    authors note costs **2-fold** in capture 🟢 ("Future directions for optimization")
    — see [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). Undigested genomic DNA
    lacks the round-3 handle and is not amplified 🟢.

### Variants in the paper (Supp. Table 4)

| Variant | Difference | Experiment |
|---|---|---|
| "Original" (1×) | 1× ligation adaptor/blocker concentrations, no detergent | 1.06SaEc, 1.10, 1.06 |
| **Optimised (4×)** | 4× annealed adaptors and blockers, Tween-20 during washes | 2.01 |
| Round-3 test | half plate 1× / half 4× round-3 adaptor (4× gives 2.7× more mRNA) | 1.12 |
| **Template switching** instead of second strand | AMPure-purified lysate + **SB14** TSO, Maxima H Minus, betaine, 42 °C 90 min, 85 °C 5 min; AMPure; 10 cycles PrimeSTAR GXL with **SB86 + SB15**; then Nextera XT as usual. Second-strand synthesis gave ~3× more mRNAs/cell, so it was kept ([template switching](../ref/concepts/template-switching.md)) | 1.08 (also 10× RT primer, 5×/3.5× ligation, 2 h RT) |
| "RT clean-up" | a second RT on purified lysate (no TSO) before second strand; no yield change | 1.06SaEc-replicate |

The template-switch product would be `5′-SB86 site … bc1 · N6 · cDNA · CCC-complement ·
revcomp(SB14 DNA part)-3′` before tagmentation 🟡 (standard mechanism; not drawn in the
paper).

## 5. Final library — 🟡 (assembled from the oligos above)

```
5'- P5 · i5(8) · TruSeq Read 1 (33) · UMI(7) · bc3(7) · GGTCCTTG · GCTTCGC · bc2(7) · CCTCCTAC · GCCAGA · bc1(7) · cDNA (starts with the N6 primer) · ME' · s7' · i7'(8) · P7' -3'
```

`ME'` = `nextera.ME_RC` (`CTGTCTCTTATACACATCT`), `s7'` = `nextera.S7_RC`, `P7'` =
`illumina.P7_RC`. The first 4 nt of the round-3 oligo (`AGAA`) are replaced by TruSeq's
`TCCC` from the i50x primer.

## 6. Sequencing

🟢 NextSeq 500/550 High Output v2.5, **75-cycle kit**: **Read 1 58 cycles** (UMI +
barcodes), **Read 2 17 cycles** (cDNA), **Index 1 8**, **Index 2 8**.

- Read 1 (TruSeq Read 1 primer, `illumina.TRUSEQ_READ1`) reads UMI(7) · bc3(7) ·
  15-nt linker `GGTCCTTGGCTTCGC` · bc2(7) · 14-nt linker `CCTCCTACGCCAGA` · bc1(7) = **57
  nt** 🟡 (computed); the 58th base is the first hexamer base. The demultiplexing
  description (barcode 3 + downstream linker `GGTCCTTGGCTTCGC`, 21 positions; barcode 2
  + linker, 20) matches this layout 🟢.
- Read 2 (Nextera Read 2 primer = `nextera.READ2_PRIMER`) reads the cDNA from the Tn5
  end, **sense to the RNA** 🟡 (the first strand is antisense; read 2 copies the
  opposite strand). The pipeline aligns read 2 only and counts sense alignments 🟢;
  this also explains the mild 3′ coverage bias (Supp. Fig. 2: tagmentation must fall
  5′ of the hexamer priming site) 🟢.
- Index 1 (i7) uses the Nextera index-1 primer `nextera.INDEX1_PRIMER`; index 2 (i5) the
  TruSeq index-2 primer — 🟡 (upstream; the paper does not name custom primers, and a
  mixed TruSeq-i5 / Nextera-i7 library needs both).

## 7. Upstream scg_lib_structs page vs. the paper

| Point | Upstream | Paper | Verdict |
|---|---|---|---|
| Round 1/2/3 oligo templates | `/5Phos/ GCCAGA[7]NNNNNN`, `/5Phos/ GCTTCGC[7]CCTCCTAC`, `…GATCT[7 UMI][7]GGTCCTTG` | identical (Supp. Table 3) | agree |
| SB80–SB85 | as §2.2 | identical (Supp. Table 2) | agree |
| S/W rule | first R1 base S, last R2/R3 base W | Table 3 obeys it (computed) | agree |
| Second strand | NEBNext Second Strand Synthesis | same | agree |
| Amplifiable product | only round-3-end + s7 | same, and the paper quantifies the 2-fold loss | agree |
| PCR primers | NEB i50x + Nextera N7xx | same kits named | agree (sequences 🟡) |
| Read 1 length | "at least 57 cycles" | 58 cycles used | consistent |
| Read 2 / index | 17 cycles cDNA, 8 + 8 | same | agree |
| Hairpin annealing of SB84/SB81, EDTA in round-3 block, 1× vs 4×, TSO variant | not mentioned | in Methods | upstream omits |
| Source | authors' website (2024) | paper (2020) | website not checked 🔴 |

No disagreement in sequence was found.

## 8. Open questions

- 🔴 Why a **hairpin decoy** (SB84, SB81) in addition to the full-complement blocker
  (SB85, SB82); which of the two does the work is not tested in the paper.
- 🔴 Why the round-3 oligo begins `AGAA…` rather than with full TruSeq Read 1 (shared
  with SPLiT-seq; probably inherited).
- 🔴 The `GAAT` between the SMART handle and `rGrGrG` in SB14 is not explained.
- 🟡 Whether the authors' later website protocol (upstream's source) changed reagents or
  concentrations relative to the paper — not fetched.
- SB5/SB6 pairing descriptions in Supp. Table 2 are self-referential (typo); resolved by the Methods, which use SB5 + SB6 together 🟢.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` got the PMC full text, all 16 Springer ESM files and the upstream
page. The three PMC `NIHMS…supplement` files it reports as **saved** are the PMC
download-gate HTML (≈21 kB each), not the PDF/XLSX/GZ — they should have been listed as
`(manual)`; the Springer copies cover them. `tools/doctext.py` flattened Supp. Table 2 so
that names, descriptions and sequences fell out of order (SB13 and SB14 were run together
into one line); a `pdftotext -layout` pass fixed it. `tools/scrape_primers.py` ranked the
round-1/2 barcode table first; it merged SB13 and SB14 into one 47-nt "oligo" because of
that text. `--find` reported **no location** for SB84 (`GCCAGASACGTTAGGCAGGACCTAACGT`) and
for full round-3 oligos (`AGAATACACGACGCTCTTCCGATCTNNNNNNN…GGTCCTTG`), although both
occur verbatim (plain `grep`: `MOESM1_ESM.pdf.txt` / `upstream_PETRI-seq.html.txt`, and
`MOESM3_ESM.xlsx.txt` line 4); `--find` on the 15-nt barcode tail found the latter.
