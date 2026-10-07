# Quartz-Seq and Quartz-Seq2 — poly(A)-tagging whole-transcript amplification

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note. Claims taken only from
> the upstream scg_lib_structs page are 🟡 (secondary source).

This note covers the whole family: **Quartz-Seq** (2013, one tube per cell, full-length
cDNA, TruSeq ligation library; §A) and **Quartz-Seq2** (2018, cell-barcoded and UMI-tagged
RT primer, pooled 384 / 1536 cells, 3'-end library; §B). They share one amplification
handle and one WTA chemistry, so the oligos interlock across both (§C).

**Quartz-Seq** — Sasagawa Y, Nikaido I, Hayashi T, Danno H, Uno KD, Imai T, Ueda HR.
"Quartz-Seq: a highly reproducible and sensitive single-cell RNA sequencing method, reveals
non-genetic gene-expression heterogeneity." *Genome Biology* 14:R31 (2013).
doi:[10.1186/gb-2013-14-4-r31](https://doi.org/10.1186/gb-2013-14-4-r31) · PMID 23594475 ·
PMC4054835. Data: GEO GSE42268.

**Quartz-Seq2** — Sasagawa Y, Danno H, Takada H, Ebisawa M, Tanaka K, Hayashi T, Kurisaki A,
Nikaido I. "Quartz-Seq2: a high-throughput single-cell RNA-sequencing method that
effectively uses limited sequence reads." *Genome Biology* 19:29 (2018).
doi:[10.1186/s13059-018-1407-3](https://doi.org/10.1186/s13059-018-1407-3) · PMID 29523163 ·
PMC5845169. Data: GEO GSE99866, DRA002954.

Cited by upstream for the suppression-PCR idea (not fetched, not read): Plessy C *et al.*
"Linking promoters to functional transcripts in small samples with nanoCAGE and CAGEscan."
*Nat Methods* 7:528 (2010). doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470).

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/quartz-seq-family__10.1186+gb-2013-14-4-r31/`, never committed):

| File | What | Used for |
|---|---|---|
| `gb-2013-14-4-r31_PMC4054835.html.txt` | Quartz-Seq full text (PMC) | **Table 1 (all Quartz-Seq oligos)**, WTA and LIMprep methods, design rationale |
| `gb-2013-14-4-r31_supp_gb-2013-14-4-r31-S1.PDF.txt` | Additional file 1, figure legends | suppression-primer variants A/B/C/D (Figs S3, S4) |
| `gb-2013-14-4-r31_supp_gb-2013-14-4-r31-S2.DOCX.txt` | Additional file 2, supplementary note | polymerase choice, single-tube carry-over |
| `s13059-018-1407-3_PMC5845169.html.txt` | Quartz-Seq2 full text (PMC) | five-step overview, WTA, adaptor and library methods, sequencing |
| `s13059-018-1407-3_supp_13059_2018_1407_MOESM5_ESM.xlsx.txt` | Additional file 5, Table S4 | **every Quartz-Seq2 oligo**, 384 v3.1 + 1536 v3.2 RT primers |
| `s13059-018-1407-3_supp_13059_2018_1407_MOESM6_ESM.xlsx.txt` | Additional file 6, Table S5 | TdT buffer compositions (T55 etc.) |
| `s13059-018-1407-3_supp_13059_2018_1407_MOESM1_ESM.pdf.txt` | Additional file 1, supplementary note + legends | why exonuclease I was dropped, ligation vs Nextera |
| `upstream_Quartz-seq_family.html.txt` | scg_lib_structs page | second source, checked against the above (§D) |

Not read / not usable:

| File or URL | Why |
|---|---|
| `gb-2013-14-4-r31_supp_gb-2013-14-4-r31-S3.PDF` | image-only PDF; text twin is empty (Figure S7, scatter plots — not chemistry) |
| `…-S4.XLS`, `…-S7.XLS`, `…-S8.XLS` | legacy Excel, no `.txt` twin made. S8 holds qPCR primers (seen with `strings`), S7 the sample table; not library chemistry |
| nanoCAGE paper, doi:10.1038/nmeth.1470 (PMC2906222) | cited only; not fetched |
| RIKEN step-by-step protocols (linked from the upstream page) | not fetched |

---

## 1. What it is, and what is new

Quartz-Seq is a **single-tube, purification-free** whole-transcript amplification (WTA)
for one cell in 0.4 µL lysis buffer, built on **poly(A) tailing** of first-strand cDNA by
terminal transferase rather than on template switching. 🟢 (main text, "Whole-transcript
amplification for single-cell Quartz-Seq"; the paper places itself against both families
in its Background.) Five steps in one tube: RT with an oligo-dT primer → exonuclease I
→ poly(A) tailing → second strand from a poly(dT) "tagging primer" → PCR with a single
"suppression PCR primer". 🟢

The point of the method is **killing the RT-primer byproduct** that poly(A)-tailing WTA
had always made (RT primer gets tailed, primes second strand, gives short junk). Three
fixes, used together 🟢 (main text, Figs S2–S4):

| | Fix | Mechanism |
|---|---|---|
| 1 | Minimal RT primer, then **exonuclease I** digests what is left | removes single-stranded primer (incompletely) |
| 2 | **Restricted TdT time** (50 s) | leftover tailed primers stay short |
| 3 | **Suppression PCR** with one primer (the "B primer") | a short template with the same handle at both ends folds into a "pan-like" hairpin and outcompetes primer binding; long cDNA does not |

Quartz-Seq2 keeps the WTA chemistry but moves the cell identity into the RT primer
(14- or 15-nt cell barcode + 8-nt UMI), so cDNA from up to **1536 wells is pooled right
after RT**; column purification of the pool replaces exonuclease I, which works because
the RT primer is kept short (v3.1 is a 73-mer, about the length of the 70-mer Quartz-Seq
primer; a 126-mer gave heavy byproduct) 🟢 (MOESM1 Supplemental Note). The v3.2 primer is
74 nt 🟡 (counted from Table S4). It also improves
tagging efficiency 3.6-fold (T55 TdT buffer + an "Increment" ramp into second-strand
synthesis), uses 4× less RT enzyme (RT25), and builds a **3'-end library** whose Read 1
is a custom primer reading barcode + UMI directly. 🟢 (main text, "Improvement of poly(A)
tagging efficiency")

Concepts: oligo-dT priming as in [reverse transcription](../ref/concepts/reverse-transcription.md);
the poly(A)-tail route here is the alternative to
[template switching](../ref/concepts/template-switching.md) for putting a handle on the
cDNA 3' end. No Tn5 in the main protocols (Quartz-Seq2 mentions a Nextera XT library as an
alternative; see [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)).

---

## A. Quartz-Seq (2013)

### A.2 Oligos

🟢 Verbatim from Table 1 of the main text (`*` = phosphorothioate bond, as the table
footnote says; `(NH2)` = 5' amino; `(5'-phosphate)` as written).

```
RT primer (WTA)      TATAGAATTCGCGGCCGCTCGCGATAATACGACTCACTATAGGGCGTTTTTTTTTTTTTTTTTTTTTTTT
Tagging primer       TATAGAATTCGCGGCCGCTCGCGATTTTTTTTTTTTTTTTTTTTTTTT
Suppression primer   (NH2)-GTATAGAATTCGCGGCCGCTCGCGAT
TRSU                 AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T
TRSI-2               (5'-phosphate)GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTT*G
TRSI-4               (5'-phosphate)GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTT*G
TRSI-5               (5'-phosphate)GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTT*G
TRSI-6               (5'-phosphate)GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTT*G
TRSI-7               (5'-phosphate)GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTT*G
TRSI-12              (5'-phosphate)GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTT*G
TPC1                 AATGATACGGCGACCACCGA*G
TPC2                 CAAGCAGAAGACGGCATACGA*G
RT primer (qPCR)     TATAGAATTCGCGGCCGCTCGCGATAATACGACTCACTATAGGGCGTTTTTTTTTTTTTTTTTTTTTTTT
```

The suppression-primer candidates tested (Additional file 1, Figs S3c / S4c) 🟢:
A `TATAGAATTCGCGGCCGCTCGCGAT`, **B `GTATAGAATTCGCGGCCGCTCGCGAT` (used)**,
C `CTATAGAATTCGCGGCCGCTCGCGAT`, D `TGTATAGAATTCGCGGCCGCTCGCGAT`. The figure legends give
them without the 5' amino group that Table 1 shows on the suppression primer.

### A.3 How they interlock — 🟡 (computed)

One 25-nt handle, called **M** here (the paper's "PCR target region (M)", Fig. 1 legend 🟢):
**`TATAGAATTCGCGGCCGCTCGCGAT`** (= primer A).

- **RT primer (WTA)** = M + `AATACGACTCACTATAGGGCG` (21 nt) + T₂₄ — 70 nt. The middle
  contains the **T7 promoter** `TAATACGACTCACTATAG` starting at the last T of M (offset 24),
  followed by `GGCG`; this is what Quartz-Chip uses for cRNA labelling. No `VN` anchor.
- **Tagging primer** = M + T₂₃ — 48 nt (the T of M's 3' end makes the run 24 T in a row).
- **Suppression primer** = `G` + M — 26 nt. So after second-strand synthesis **both ends of
  every product carry M**, one primer amplifies everything, and the 5' extra G is the only
  difference between primers A and B.
- **TPC1** = `illumina.P5[:21]`; **TPC2** = `illumina.P7[:22]`; each has a phosphorothioate
  before its last base.
- **TRSU** = `illumina.TRUSEQ_P5_FULL` (P5 + TruSeq Read 1, 58 nt) with a phosphorothioate
  before the 3' T overhang.
- **TRSI-n** = `illumina.INDEX1_PRIMER` (33 nt) + 6-nt index + `illumina.P7_RC` — 63 nt;
  starts with `illumina.STEM` (`GATCGGAAGAGC`), which pairs with the `GCTCTTCCGATC` just
  before TRSU's 3' T: a **12-bp TruSeq Y-adapter stem**. The index is written in **read
  orientation** (the Index 1 primer is TRSI's own 5' end, so Index 1 reads the bases as
  written): TRSI-2 `CGATGT`, -4 `TGACCA`, -5 `ACAGTG`, -6 `GCCAAT`, -7 `CAGATC`,
  -12 `CTTGTA` — the TruSeq LT indices of the same numbers.

Why the byproduct is suppressed 🟡 (mechanism from Fig. S4 legend 🟢, sequence check
computed): a leftover RT primer that got tailed and then copied by the tagging primer
becomes `M-T7-T24 · A-tail · (complement) · M'` — a ~150-nt molecule whose two ends are
M and its complement, so it snaps shut into a hairpin before the suppression primer can
bind. Real cDNA has the same ends but is long enough for the primer to win.

### A.4 Step by step

All 🟢 from Methods; one tube from lysis to amplified cDNA, no purification.

1. **Sort** one cell (Hoechst-gated FACS) into **0.4 µL lysis buffer (0.5 % NP-40)** with
   an air bubble added by pipette (improves success rate, Fig. S13), on a 0 °C rack.
2. **Priming**: + 0.8 µL priming buffer (1.5× PCR buffer with MgCl₂, "41.67 pmol/l" RT
   primer, 4 U/µL RNasin Plus, 50 µM dNTPs); 70 °C 90 s, 35 °C 15 s.
3. **RT**: + 0.8 µL (1× PCR buffer, 25 U/µL SuperScript III, 12.5 mM DTT); 35 °C 5 min,
   45 °C 20 min, 70 °C 10 min. First strand: `5'-M-T7-GGCG-T24-cDNA…-3'`.
4. **Primer digestion**: + 1 µL exonuclease I (1.5 U/µL); 37 °C 30 min, 80 °C 10 min.
5. **Poly(A) tailing + RNA removal**: + 2.5 µL (1× PCR buffer, 3 mM dATP, 33.6 U/µL terminal
   transferase, 0.048 U/µL RNase H); **37 °C 50 s** (deliberately short), 65 °C 10 min.
6. **Second strand**: + 23 µL (1.09× MightyAmp Buffer v2, "70 pmol/l" tagging primer,
   0.054 U/µL MightyAmp polymerase); 98 °C 130 s, 40 °C 1 min, 68 °C 5 min. The tagging
   primer's T₂₄ anneals to the new poly(A); product: `M-T24-…(cDNA')…-T24-GGCG-T7'…`
   with M at one end and the RT-primer handle at the other. 🟡 (orientation)
7. **Suppression PCR**: + 25 µL (1× MightyAmp Buffer v2, 1.9 µM suppression primer);
   21 × (98 °C 10 s, 65 °C 15 s, 68 °C 5 min), then 68 °C 5 min. MinElute or AMPure XP.
8. **LIMprep library** (ligation-based, KAPA reagents, home-made TruSeq adapters):
   20 ng cDNA in 130 µL TE → Covaris S220 (duty 10 %, 175 W, 100 cycles/burst, 600 s) →
   Zymo DCC-5 → KAPA end repair 20 °C 30 min → DCC-5 → A-tailing 30 °C 30 min → DCC-5 →
   ligation with 10 pmol adapter, 20 °C 15 min → two rounds of AMPure XP with a high-salt
   "binding support buffer" (1 M NaCl, 20 mM MgCl₂) to remove adapter dimer → KAPA HiFi PCR
   with 17.5 pmol each TPC1 / TPC2: 98 °C 45 s, then typically 10–12 × (98 °C 15 s,
   60 °C 30 s, 72 °C 30 s) → AMPure XP.
9. **Adapter annealing**: TRSU + one TRSI, 100 µM each in adaptor buffer (10 mM Tris pH 7.8,
   0.1 mM EDTA, 50 mM NaCl), 95 °C 2 min then −0.5 °C per cycle × 170; dilute to 10 µM,
   store at −80 °C.

### A.5 Final library — 🟡 (assembled from the oligos)

A standard TruSeq single-index library of sheared full-length cDNA:

```
5'- P5 · TruSeq Read 1 · <cDNA fragment> · TruSeq Read 2' · i7' · P7' -3'
```

with `P5 · TruSeq Read 1` = TRSU and `TruSeq Read 2' · i7' · P7'` = the TRSI-n top strand
(`AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` · index as written · `ATCTCGTATGCCGTCTTCTGCTTG`).
Fragments from the two cDNA ends carry M (plus T7 or T₂₄) next to the insert; the paper
counts **7.68 ± 0.66 % of reads as WTA adaptor** 🟢 and trims them with Trimmomatic 🟢 —
so end fragments do reach the sequencer (see §D on upstream's claim).

### A.6 Sequencing

🟢 HiSeq 1000/2000, **paired-end**; comparisons use "PE sequences of 50 bp" and
single-end 50-bp read subsets. Standard TruSeq read primers (`illumina.TRUSEQ_READ1`,
`illumina.TRUSEQ_READ2`, `illumina.INDEX1_PRIMER`); 6-nt Index 1. 🟡 (primers not named
in the paper; implied by the TruSeq adapters)

---

## B. Quartz-Seq2 (2018)

### B.2 Oligos

🟢 Verbatim from Table S4 (`MOESM5_ESM.xlsx`). Sheet "Quartz-Seq2 related oligo DNA":

```
Tagging primer    TATAGAATTCGCGGCCGCTCGCGATTTTTTTTTTTTTTTTTTTTTTTT
gM_primer         GTATAGAATTCGCGGCCGCTCGCGAT
rYshapeP5         GATCGGAAGAGCGTCGTGTA
rYshapeP7LT06     CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT12     CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT05     CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT19     CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT02     CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT04     CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT07     CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
rYshapeP7LT16     CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
P5-gMac_hybrid    aatgatacggcgaccaccgagatctacattgtatagaattcgcggccgctcgcgaT*A*c
TPC2              CAAGCAGAAGACGGCATACGA*G
Read1DropQuartz   ACATTGTATAGAATTCGCGGCCGCTCGCGATAC
```

(P5-gMac_hybrid is given in lower case with upper-case `T*A*` before the last base, and
"DNA, * S" as its type — `*` = phosphorothioate. Case as in the source; it carries no
meaning the table explains. 🟢 / 🔴)

RT primers, two sets (sheets "384 RT primer with SeqLv14mer" and "1536 RT primer with
SeqLv15mer"); first rows verbatim 🟢:

```
MDRT001     TATAGAATTCGCGGCCGCTCGCGATACCATATTCCTGGTGGNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTT    v3.1, 384
MDRT384     TATAGAATTCGCGGCCGCTCGCGATACGGCGAGTACCTCTTNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTT
eMDRT0001   TATAGAATTCGCGGCCGCTCGCGATACATCAATCCTATCTGCNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTT   v3.2, 1536
eMDRT1536   TATAGAATTCGCGGCCGCTCGCGATACGAGTTGGCTACTATANNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTT
```

General form: `TATAGAATTCGCGGCCGCTCGCGATAC <cbc> NNNNNNNN T24`, OPC purification grade.

### B.3 How they interlock — 🟡 (computed)

- **Every RT primer** = M + `AC` + cell barcode + N₈ (UMI) + T₂₄. Parsed for all rows of
  both sheets: v3.1 = **384 unique 14-nt** barcodes (73 nt), v3.2 = **1536 unique 15-nt**
  barcodes (74 nt); each parsed barcode equals the table's "cell barcode sequence" column;
  minimum pairwise Hamming distance 5 in both sets (the paper designs for Sequence–Levenshtein
  distance ≥ 5 🟢). No T7 promoter any more, no `VN` anchor. (The v3.2 sheet's caption still
  says "14 mer", its title says 15; the sequences are 15. 🟢 / 🟡)
- **Tagging primer** and **gM_primer** are identical to Quartz-Seq's tagging and
  suppression ("B") primers; the gM_primer is listed **without** the 5' amino group.
- **The `AC` is the discriminator.** After WTA the RT-primer end of a cDNA reads
  `G·M·AC·<cbc>·<umi>·T24`, the tagging-primer end `G·M·T24`. **P5-gMac_hybrid** =
  `illumina.P5[:28]` + `TT` + gM_primer + `AC` (58 nt): its 3' end is M + `AC`, so it
  extends only from the RT-primer (3') end; at a tagging-primer end its last two bases
  face `TT` and mismatch. The two phosphorothioates before the 3' end plausibly stop KAPA
  HiFi's 3'→5' proofreading from trimming that mismatch away. 🟡 (inferred; the paper does
  not explain the design)
- **P5 is truncated in P5-gMac_hybrid**: the first 28 nt of P5 are followed by `T`, where
  P5's 29th base is `C`. Whether this matters for cluster generation is not discussed. 🔴
- **Read1DropQuartz** = `ACATT` + gM_primer + `AC` (33 nt) = the last 33 nt of
  P5-gMac_hybrid. It ends exactly at the `AC`, so **Read 1 starts on the first base of the
  cell barcode**: 14 + 8 = 22 cycles (v3.1) and 15 + 8 = 23 cycles (v3.2) — exactly the
  Read 1 lengths the paper uses 🟢.
- **rYshapeP7LTnn** = `illumina.P7` + 6-nt index + `illumina.TRUSEQ_READ2` (64 nt); the
  index is written as the **reverse complement** of what Index 1 reads: LT06 `ATTGGC`
  reads `GCCAAT` (= TRSI-6 of 2013), LT12 → `CTTGTA`, LT05 → `ACAGTG`, LT19 → `GTGAAA`,
  LT02 → `CGATGT`, LT04 → `TGACCA`, LT07 → `CAGATC`, LT16 → `CCGTCC`. Its 3' end is
  `GCTCTTCCGATC` + `T`.
- **rYshapeP5** = `illumina.INDEX2_PRIMER_RC[1:21]` (20 nt); its first 12 nt are
  `illumina.STEM` and pair with the `GCTCTTCCGATC` before the LT's 3' T — the same 12-bp
  stem as the 2013 adapter. The other 8 nt (`GTCGTGTA`) are the start of a TruSeq Read 1
  arm, unpaired: a **truncated Y adapter** with no P5. **No 5' phosphate is listed** on
  rYshapeP5, so only the LT strand is ligated (its 3' T to the insert's 5'-phosphate); P5
  enters later from P5-gMac_hybrid. 🟡
- **TPC2** is the 2013 oligo, unchanged.

### B.4 Step by step

🟢 from Methods unless marked.

1. **Plate**: 384-well PCR plates, 1 µL lysis buffer per well (0.1111 µM of that well's RT
   primer, 0.12 mM dNTPs, 0.3 % NP-40, 1 U/µL RNasin Plus, ERCC spike-ins). Sort one cell
   per well (SH800 or MoFlo Astrios EQ; PI / Calcein-AM / Hoechst gating), spin
   10,000 g 1 min, mix, freeze at −80 °C. No air bubble is mentioned for the 1 µL wells
   (the 2013 protocol used one with 0.4 µL); the paper does not say why. 🟡
2. **Denature / prime**: 70 °C 90 s, 35 °C 15 s.
3. **RT (RT25)**: + 1 µL (2× Thermopol buffer, 5 U/µL SuperScript III, 0.55 U/µL RNasin
   Plus); 35 °C 5 min, **50 °C 50 min**, 70 °C 15 min. (RT100 = 20 U/µL.) Every first strand
   now carries `M·AC·<cbc>·<umi>·T24`.
4. **Pool**: plates inverted onto a collector and spun (3010 g, 3 min); ~650–700 µL per
   384-well plate. Purify with Zymo DNA Clean & Concentrator-5 (3 columns per plate for
   v3.1; 8 columns per four plates for v3.2), 20 µL eluate per column. **No exonuclease I.**
5. **Poly(A) tailing + RNase H**: + 25 µL (1× Thermopol, 2.4 mM dATP, 0.0384 U/µL RNase H,
   26.88 U/µL TdT); 37 °C **75 s**, 65 °C 10 min. (25 µL of 1× Thermopol mix on the 20 µL
   eluate gives ≈ 0.56× Thermopol in 45 µL, i.e. the "T55" buffer of Table S5, which is
   0.55× Thermopol/T100 🟡 computed.)
6. **Second strand**: split ~11 µL into four tubes, + 46.16 µL PCR I premix (MightyAmp v2,
   0.06932 µM tagging primer, MightyAmp); 98 °C 130 s, 40 °C 1 min, then **"Increment"**:
   up to 68 °C at 0.2 °C per second, 68 °C 5 min.
7. **Suppression PCR**: + 50.232 µL PCR II premix (1.8952 µM gM primer); **11 cycles**
   (v3.1) or **9 cycles** (v3.2) of 98 °C 10 s, 65 °C 15 s, 68 °C 5 min; 68 °C 5 min.
   Pool, MinElute (with sodium acetate + PB), then AMPure XP. Average cDNA ~1400 bp.
8. **Truncated adapter**: 5 µL each of 100 µM rYshapeP5 and one rYshapeP7LTxx; 90 °C 90 s,
   cool to 10 °C at 0.5 °C per 30 s; dilute to 10 µM.
9. **Library**: 5–10 ng cDNA → Covaris LE220 (duty 15 %, 450 W, 200 cycles/burst, 80 s) →
   DCC-5 into 10 µL → KAPA End repair & A-tailing 37 °C 60 min, 65 °C 30 min → + 2 µL
   1.5 µM truncated adapter + KAPA ligase, 20 °C 15 min → AMPure XP (18 µL to 22 µL) →
   KAPA HiFi with TPC2 and P5-gMac_hybrid (1.75 µL of 10 µM per the text, ambiguous for
   P5-gMac): 98 °C 45 s, **8 ×** (98 °C 15 s, 60 °C 30 s, 72 °C 1 min), 72 °C 5 min →
   AMPure XP (40 µL to 50 µL).
10. **Why only 3' ends amplify** 🟡: after ligation each strand carries the LT arm
    (P7 · i7' · TruSeq Read 2) on its 5' end. TPC2 primes on the copy of that arm; the
    other primer, P5-gMac_hybrid, needs `M·AC` at the far end, which only the RT-primer
    end of a cDNA has. Internal fragments and tagging-primer-end fragments get TPC2 alone
    (linear at best).

Alternative library (mentioned, not the main one) 🟢: Nextera XT on 0.75 ng cDNA,
amplified with P5-gMac_hybrid and a Nextera XT P7 primer. Ligation was preferred because
familiar genes were easier to see (MOESM1).

### B.5 Final library — 🟡 (assembled from the oligos; agrees with upstream)

```
5'- P5[:28] · TT · gM_primer (G·M) · AC · <cbc 14|15> · <umi 8> · T24 · <cDNA, 3' end> · TruSeq Read 2' · i7' · P7' -3'
```

with `TruSeq Read 2' · i7' · P7'` = `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` · index ·
`ATCTCGTATGCCGTCTTCTGCTTG` (the reverse complement of the LT oligo's 3' part). Nextera
variant: same 5' half, then `<cDNA> · ME' · s7' · i7' · P7'`. 🟡

### B.6 Read layout

🟢 NextSeq 500 (mainly; High Output v2 75-cycle kit) or HiSeq 2500.

| Read | Primer | Length | Content |
|---|---|---|---|
| Read 1 | **Read1DropQuartz** (custom) | 22 (v3.1) / 23 (v3.2) | cell barcode + UMI 🟡 (from the oligo layout) |
| Index 1 | standard (`illumina.INDEX1_PRIMER` 🟡) | 6 | pool (sample) barcode from rYshapeP7LT |
| Read 2 | standard TruSeq Read 2 🟡 | 64–118 (v3.1), 63 (v3.2); mostly trimmed to 62 | cDNA, from the fragment end towards the poly(A) |

Pool barcodes allow several 384/1536-cell pools per run; data processing follows Drop-seq
tools, with a custom barcode corrector for up to two substitutions/indels. 🟢

---

## C. Across the family — 🟡 (computed)

| | Quartz-Seq | Quartz-Seq2 |
|---|---|---|
| RT primer | M · T7 promoter · `GGCG` · T₂₄ (70 nt) | M · `AC` · cbc · N₈ · T₂₄ (73/74 nt) |
| Cell identity | the tube; TruSeq i7 per cell | 14/15-nt cbc in RT primer; i7 = pool |
| Leftover primer removal | exonuclease I in the tube | pooling + column purification |
| TdT | 37 °C 50 s, PCR buffer | 37 °C 75 s, Thermopol (T55) |
| Second strand | 98 / 40 / 68 °C | 98 / 40 °C, ramp 0.2 °C/s to 68 °C |
| PCR primer | (NH2)-gM, 21 cycles | gM, 11 / 9 cycles |
| Library | full-length, TruSeq Y adapter, TPC1 + TPC2 | 3' end, truncated Y adapter, P5-gMac + TPC2 |
| Read 1 | TruSeq Read 1 into cDNA | Read1DropQuartz into cbc + UMI |

The 2018 P7 arm is the 2013 one turned around: for every index number both papers use
(2, 4, 5, 6, 7, 12), **rYshapeP7LTnn = reverse complement of TRSI-n + a 3' `T`**. In 2013
the P7 strand is the phosphorylated short-overhang strand and P5 comes from TRSU; in 2018
the P7 strand carries the 3' T overhang and is the one that ligates, while the P5 side is
reduced to the 20-nt rYshapeP5.

Shared verbatim: Tagging primer, the 26-nt gM / suppression ("B") primer sequence, TPC2,
the 12-bp `GATCGGAAGAGC` adapter stem, and TruSeq LT index numbering (TRSI-6 and
rYshapeP7LT06 give the same Index 1 read `GCCAAT`).

---

## D. Upstream (scg_lib_structs) checked against the papers

| Upstream claim | Papers | Verdict |
|---|---|---|
| RT primer, Tagging primer, Suppression primer (NH2), TRSU, TPC1, TPC2 sequences | Table 1 | **agree** |
| TRSI written generically with `NNNNNN` | Table 1 lists six specific TRSI (2, 4, 5, 6, 7, 12) | agree in form; upstream gives no list |
| TRSI has no 5' modification | Table 1: **(5'-phosphate)** on every TRSI | **disagree** — upstream omits the phosphate |
| Quartz-Seq steps go RT → RNase H/poly(A) tailing → tagging | Methods: **exonuclease I** digestion between RT and tailing | **upstream omits a step**; its drawing of tailed leftover RT primer is still correct, since the paper says Exo I is incomplete |
| Fragments from both cDNA ends "cannot be fully ligated" because of the 5' NH2, so only internal fragments are sequenceable | Paper: 7.68 % of reads are WTA adaptor; trimmed bioinformatically. One strand of an NH2-end fragment can still receive both adapters (TRSI ligates to the 3' end of the complementary strand) | **doubtful** 🟡 — end fragments are reduced, not excluded; the NH2 rationale is not in the paper |
| "T7 promoter only for Quartz-chip" | Same RT primer for Quartz-Seq and Quartz-Chip; T7 used for cRNA labelling | agree (the T7 is present in Quartz-Seq libraries, just unused) |
| Quartz-Seq2: eMDRT0001–1536, 15-bp cbc, 8-bp UMI | Table S4: also a **v3.1 set of 384 × 14-nt** (MDRT001–384) | agree; upstream covers only v3.2 |
| gM primer = suppression primer | sequences identical; Table S4 shows no NH2 | agree |
| rYshapeP5, rYshapeP7LT, TPC2, Read1DropQuartz | Table S4 | **agree** (sequences) |
| P5-gMac_hybrid `…CGCGATAC` | Table S4 `…cgcgaT*A*c` (two phosphorothioates, mixed case) | sequence agrees; **upstream drops the phosphorothioates** |
| Leftover RT primer gets tailed in Quartz-Seq2 | Pool is column-purified first; byproduct risk discussed in Fig. S2 | consistent, but upstream omits the pooling / column step |
| Read 1 23 cycles, Index 6 | v3.2: 23 / 6 / 63; v3.1: 22 / 6 / 64–118 | agree for v3.2 |
| Index 1 primer and TruSeq Read 2 primer | not in Table S4 | upstream fills in standard primers 🟡 |

---

## E. Open questions

- 🔴 Units: Quartz-Seq gives the RT primer at "41.67 pmol/l" and the tagging primer at
  "70 pmol/l"; Quartz-Seq2's tagging primer is 0.06932 µM (≈ 69 nM). The 2013 "pmol/l"
  may be a typo for nmol/l. 🟡
- 🔴 Why the Quartz-Seq suppression primer carries a 5' amino group (Table 1) when the
  Quartz-Seq2 gM_primer, same sequence, does not.
- 🔴 Why P5-gMac_hybrid uses only P5[:28] followed by `T`; and whether its mixed case means
  anything.
- 🔴 The concentration of P5-gMac_hybrid in the library PCR (the sentence gives 1.75 µL of
  10 µM only for TPC2 unambiguously).
- 🟡 Read 2 primer for the Nextera variant of Quartz-Seq2 is not stated.
- 🔴 Function of the `GGCG` between the T7 promoter and T₂₄ beyond the T7 +1 `GGG` (not
  explained).

## F. How this note was made (tool evaluation)

`get_sources.py` fetched both papers (PMC full text + XML) and both supplementary zips from
Europe PMC, and `doctext.py` made text twins for the PDF/DOCX/XLSX/HTML files.
`scrape_primers.py --max-hits 80` ranked Quartz-Seq2's Table S4 first and found Quartz-Seq's
Table 1 in the PMC text, recognising P5/P7/TruSeq parts. Misses, worked around by reading
the text: the lower-case **P5-gMac_hybrid** (`aatgatacgg…T*A*c`) was not reported as a hit,
only in other hits' context; the legacy `.XLS` supplements (S4, S7, S8) and the image-only
`S3.PDF` got no usable text. The order of reactions came from reading the Methods.

Verification of the 🟢 sequences with `scrape_primers.py --find`: most were located, but
`--find` reported no location (or only the upstream page) for oligos written with a `*`
inside the sequence (TRSU, all TRSI-n — it does not bridge the phosphorothioate mark), with
literal `N` (the MDRT / eMDRT rows), or in lower case (P5-gMac_hybrid in Table S4). Each of
these was confirmed by an exact `grep` of the written form in the `.txt` source instead
(Table 1 lines 1386–1428 of the 2013 PMC text; Table S4 rows 5, 388, 393, 1928, 1943).
