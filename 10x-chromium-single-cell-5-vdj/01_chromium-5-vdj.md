# 10x Chromium Single Cell 5' V(D)J — barcoded-TSO capture and nested TCR/BCR enrichment

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`). Claims taken only from
> the upstream scg_lib_structs page are 🟡 (secondary source), even where they quote
> sequences.

**Chromium Single Cell 5' Immune Profiling (V(D)J)** — a commercial kit from 10x
Genomics. There is **no defining publication** (the catalogue row says "vendor/kit
protocol"); the primary source is the vendor user guide:

- 10x Genomics, *Chromium Next GEM Single Cell 5' Reagent Kits v2 (Dual Index) User
  Guide*, CG000331, Rev F (© 2024; Rev E cross-checked). Kits: Single Cell 5' Kit v2
  PN-1000263/1000265, Gel Bead Kit v2 PN-1000264/1000267, Library Construction Kit
  PN-1000190, Human/Mouse TCR/BCR Amplification Kits PN-1000252…1000255, Chip K,
  Dual Index Kit TT Set A PN-1000215. 🟢
- Upstream drawing: Teichlab scg_lib_structs, "10x Chromium 5' Immune Profiling Feature
  Barcoding" (based on the v2 Feature Barcoding user guide, Rev A, which upstream names
  but which was **not** fetched here).

Sources read (fetched into `$CHEM_DATA/sources/10x-chromium-single-cell-5-vdj/`, never
committed):

| File | What | Used for |
|---|---|---|
| `upstream_10xChromium5vdjfb.html.txt` | scg_lib_structs method page (fetched by `tools/get_sources.py`) | second account of oligos, steps, all three libraries (GEX, Feature, V(D)J), read layout; the **only** source for Feature Barcode and v3 |
| `CG000331_5v2_UserGuide_RevF.pdf(.txt)` | 10x user guide CG000331 Rev F (fetched by hand, see below) | **all oligos** (Appendix "Oligonucleotide Sequences"), every reaction step and condition, sequencing |
| `CG000331_5v2_UserGuide_RevE.pdf(.txt)` | the same guide, Rev E | cross-check: the set of sequences ≥12 nt is identical to Rev F |

The two PDFs were not found by `get_sources.py` (no DOI); they were downloaded by hand
from public mirrors at the NCI CRTP site
(`crtp.ccr.cancer.gov/wp-content/uploads/fsgc-protocols/Chromium_Next_GEM_Single_Cell_5_v2_UserGuide_RevF.pdf`,
`crtp.ccr.cancer.gov/wp-content/uploads/2023/10/CG000331_ChromiumNextGEMSingleCell5-v2_UserGuide_RevE.pdf`).

Not fetched (get by hand into the same directory):

| Document | Why it matters | URL |
|---|---|---|
| CG000330, 5' v2 Dual Index **with Feature Barcode** (cell surface protein + immune receptor), the guide upstream drew from | Feature Barcode oligo, "Feature cDNA Primers 4", Dual Index Kit TN Set A, feature library steps — all 🟡 here | https://support.10xgenomics.com/permalink/user-guide-chromium-single-cell-5-reagent-kits-user-guide-v2-chemistry-dual-index-with-feature-barcoding-technology-for-cell-surface-protein-and-immune-receptor-mapping |
| GEM-X Universal 5' v3 user guide | the 12-nt UMI and new barcode whitelist claimed by upstream | https://www.10xgenomics.com/support/universal-five-prime-gene-expression |
| Dual Index Kit TT Set A index table | the 10-nt i5/i7 sequences and their written orientation | 10x support site, "Dual Index Kit TT Set A" sample-index CSV |

---

## 1. What it is

A droplet (GEM) method in which the **cell barcode sits on the template-switching
oligo**, not on the oligo-dT. Each Gel Bead releases a barcoded TSO
(partial TruSeq Read 1 · 16-nt 10x barcode · 10-nt UMI · 13-nt TSO ending `rGrGrG`); a
free, unbarcoded poly(dT)VN primer starts reverse transcription. The barcode is therefore
copied onto the **5' end of the transcript** 🟢 (guide Appendix, GEM-RT product). That
is the end nearest the V(D)J rearrangement, which is why a 5' kit rather than the 3' kit
is used for TCR/BCR 🟡 (inference, not stated in this form in the guide). The mechanism is
[template switching](../ref/concepts/template-switching.md); upstream calls it
"conceptually very similar to STRT-seq". 🟡

From one pooled, amplified full-length cDNA three libraries can be made: **5' gene
expression** (GEX), **TCR** and/or **BCR** (V(D)J), and — with the Feature Barcode kit —
**cell-surface protein** (FB). This note concentrates on V(D)J; GEX and FB differences are
in §4–5.

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Barcode + UMI on the TSO**, delivered by a dissolving Gel Bead | 10x 3' puts them on the oligo-dT; STRT-seq / SMART-seq put a handle (no barcode) on the TSO |
| 2 | **Nested, constant-region-specific PCR** (outer Mix 1, inner Mix 2) from whole-transcriptome cDNA, with one universal forward primer in the partial-Read-1 handle | targeted 5'-RACE-style repertoire protocols |
| 3 | **Fragment + A-tail + ligate** a TruSeq-Read-2 half adaptor, so only fragments keeping the barcoded 5' end become sequenceable; Read 2 tiles across V(D)J from random break points | the 10x 3' v3 / 5' GEX library construction |

## 2. Oligos

🟢 Verbatim from CG000331 Rev F, Appendix "Oligonucleotide Sequences" (pp. 73–75;
`.txt` lines ~6790–7190). `rG` = ribo-G; `N16` = 10x barcode, `N10` = UMI, `N10` in index
primers = sample index; dashes are the guide's segment separators.

```
Gel Bead Primer (barcoded TSO)
  5'-CTACACGACGCTCTTCCGATCT-NNNNNNNNNNNNNNNN-NNNNNNNNNN-TTTCTTATATrGrGrG-3'

Poly-dT RT Primer  PN-2000007
  5'-AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3'

cDNA Primers  PN-2000089
  Forward (Partial Read 1)   5'-CTACACGACGCTCTTCCGATCT-3'
  Reverse (Non-poly(dT))     5'-AAGCAGTGGTATCAACGCAGAG-3'

V(D)J Amplification forward primer, identical in all eight Mix 1 / Mix 2 tubes
  5'-GATCTACACTCTTTCCCTACACGACGC-3'

Adaptor Oligos  PN-2000094
  5'-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3'
  5'-GCTCTTCCGATCT-3'   (the partner strand, written 5'->3'; upstream draws it 3'->5')

Dual Index Kit TT Set A  PN-1000215
  Forward  5'-AATGATACGGCGACCACCGAGATCTACAC-N10-ACACTCTTTCCCTACACGACGCTC-3'
  Reverse  5'-CAAGCAGAAGACGGCATACGAGAT-N10-GTGACTGGAGTTCAGACGTGT-3'
```

**Reverse outer primers (V(D)J Amplification 1)** 🟢

```
Human T Cell Mix 1 v2  PN-2000242
  TGAAGGCGTTTGCACATGCA
  TCAGGCAGTATCTGGAGTCATTGAG
Human B Cell Mix 1 v2  PN-2000254
  CAGGGCACAGTCACATCCT
  TGCTGGACCACGCATTTGTA
  GGTTTTGTTGTCGACCCAGTCT
  TTGTCCACCTTGGTGTTGCT
  CATGACGTCCTTGGAAGGCA
  TGTGGGACTTCCACTG
  TTCTCGTAGTCTGCTTTGCTCAG
Mouse T Cell Mix 1 v2  PN-2000256
  CTGGTTGCTCCAGGCAATGG
  TGTAGGCCTGAGGGTCCGT
Mouse B Cell Mix 1 v2  PN-2000258
  TCAGCACGGGACAAACTCTTCT
  GCAGGAGACAGACTCTTCTCCA
  AACTGGCTGCTCATGGTGT
  TGGTGCAAGTGTGGTTGAGGT
  TGGTCACTTGGCTGGTGGTG
  CACTTGGCAGGTGAACTGTTTTCT
  AACCTTCAAGGATGCTCTTGGGA
  GGACAGGGATCCAGAGTTCCA
  AGGTGACGGTCTGACTTGGC
  GCTGGACAGGGCTCCATAGTT
  GGCACCTTGTCCAATCATGTTCC
  ATGTCGTTCATACTCGTCCTTGGT
```

**Reverse inner primers (V(D)J Amplification 2)** 🟢

```
Human T Cell Mix 2 v2  PN-2000246
  AGTCTCTCAGCTGGTACACG
  TCTGATGGCTCAAACACAGC
Human B Cell Mix 2 v2  PN-2000255
  GGGAAGTTTCTGGCGGTCA
  GGTGGTACCCAGTTATCAAGCAT
  GTGTCCCAGGTCACCATCAC
  TCCTGAGGACTGTAGGACAGC
  CACGCTGCTCGTATCCGA
  TAGCTGCTGGCCGC
  GCGTTATCCACCTTCCACTGT
Mouse T Cell Mix 2 v2  PN-2000257
  AGTCAAAGTCGGTGAACAGGCA
  GGCCAAGCACACGAGGGTA
Mouse B Cell Mix 2 v2  PN-2000259
  TACACACCAGTGTGGCCTT
  CAGGCCACTGTCACACCACT
  CAGGTCACATTCATCGTGCCG
  GAGGCCAGCACAGTGACCT
  GCAGGGAAGTTCACAGTGCT
  CTGTTTGAGATCAGTTTGCCATCCT
  TGCGAGGTGGCTAGGTACTTG
  CCCTTGACCAGGCATCC
  AGGTCACGGAGGAACCAGTTG
  GGCATCCCAGTGTCACCGA
  AGAAGATCCACTTCACCTTGAAC
  GAAGCACACGACTGAGGCAC
```

All written 5'→3' in the guide. The assignment of primers to mixes was read from the
two-column PDF layout (`pdftotext -layout`), because the plain `.txt` twin interleaves
the columns of the Mix 2 table. 🟡 Which gene each primer targets is **not stated** in the
guide 🔴; upstream says only that they "anneal to the constant region". The counts — 2
per T-cell mix (consistent with TRA + TRB constant regions), 7 per human B mix, 12 per
mouse B mix — are 🟡 suggestive of one primer per isotype / light-chain constant region,
not verified by alignment.

**Feature Barcode variant — upstream only** 🟡 (the FB user guide was not fetched):

```
Barcoded oligo on antibody (TotalSeq-C-type)
  5'-CGGAGATGTGTATAAGAGACAGNNNNNNNNNN[15-bp FB]NNNNNNNNNCCCATATAAGAAA-3'
Feature cDNA Primers 4  PN-2000277 (three primers)
  Forward            5'-CTACACGACGCTCTTCCGATCT-3'
  Reverse (cDNA)     5'-AAGCAGTGGTATCAACGCAGAG-3'
  Reverse (FB)       5'-CTCGTGGGCTCGGAGATGTG-3'
Dual Index Kit TN Set A  PN-3000510
  Forward  5'-AATGATACGGCGACCACCGAGATCTACAC[10-bp i5]ACACTCTTTCCCTACACGACGCTC-3'
  Reverse  5'-CAAGCAGAAGACGGCATACGAGAT[10-bp i7]GTCTCGTGGGCTCGG-3'
```

### How the oligos interlock — 🟡 (computed)

- **Gel Bead Primer** (61 nt with v2's 10-nt UMI; 63 nt with v3's 12-nt UMI per upstream)
  = `illumina.TRUSEQ_READ1[11:]` (the last 22 nt of TruSeq Read 1 — "partial Read 1") +
  N16 + N10 + the 13-nt TSO `TTTCTTATATGGG` (last three bases ribo).
- **Poly-dT RT Primer** (57 nt) = **`rt.SMART_HANDLE`** (23 nt, `…CGCAGAGT`) + `AC` +
  T₃₀ + `VN`. So the 3' end of the cDNA carries the classic SMART handle, and the 5' end
  carries partial Read 1 — the barcoded end and the handle end use different primers.
- **cDNA reverse primer** = `rt.SMART_HANDLE[:22]`; **cDNA forward primer** = the 22-nt
  partial Read 1 itself. Amplified cDNA top strand therefore ends
  `…GTACTCTGCGTTGATACCACTGCTT`, = revcomp of the first 25 nt of the RT primer. The
  guide's GEM-RT product bottom strand (3'→5' `GATGTGCTGCGAGAAGGCTAGA…AAAGAATATACCC`) is
  exactly the plain complement of partial Read 1 and of the TSO.
- **V(D)J forward primer** (27 nt) = `illumina.NEBNEXT_UNIVERSAL_PRIMER[20:47]`, i.e.
  the last 9 nt of **P5** (`GATCTACAC`, whose `ACAC` is also the first 4 of Read 1)
  running into `TRUSEQ_READ1[:22]`. Its 3' 11 nt (`CTACACGACGC`) are the 5' 11 nt of the
  partial Read 1 on the cDNA, so it anneals by 11 nt in the first cycle and adds a 16-nt
  5' extension; the V(D)J product begins
  `GATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` (38 nt) = `GATCT` + the **complete** 33-nt
  TruSeq Read 1. Both PCRs (outer and inner) use the same forward primer, so nesting is
  only on the constant-region side.
- **Adaptor Oligos**: top strand (33 nt) = `illumina.INDEX1_PRIMER` =
  `illumina.CANONICAL_ARM_READ2`; bottom strand read 5'→3' is `GCTCTTCCGATCT` (13 nt) =
  `illumina.STEM_COMPLEMENT` + a 3' **T overhang**. The 12-nt stem pairs with the top
  strand's first 12 nt (`illumina.STEM`); the T pairs with the A-tail of the fragment. A
  half adaptor, not a Y adaptor: it carries no Read-1 side at all.
- **TT index primers**: forward = `illumina.P5` + i5 + `TRUSEQ_READ1[:24]` (63 nt); its
  Read-1 part matches the V(D)J product from position 5, after the `GATCT` that the
  V(D)J forward primer left. Reverse = `illumina.P7` + i7 + `TRUSEQ_READ2[:21]` (55 nt);
  revcomp of its 21-nt 3' part is the **last** 21 nt of the adaptor top strand. Library
  insert-side junction `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` = `A` (tail) + adaptor top
  = revcomp(`TRUSEQ_READ2`).
- **Feature Barcode** (upstream oligos): the antibody oligo's 3' 13 nt
  `CCCATATAAGAAA` = **revcomp of the 13-nt TSO** — it hybridises to the Gel Bead Primer
  directly and serves as template for its extension, no RT of RNA involved. Its 5' 22 nt
  `CGGAGATGTGTATAAGAGACAG` = `nextera.ADAPTOR_S7[12:]` (last 3 nt of s7 + the full
  `nextera.ME`). The FB cDNA reverse primer = `nextera.ADAPTOR_S7[2:22]`; the TN-set
  reverse index primer's 3' part = `nextera.S7` (15 nt). So the protein library is read
  with the **Nextera** Read 2 / i7 primers, the V(D)J and GEX libraries with **TruSeq**
  ones.

## 3. Step by step (V(D)J; CG000331 Rev F)

Volumes and conditions 🟢 from the guide.

1. **GEM generation** (Chip K, Chromium Controller): cells + Master Mix (RT Reagent B
   18.8 µl, Poly-dT RT Primer 7.3 µl, Reducing Agent B 1.9 µl, RT Enzyme C 8.3 µl per
   reaction) + Gel Beads + Partitioning Oil. The Gel Bead dissolves and releases the
   barcoded TSO; the cell lyses.
2. **GEM-RT** 53 °C 45 min, 85 °C 5 min (125 µl). Poly(dT)VN primes on poly(A); the
   RT adds non-templated C's, the TSO's `rGrGrG` pairs and the RT switches onto the Gel
   Bead Primer, copying barcode, UMI and partial Read 1 onto the cDNA 3' end. 🟢 product;
   mechanism 🟡 per [template switching](../ref/concepts/template-switching.md) and
   [reverse transcription](../ref/concepts/reverse-transcription.md). First strand,
   5'→3': `SMART handle · AC · T30 · NV · <cDNA, antisense> · CCC · TSO' · UMI' · BC' · partial R1'`.
3. **Break GEMs** (Recovery Agent), **Dynabeads MyOne SILANE** cleanup, elute 35 µl.
4. **cDNA amplification**, 65 µl mix (Amp Mix + cDNA Primers) + 35 µl: 98 °C 45 s;
   98 °C 20 s / 63 °C 30 s / 72 °C 1 min; 72 °C 1 min. Total cycles 11–16 by cell
   number and RNA content (table in guide). **SPRIselect 0.6×**, elute 45 µl.
5. **V(D)J Amplification 1** (outer): 2 µl cDNA + 50 µl Amp Mix + 48 µl T- or B-Cell
   Mix 1; 98 °C 45 s; 98 °C 20 s / 62 °C 30 s / 72 °C 1 min — **12 cycles for T, 8 for B**;
   72 °C 1 min. T and B are separate reactions from the same cDNA.
6. **Double-sided SPRIselect** 0.5× then 0.8× (keep supernatant of the first, beads of
   the second), elute 35 µl.
7. **V(D)J Amplification 2** (inner): 35 µl + 50 µl Amp Mix + 15 µl Mix 2; same profile,
   **10 cycles for T, 8 for B**. Double-sided SPRIselect 0.5×/0.8×, elute 45 µl.
   Product (guide): `GATCT · TruSeq R1 · BC · UMI · TSO · 5'UTR-V(D)J-C · inner primer`.
8. **QC**: Bioanalyzer HS, region ~200–9000 bp; carry ≤50 ng in 20 µl.
9. **Fragmentation, end repair, A-tailing** (one tube, 50 µl): 32 °C **2 min**, then
   65 °C 30 min. (GEX uses 32 °C **5 min** — shorter fragmentation for V(D)J.)
10. **Adaptor ligation**: + 50 µl (Ligation Buffer 20, DNA Ligase 10, Adaptor Oligos 20),
    20 °C 15 min. SPRIselect 0.8×, elute 30 µl.
11. **Sample index PCR**: 30 µl + 50 µl Amp Mix + 20 µl Dual Index TT Set A; 98 °C 45 s;
    8 × (98 °C 20 s, 54 °C 30 s, 72 °C 20 s); 72 °C 1 min. SPRIselect 0.8×, elute 35 µl.

Which fragments become library 🟡 (inferred from primer positions): after step 9 the
V(D)J amplicon is cut randomly; every piece is A-tailed and gets the half adaptor on both
ends. Only the **5'-most piece** carries the Read-1 handle the forward index primer needs,
so only it acquires P5. Internal pieces have the adaptor on both ends and can at most be
copied by the reverse primer alone into P7…P7 molecules that cannot cluster (upstream
calls them "not amplifiable"). The constant-region end (inner primer) is likewise lost
unless no break occurred. Read 2 therefore starts at a random point inside the
amplicon and reads back toward the 5' end — collectively tiling C → J → D → V.

## 4. Final library structures — 🟡 (assembled from the guide's appendix)

**V(D)J (TCR or BCR)** — the guide's appendix writes exactly this:

```
5'- P5 · i5(10) · TruSeq Read 1 (33) · BC(16) · UMI(10) · TSO TTTCTTATATGGG (13) · <5'UTR-V-D-J-C fragment> · A · adaptor top (A + adaptor top = revcomp of TruSeq Read 2) · i7(10) · P7' -3'
```

Fixed length outside the insert: 178 bp (v2), 180 bp with a 12-nt UMI (v3), not counting
the A-tail base (179 / 181 with it). 🟡 (computed: 29+10+33+16+10+13 on the P5 side,
33+10+24 on the P7 side)

**5' GEX** — identical layout; the insert is any fragmented cDNA rather than a TCR/BCR
amplicon (guide steps 5.1–5.5, index PCR cycles by cDNA input). 🟢 structure.

**Feature Barcode** (upstream only) 🟡:

```
5'- P5 · i5(10) · TruSeq Read 1 · BC(16) · UMI(10) · TSO(13) · N9 · FB(15) · N10 · ME' · s7' · i7(10) · P7' -3'
```

## 5. Sequencing

🟢 Paired end, dual index: **Read 1 26 cycles** (16 BC + 10 UMI), **i7 10**, **i5 10**,
**Read 2 90 cycles**. Depth ≥5,000 read pairs per cell for V(D)J, ≥20,000 for GEX.
Primers 🟡 (computed): Read 1 = `illumina.TRUSEQ_READ1`; i7 read = `illumina.INDEX1_PRIMER`
(= adaptor top strand); i5 read primed by `illumina.INDEX2_PRIMER_RC` (= revcomp of `TRUSEQ_READ1`, upstream's "Truseq i5 index sequencing primer", on the regenerated top strand); Read 2 =
`illumina.TRUSEQ_READ2`. Feature library: i7 read with `nextera.INDEX1_PRIMER`, Read 2 with
`nextera.READ2_PRIMER` (upstream).

Read 1 at 26 cycles ends exactly at the end of the UMI; the 13-nt TSO and the 5' UTR are
never read in Read 1. Read 2 (90 nt) covers V(D)J only through overlapping fragments, so
full-length contigs are assembled computationally (Cell Ranger `vdj`). 🟡

## 6. Upstream scg_lib_structs vs the user guide

| Point | Upstream page | CG000331 Rev F | Verdict |
|---|---|---|---|
| Gel Bead Primer, Poly-dT RT primer, cDNA primers, V(D)J forward primer, adaptor, TT index primers | sequences given | same sequences | **agree** (all found by `--find` in both files) |
| Outer / inner reverse primers | "check the user guide" | 2/7/2/12 primers per mix, sequences listed | upstream incomplete; guide fills it |
| Mix part numbers | 2000242/54/56/58, 2000246/55/57/59 | same | agree |
| cDNA amplification primers | "Feature cDNA Primers 4", PN-2000277, three primers (adds FB reverse) | "cDNA Primers", PN-2000089, two primers | different kits (FB vs non-FB guide), not a conflict |
| Dual Index Kit TT Set A part number | PN-3000431 | PN-1000215 (kit) | **disagree** on PN (probably plate vs kit number 🟡); primer sequences identical |
| Adaptor part number | PN-220026 | PN-2000094 (DNA Ligase is 220110/220131) | **disagree** on PN; sequence identical |
| Index kit for the protein library | TN Set A, PN-3000510 | not in this guide | upstream only 🟡 |
| UMI length | 10 nt v1/v2, 12 nt v3; Read 1 26 / 28 cycles | 10 nt, Read 1 26 cycles (v2 only) | agree for v2; v3 unverified |
| Sequencing | R2 90, i7 10, i5 10 | same | agree |
| Fragmentation | "fragment and A-tailing" | 32 °C 2 min (V(D)J) vs 5 min (GEX), then 65 °C 30 min | guide adds detail |

## 7. Open questions

- 🔴 Target gene of each outer/inner primer, and how far inner primers sit upstream
  (5') of outer ones in each constant region — not published; would need alignment to
  TRAC/TRBC/IGH*/IGK/IGL references.
- 🔴 5'-phosphorylation of the adaptor top strand (needed for ligation) — not shown in
  the guide's sequences.
- 🔴 The 10-nt i5/i7 sequences of Dual Index Kit TT Set A, and whether they are written
  as read or as reverse complement in the primers.
- 🟡 Feature Barcode oligo, primers and TN index kit — upstream only; fetch CG000330.
- 🟡 v3 (GEM-X) 12-nt UMI and 28-cycle Read 1 — upstream only.
- 🟡 What the "13 nt TSO" of the guide is chemically: the appendix writes the 3' three
  bases as ribo-G (`rGrGrG`) and the rest as DNA; whether any LNA is used is not stated.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` found only the upstream page (kit protocol, no DOI, so no paper
to resolve). The vendor guide was located by web search and downloaded by hand from a
public mirror; `tools/doctext.py` converted both revisions. `tools/scrape_primers.py`
found every appendix oligo and recognised `rt.SMART_HANDLE`, `illumina.TRUSEQ_READ1`,
`INDEX1_PRIMER`, `P5`, `P7`, but the `.txt` twin interleaves the two-column table of
V(D)J mixes so that mix membership of the inner primers is ambiguous; `pdftotext -layout`
was needed to assign them. Conditions were read from the guide, not scraped.
