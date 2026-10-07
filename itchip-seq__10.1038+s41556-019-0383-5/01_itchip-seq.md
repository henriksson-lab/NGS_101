# itChIP-seq — simultaneous indexing and tagmentation-based ChIP-seq

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available to us. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note; they are not yet asserted in a self-test, because
> this protocol has no `tools/` module yet.

**itChIP-seq** — Ai S, Xiong H, Li CC, Luo Y, Shi Q, Liu Y, Yu X, Li C, He A. "Profiling
chromatin states using single-cell itChIP-seq." *Nature Cell Biology* 21, 1164–1172
(2019; Technical Report, published 3 Sep 2019). doi:[10.1038/s41556-019-0383-5](https://doi.org/10.1038/s41556-019-0383-5),
PMID 31481796. Data: GEO GSE109762. Code: github.com/Helab-bioinformatics/itChIP.

Other papers:

- **Protocol Exchange** — Ai S et al. "itChIP — simultaneous indexing and
  tagmentation-based ChIP-seq" (2019), doi 10.21203/rs.2.11366/v1 as cited in the paper's
  reference list (the catalogue holds pex-555/v1 with an unconfirmed DOI). This is the
  bench protocol, and the source upstream drew from. **Not fetched.** 🔴
- **Chemistry origin**: Amini S et al. 2014 (CPT-seq, *Nat Genet* 46:1343) — the T5/T7
  barcoded transposon design, cited by the paper for Supplementary Fig. 4a. Covered in
  [sci-ATAC-seq](../sci-atac-seq__10.1126+science.aab1601/01_sci-atac-seq.md).
- Cited by upstream for "semi-suppressive PCR": nanoCAGE / CAGEscan (Plessy et al. 2010,
  *Nat Methods*, doi 10.1038/nmeth.1470).

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/itchip-seq__10.1038+s41556-019-0383-5/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41556-019-0383-5_41556_2019_383_MOESM3_ESM.xlsx` | Supplementary Table 1 "Custom primers" | **every oligo** (Nextera / sci-style variant) |
| `s41556-019-0383-5_41556_2019_383_MOESM1_ESM.pdf` | Supplementary Information (figure legends) | Tn5 amount, SDS, opening and tagmentation conditions; Supp. Fig. 4 = custom Nextera layout, HiSeq 2500 |
| `s41556-019-0383-5_41556_2019_383_MOESM2_ESM.pdf` | Reporting summary | paired-end 150 bp |
| `s41556-019-0383-5_41556_2019_383_MOESM4_ESM.xlsx` | Supplementary Table 2 | Tn5 amount per cell number |
| `s41556-019-0383-5_41556_2019_383_MOESM5_ESM.xlsx` | Supplementary Table 3 | library type per sample ("Custom Nextera") |
| `s41556-019-0383-5_41556_2019_383_MOESM6..9_ESM.xlsx` | comparison, gene lists, source data | not chemistry |
| `upstream_itChIP-seq.html` | scg_lib_structs page | **TruSeq-connector variant** (from Protocol Exchange), step-by-step, read layout |
| `article_s41556-019-0383-5.html` | nature.com article page, fetched by hand with curl | abstract, figure titles, reference list only (paywalled) |

Not obtained 🔴 (fetch by hand):

| What | URL |
|---|---|
| Main text + **Methods** (paywalled) | https://www.nature.com/articles/s41556-019-0383-5 |
| Protocol Exchange bench protocol (incl. its "Table 3" of Tn5 barcodes, TruSeq connector primers) | https://protocolexchange.researchsquare.com/article/pex-555/v1 |

So: **no reaction conditions from the Methods** are available to us. What we have is the
complete oligo table, the optimisation conditions in the Supplementary Figure 1 legend,
and upstream's reading of the Protocol Exchange.

---

## 1. What it is

A ChIP-seq that **tagments first and immunoprecipitates second** 🟢 (abstract; upstream):
fixed cells are permeabilised, chromatin is "opened" with SDS, and **barcoded Tn5**
tagments the chromatin **inside each cell or well** — that is the "simultaneous indexing
and tagmentation". Wells are then pooled, the chromatin is immunoprecipitated together,
crosslinks are reversed, and the already-adapted, already-barcoded DNA goes straight to
PCR. One IP serves "tens of single cells to … thousands" 🟢 (abstract). The abstract
claims ~9,000 unique reads per cell 🟢.

Cell identity = the **combination of a T5 and a T7 Tn5 barcode**, both written into the
transposon oligos and so inserted by Tn5 itself 🟡 (Supp. Table 1; Supp. Fig. 4a legend:
"T5 and T7 barcodes are introduced during Tn5 complex assembly"). A second, sample-level
index comes from the library PCR (i5/i7 in the Nextera variant; i7 only in the TruSeq
variant). There is **no UMI** 🟡 (none appears in any oligo).

| | New thing here | Builds on |
|---|---|---|
| 1 | **Barcoded tagmentation of fixed chromatin, then a pooled ChIP** | ChIPmentation (tagment *after* IP, on beads); Drop-ChIP (MNase + droplet barcoding *before* pooled IP) — upstream makes the Drop-ChIP comparison; [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) |
| 2 | Transposon oligos **identical in layout to the CPT-seq / sci-ATAC-seq "Universal Connector"** oligos, extended from 8×12 to 24×25 barcodes | Amini 2014; [sci-ATAC-seq](../sci-atac-seq__10.1126+science.aab1601/01_sci-atac-seq.md) |
| 3 | A **TruSeq "connector" PCR** that grafts TruSeq Read 1/2 onto the sci-style ends so a standard (non-custom-primer) run works | upstream, from the Protocol Exchange and Fig. 5 ("Truseq library preparation method for low-input and single-cell itChIP", title only 🟢) |

## 2. Oligos

🟢 Verbatim from Supplementary Table 1 (`MOESM3_ESM.xlsx`, sheet "Custom primer"), in
table order. Names as given; the table's third column, where present, gives the vendor /
Amini name (e.g. `T5-1_Universal_Connector_A`, "equally to TnExtMErev_primer").
`5Phos/` = 5' phosphate, written exactly as the table writes it.

```
i5-N501                                            AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC
i5-N502                                            AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC
i5-N503                                            AATGATACGGCGACCACCGAGATCTACACTATCCTCTTCGTCGGCAGCGTC
i5-N504                                            AATGATACGGCGACCACCGAGATCTACACAGAGTAGATCGTCGGCAGCGTC
i5-N505                                            AATGATACGGCGACCACCGAGATCTACACGTAAGGAGTCGTCGGCAGCGTC
i5-N506                                            AATGATACGGCGACCACCGAGATCTACACACTGCATATCGTCGGCAGCGTC
i5-N507                                            AATGATACGGCGACCACCGAGATCTACACAAGGAGTATCGTCGGCAGCGTC
i5-N508                                            AATGATACGGCGACCACCGAGATCTACACCTAAGCCTTCGTCGGCAGCGTC
i7-N701                                            CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG
i7-N702                                            CAAGCAGAAGACGGCATACGAGATCTAGTACGGTCTCGTGGGCTCGG
i7-N703                                            CAAGCAGAAGACGGCATACGAGATTTCTGCCTGTCTCGTGGGCTCGG
i7-N704                                            CAAGCAGAAGACGGCATACGAGATGCTCAGGAGTCTCGTGGGCTCGG
i7-N705                                            CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGG
i7-N706                                            CAAGCAGAAGACGGCATACGAGATCATGCCTAGTCTCGTGGGCTCGG
i7-N707                                            CAAGCAGAAGACGGCATACGAGATGTAGAGAGGTCTCGTGGGCTCGG
i7-N708                                            CAAGCAGAAGACGGCATACGAGATCCTCTCTGGTCTCGTGGGCTCGG
i7-N709                                            CAAGCAGAAGACGGCATACGAGATAGCGTAGCGTCTCGTGGGCTCGG
i7-N710                                            CAAGCAGAAGACGGCATACGAGATCAGCCTCGGTCTCGTGGGCTCGG
i7-N711                                            CAAGCAGAAGACGGCATACGAGATTGCCTCTTGTCTCGTGGGCTCGG
i7-N712                                            CAAGCAGAAGACGGCATACGAGATTCCTCTACGTCTCGTGGGCTCGG
T5-1                                               TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-2                                               TCGTCGGCAGCGTCTCCACGCATAGAGGCGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-3                                               TCGTCGGCAGCGTCTCCACGCCCTATCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-4                                               TCGTCGGCAGCGTCTCCACGCGGCTCTGAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-5                                               TCGTCGGCAGCGTCTCCACGCAGGCGAAGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-6                                               TCGTCGGCAGCGTCTCCACGCTAATCTTAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-7                                               TCGTCGGCAGCGTCTCCACGCCAGGACGTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-8                                               TCGTCGGCAGCGTCTCCACGCGTACTGACGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-9                                               TCGTCGGCAGCGTCTCCACGCTAGGTATGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-10                                              TCGTCGGCAGCGTCTCCACGCAACACCTAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-11                                              TCGTCGGCAGCGTCTCCACGCCTCCGAACGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-12                                              TCGTCGGCAGCGTCTCCACGCCAACGGCAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-13                                              TCGTCGGCAGCGTCTCCACGCCAATGTAGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-14                                              TCGTCGGCAGCGTCTCCACGCGGCTACCCGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-15                                              TCGTCGGCAGCGTCTCCACGCAAAGTCCGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-16                                              TCGTCGGCAGCGTCTCCACGCTTCCGCGGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-17                                              TCGTCGGCAGCGTCTCCACGCAGGCACTTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-18                                              TCGTCGGCAGCGTCTCCACGCCTTCAGTGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-19                                              TCGTCGGCAGCGTCTCCACGCGCCGGTAGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-20                                              TCGTCGGCAGCGTCTCCACGCTTCAATCCGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-21                                              TCGTCGGCAGCGTCTCCACGCCCACACACGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-22                                              TCGTCGGCAGCGTCTCCACGCATATTATCGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-23                                              TCGTCGGCAGCGTCTCCACGCCCGAAGCAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T5-24                                              TCGTCGGCAGCGTCTCCACGCGTATCGGTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
T7-1                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-2                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCTCTCCGGACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-3                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCAATGAGCGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-4                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCGGAATCTCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-5                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCTTCTGAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-6                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCACGAATTCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-7                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCAGCTTCAGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-8                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCGCGCATTACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-9                                               GTCTCGTGGGCTCGGCTGTCCCTGTCCCATAGCCGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-10                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCTTCGCGGACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-11                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCGCGCGAGACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-12                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCCTATCGCTCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-13                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCTGTCCGCGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-14                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCTTAAACTTCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-15                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCACCACAACCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-16                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCGCCTCTGGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-17                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCTCGCCCACCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-18                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCCACTAGGCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-19                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCTCGAAGCCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-20                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCGCATGTACCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-21                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCGTTCGAGTCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-22                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCCCGGGCGCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-23                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCAGATTTAACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-24                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCCACCATTGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
T7-25                                              GTCTCGTGGGCTCGGCTGTCCCTGTCCAATAAGACCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
Mosaic-A_primer                                    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Mosaic-B_primer                                    GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
common annealing primer                            5Phos/CTGTCTCTTATACACATCT
Forward i5/i7 PCR primer (Universal_Connector_A)   TCGTCGGCAGCGTCTCCACGC
Reverse i5/i7 PCR primer (part of Universal_Connector_B) GTCTCGTGGGCTCGGCTGTCC
Forward MEA/B PCR primer (part of Tn5ME-A_FC-121-1030) TCGTCGGCAGCGTCAGAT
Reverse MEA/B PCR primer (part of Tn5ME-B_FC-121-1031) GTCTCGTGGGCTCGGAGA
Read 1 Sequencing Primer (Extended Mosaic End A)   GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
Read 2 Sequencing Primer (Extended Mosaic End B)   CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
Index 1 Sequencing Primer (Reverse Complement of Extended Mosaic End B) CTGTCTCTTATACACATCTGAGGCGGAGACGGTG
Index 2 Sequencing Primer (Grafted P5)             AATGATACGGCGACCACCGAGATCTACAC
AATGATACGGCGACCACCGAGATCTACACTTGCCTAGNNNNNNNNNNNNNNNNNNNNNGAATCTTAGCGATCGAGGACGGC
CAAGCAGAAGACGGCATACGAGATCTTGGGTANNNNNNNNNNNNNNNNNNNNNNNNNNNGCGAATTCCACCGTCTCCGCCTC
```

The last two lines sit in the table under "Primer sequece for amplification spiki-in
library" with no names 🟢; they are listed here as written.

Upstream adds the TruSeq-variant oligos 🟡 (secondary: upstream quotes the Protocol
Exchange, which we could not fetch):

```
Connector Primer F      ACACTCTTTCCCTACACGACGCTCTTCCGATCTTCGTCGGCAGCGTCTCCACGC
Connector Primer R      GACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTCTCGTGGGCTCGGCTGTCCCTGT
Universal P5 primer     AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Indexed P7 primer       CAAGCAGAAGACGGCATACGAGAT<i7>GACTGGAGTTCAGACGTGTGCTCTTCCGATCT
ME bottom (upstream)    /Phos/AGATGTGTATAAGAGACAG
```

### How they interlock — 🟡 (computed with `lib/`)

**Transposon oligos (the cell barcodes).**

- All 24 **T5-N** = `TCGTCGGCAGCGTCTCCACGC` (21 nt, "Universal Connector A") + 8-nt
  T5 barcode + `GCGATCGAGGACGGCAGATGTGTATAAGAGACAG` (34 nt, "extended ME A") — 63 nt.
  The connector starts with **`nextera.S5`** (14 nt) and adds `TCCACGC`; the extended ME
  ends in **`nextera.ME`** (19 nt) behind a 15-nt spacer `GCGATCGAGGACGGC`.
- All 25 **T7-N** = `GTCTCGTGGGCTCGGCTGTCCCTGTCC` (27 nt, "Universal Connector B")
  + 8-nt T7 barcode + `CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG` (34 nt, "extended ME B")
  — 69 nt. The connector starts with **`nextera.S7`** (15 nt).
- The **common annealing primer** `CTGTCTCTTATACACATCT` = **`nextera.ME_RC`**: the
  19-nt non-transferred strand, shared by every T5 and T7 transposon.
- Barcodes within each set differ at ≥ 3 of 8 positions (minimum pairwise Hamming
  distance 3 for T5 and for T7).
- **T5-1…T5-8 and T7-1…T7-12 are, number for number, Amini 2014's
  `P5_i5_1..8_Universal_Connector_A_C15_ME` and `P7_i7_1..12_Universal_Connector_B_D15_ME`**
  (checked against Amini Supplementary Table 4, in the sci-ATAC-seq sources). T5-9…24 and
  T7-13…25 are additions. 24 × 25 = **600 barcode combinations** per Tn5 round.
- Amini's non-transferred strand ("pMENTS", `5Phos/CTGTCTCTTATACACATCT`) is the same oligo
  as itChIP's "common annealing primer".

**Library PCR primers (Nextera / sci-style variant).**

- **i5-N501…N508** = **`illumina.P5`** + 8-nt i5 + `TCGTCGGCAGCGTC` (= `nextera.S5`,
  the first 14 nt of every T5 oligo) — 51 nt. i5 written **as named** (N501 contains
  `TAGATCGC`).
- **i7-N701…N712** = **`illumina.P7`** + 8-nt i7 + `GTCTCGTGGGCTCGG` (= `nextera.S7`,
  the first 15 nt of every T7 oligo) — 47 nt. i7 written as the **reverse complement**
  of its name (N701 contains `TCGCCTTA` = rc `TAAGGCGA`).
- So these are standard Nextera index primers **truncated to end at s5 / s7** — they
  prime on the outer end of the T5/T7 connectors and would equally prime a standard
  Nextera fragment.
- "Forward i5/i7 PCR primer" = the full 21-nt Universal Connector A; "Reverse i5/i7 PCR
  primer" = the first 21 nt of Universal Connector B. Where in the workflow they are used
  is not in our sources 🔴 (presumably qPCR or a pre-amplification).
- **Mosaic-A / Mosaic-B primer** = exactly **`nextera.ADAPTOR_S5` / `nextera.ADAPTOR_S7`**
  (the table says "equally to Tn5ME-A / Tn5ME-B"): the standard, unbarcoded Nextera
  transposon oligos, used for the low-input (non-single-cell) itChIP, whose libraries
  are "Custom Nextera" in Supp. Table 3. "Forward/Reverse MEA/B PCR primer" are their
  first 18 nt.

**Custom sequencing primers.**

- Read 1 primer = the extended ME A (identical to the 3' 34 nt of every T5 oligo); Read 2
  primer = the extended ME B. So **both genomic reads start right after the ME** — no
  barcode in Read 1 or Read 2.
- Index 1 primer = reverse complement of extended ME B (stated in its name; confirmed).
  It reads outward from the insert across the T7 barcode, then the 27-nt connector B,
  then the PCR i7: **8 + 27 + 8 = 43 cycles** to reach the end of i7.
- Index 2 primer = **`illumina.P5`** (29 nt; "grafted P5"). It reads i5, then the 21-nt
  connector A, then the T5 barcode: **8 + 21 + 8 = 37 cycles**.
- This is the CPT-seq / sci-ATAC-seq custom-primer layout.

**Spike-in library primers.**

- Line 1 = `illumina.P5` + `TTGCCTAG` + N₂₁ + `GAATCTTA` + `GCGATCGAGGACGGC`;
  the N₂₁ has exactly the length of connector A, `GAATCTTA` sits in the T5-barcode
  position (1 mismatch from T5-6 `TAATCTTA`), and the tail is the first 15 nt of
  extended ME A.
- Line 2 = `illumina.P7` + `CTTGGGTA` + N₂₇ + `GCGAATTC` + `CACCGTCTCCGCCTC`; N₂₇ =
  length of connector B, `GCGAATTC` is 1 mismatch from T7-6 `ACGAATTC`, tail = first
  15 nt of extended ME B.
- So the spike-in fragments read like real itChIP fragments with a fixed, out-of-set
  barcode pair and random sequence in the constant regions. How the spike-in is made and
  used is not in our sources 🔴.

**TruSeq connector variant (upstream).**

- Connector Primer F = **`illumina.TRUSEQ_READ1`** + full connector A (54 nt).
- Connector Primer R = TruSeq Read 2 **minus its first two bases `GT`** + the first 25 nt
  of connector B (57 nt). `GT` + its first 32 nt = **`illumina.TRUSEQ_READ2`** exactly.
  Upstream argues the `GT` should be there and draws it with `GT` 🟡 — consistent with
  our computation, but the original Protocol Exchange text is unchecked 🔴.
- Universal P5 primer = **`illumina.TRUSEQ_P5_FULL`** (P5 + TruSeq Read 1).
- Indexed P7 primer = `illumina.P7` + i7 + the same `GT`-less TruSeq Read 2 stub;
  upstream again adds `GT`.

### Upstream vs. Supplementary Table 1 — agreements and disagreements

- ✔ Transposon templates: upstream's s5 and s7 templates match all 24 T5 and all 25 T7
  oligos exactly (computed: every oligo fits `connector + [ACGT]{8} + extended ME`), and
  the counts 24 / 25 agree.
- ✘ **Non-transferred strand.** Upstream writes "Mosaic End (ME) bottom" as
  `/Phos/AGATGTGTATAAGAGACAG` — that is `nextera.ME` itself, the **same** strand as the
  3' end of the transposon oligos, which cannot anneal to it. Supplementary Table 1 gives
  `5Phos/CTGTCTCTTATACACATCT` = `nextera.ME_RC`, the correct complement. Trust the
  table. 🟢 (table) / 🟡 (diagnosis)
- ✘ **Library design.** Supplementary Table 1 and Supp. Fig. 4a describe a **custom
  Nextera** library: PCR with truncated Nextera i5/i7 primers and sequencing with four
  custom primers on a HiSeq 2500. Upstream draws only the **TruSeq connector** library,
  from the Protocol Exchange. Neither connector primer nor the TruSeq P5/P7 primers
  appear in Supplementary Table 1. Both variants are given below; the paper's Fig. 5 title
  suggests the TruSeq one is also in the paper, but we could not read it.
- Upstream's "24 s5 / 25 s7 barcodes, see Table 3 of the Protocol Exchange" matches the
  supplement's counts, so the barcode set is the same 🟡.

## 3. Step by step

Main-text Methods unavailable 🔴. Conditions marked 🟢 come from the Supplementary
Figure 1 legend and Supplementary Table 2, where they are stated as the chosen optimum;
the order of steps follows upstream and the abstract.

1. **Assemble barcoded Tn5**: anneal each T5-N or T7-N with the common annealing primer
   (`5Phos/` ME_RC), then load Tn5 🟢 (Supp. Table 1 header: "annealed with
   TnExtMErev_primer, then assembled with Tn5"). Whether T5 and T7 are loaded separately
   and mixed, or co-loaded, is not stated in our sources 🔴; the Amini design mixes
   separately loaded A- and B-type complexes 🟡.
2. **Fix** cells (the optimisation used "RT 7 min fixed ESCs" 🟢; fixative presumably
   formaldehyde 🔴).
3. **Open chromatin**: 0.3 % SDS, 62 °C 10 min (37 °C 1 h was equivalent) 🟢.
4. **Distribute** — sort single nuclei, or a limited number of nuclei, into wells 🟡
   (upstream; upstream does not mention SDS opening at all, so placing sorting after
   opening is our inference, and the abstract instead describes opening, indexing and
   tagmentation "within a single tube"). The SDS is presumably quenched (e.g. Triton X-100)
   before Tn5 🔴.
5. **Tagment**, one T5/T7 combination per well: 37 °C 1 h 🟢 (an alternative condition tested in Supp. Fig. 1g–h was
   "equivalent", but its values are only in the figure 🔴). Tn5 amount: 3 µL of 12.5 µM for 10,000 cells; per Supp. Table 2, 0.15 µL
   (12.5 µM) or 0.05 µL (37.5 µM) for **1 cell per well** 🟢. Tn5 fills each fragment
   end with T5- or T7-adapter, leaving the 9-nt gap
   ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)).
   Three product types 🟡: T5…T5, T7…T7, and T5…T7. Only T5…T7 is amplified with one
   i5- and one i7-primer.
6. **Pool all wells** 🟡 (upstream), then **sonicate briefly** to release soluble
   chromatin 🟡. The Supp. Fig. 2 legend compared release conditions on 500 cells (without
   shearing, the library needed a few more PCR cycles). From that the authors *reasoned*
   that brief sonication may give better yield and complexity for single-cell itChIP 🟢.
   Where sonication sits relative to pooling is our inference.
7. **Immunoprecipitate** with the antibody of interest (H3K4me3, H3K27ac, H3K27me3,
   Pol II, EZH2, P300 via bio-ChIP in the paper 🟢); **reverse crosslinks, purify DNA**
   🟡 (upstream).
8. **Gap fill and amplify** (Nextera variant): one i5-N50x + one i7-N70x primer per
   pool 🟡. The 9-nt gaps and the non-transferred ME_RC strand must first be filled /
   displaced by a 72 °C extension, as in every Nextera PCR 🟡 (not stated). Cycle numbers
   🔴 (Supp. Fig. 2d shows "PCR cycles" but the values are only in the figure).
   **AMPure XP** size selection after PCR 🟢 (Supp. Fig. 4a legend).
   **TruSeq variant instead of step 8** (upstream 🟡): first PCR with **Connector Primers F + R**, which
   prime on connector A and B and add TruSeq Read 1 / Read 2; purify; second PCR with the
   **Universal P5 primer** and an **Indexed P7 primer**. Upstream explains that same-end
   products (T5…T5, T7…T7) do not amplify because of **semi-suppressive PCR**: their two
   ends are complementary, so each single strand folds into a panhandle that out-competes
   primer binding 🟡 (upstream, citing nanoCAGE).

## 4. Final libraries — 🟡 (assembled from the oligos above)

Segment lists, top strand 5'→3' (the strand that starts with P5). `<T5bc>`, `<T7bc>`:
8-nt Tn5 barcodes; `<i5>`, `<i7>`: 8-nt PCR indices; `ME'`, `X'`: reverse complement.

**Variant A — custom Nextera (Supplementary Table 1, Supp. Fig. 4a)**

```
P5 (29) · <i5> (8) · connector A (21; starts with s5) · <T5bc> (8) · extended ME A (34; ends in ME) · <insert> · extended ME B' (34) · <T7bc>' (8) · connector B' (27; ends in s7') · <i7 as named> (8) · P7' (24)
```

Fixed parts total 29 + 8 + 21 + 8 + 34 + 34 + 8 + 27 + 8 + 24 = **201 bp + insert**.

**Variant B — TruSeq connector (upstream, from the Protocol Exchange)**

```
P5 (29) · TruSeq Read 1 (33) · connector A (21) · <T5bc> (8) · extended ME A (34) · <insert> · extended ME B' (34) · <T7bc>' (8) · connector B' (27) · TruSeq Read 2' (34) · <i7> (6 in upstream's final drawing) · P7' (24)
```

Upstream's final library has the `GT` added and a 6-nt i7.

## 5. Read layout / sequencing

🟢 All samples paired-end **150 bp** (reporting summary); single-cell libraries on an
**Illumina HiSeq 2500** (Supp. Fig. 4a legend).

**Variant A** (custom primers 🟢 from Supp. Table 1; lengths 🟡 computed):

| Read | Primer | Reads | Cycles needed |
|---|---|---|---|
| Read 1 | extended ME A | genomic insert, from the T5 end | 150 (stated) |
| Index 1 | rc(extended ME B) | T7 barcode (rc of the oligo's barcode) · 27 nt connector B' · i7 | 43 |
| Index 2 | grafted P5 | i5 · 21 nt connector A · T5 barcode | 37 |
| Read 2 | extended ME B | genomic insert, from the T7 end | 150 (stated) |

The index cycle numbers actually run, and whether the constant stretches were skipped
with dark cycles (as in sci-ATAC-seq), are not in our sources 🔴.

**Variant B** (upstream 🟡): standard TruSeq primers. Read 1 (150) = connector A (21) ·
T5 barcode (positions 22–29) · extended ME A (34) · insert from position 64. Index 1 =
6-nt i7. Read 2 (150) = connector B (27) · T7 barcode (positions 28–35) · extended ME B
(34) · insert from position 70. Cell = T5 × T7 barcode pair 🟡 (computed positions).

## 6. Open questions

- 🔴 The Methods: fixation, SDS quench, nuclei sorting vs. bulk cell aliquots, IP and
  wash buffers, reverse crosslinking, PCR cycles — all behind the paywall.
- 🔴 The Protocol Exchange text, including whether Connector Primer R and the Indexed P7
  primer really lack the `GT` that completes TruSeq Read 2.
- 🔴 What the "Forward/Reverse i5/i7 PCR primer" (connector A / B 21-mers) are used for.
- 🔴 How the spike-in library is generated and used (normalisation?).
- 🟡 Does the cell barcode also need the PCR i5/i7 (a third/fourth combinatorial level),
  as in sci-ATAC-seq? The sources only say i5/i7 are "introduced during PCR enrichment".

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` fetched all nine Springer supplements and the upstream page, but
not the main text (paywalled) or the Protocol Exchange. The nature.com article page was
fetched by hand with curl (abstract and figure titles only). `tools/scrape_primers.py`
ranked upstream's diagram lines first, so the real oligo table (Supplementary Table 1)
appeared only after many hits of upstream's double-strand diagrams; the i5/i7 primers
were recognised as containing `illumina.P5`, but the truncated s5/s7 tails (14/15 nt)
were not annotated. The T5/T7 oligos were parsed and checked directly from the table
text with `lib/`.
