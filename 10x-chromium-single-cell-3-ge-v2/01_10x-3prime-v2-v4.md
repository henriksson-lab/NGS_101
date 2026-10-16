# 10x Chromium Single Cell 3' Gene Expression v2, v3/v3.1 and v4 (GEM-X)

> **Evidence marking.** 🟢 verbatim from a primary source (10x Genomics technical note /
> user guide, or the paper) · 🟡 derived, inferred, or taken only from the secondary
> upstream scg_lib_structs page · 🔴 not published / not in any source read. Relationships
> marked 🟡 *(computed)* were worked out with `lib/` while writing this note.

**Defining paper (catalogue)** — Zheng GXY, Terry JM, Belgrader P, … Mikkelsen TS,
Hindson BJ, Bielas JH. "Massively parallel digital transcriptional profiling of single
cells." *Nat Commun* 8, 14049 (2017). doi:[10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049),
PMID 28091601, PMC5241818 (CC-BY 4.0).

**Important:** that paper describes the **GemCode / v1** chemistry, not v2–v4 (14-bp
barcode read as the i7 index, cDNA in Read 1, UMI in Read 2 — see §8). The catalogue
lists the same paper for the separate v1 entry
(`10x-chromium-single-cell-3-ge-v1__10.1038+ncomms14049`). Everything about the
**v2–v4 oligos and steps** in this note comes from **10x Genomics' own documents**
(fetched by hand, see below), checked against the upstream scg_lib_structs page.

Vendor documents (primary for the chemistry):

- **CG000108 Rev A** — Technical Note "Assay Scheme and Configuration of Chromium Single
  Cell 3' v2 Libraries" (2017); refers to the v2 user guide CG00052.
- **CG00052 Rev D** — original 10x Single Cell 3' Reagent Kits v2 User Guide (2018),
  checked for the v2 thermal programs and sequencing parameters. Its revision is later
  than CG000108 Rev A, so exact 2017 conditions remain revision-dependent.
- **CG000183 Rev A** — Chromium Single Cell 3' Reagent Kits **v3** User Guide.
- **CG000731 Rev A** — Chromium **GEM-X** Single Cell 3' Reagent Kits **v4** User Guide.

Sources read (in `$CHEM_DATA/sources/10x-chromium-single-cell-3-ge-v2-v4__10.1038+ncomms14049/`
and `$CHEM_DATA/sources/10x-chromium-single-cell-3-ge-v2/`, never committed):

| File | What | Used for |
|---|---|---|
| `upstream_10xChromium3.html(.txt)` | scg_lib_structs drawn page | orientation only; not molecular evidence |
| `upstreamdata_CG000108_AssayConfiguration_SC3v2.pdf(.txt)` | 10x technical note, v2 | **every v2 oligo**, verbatim (fetched by hand from the upstream site's `data/10X-Genomics/`) |
| `CG00052_SingleCell3_ReagentKitv2UserGuide_RevD.pdf(.txt)` | official 10x v2 Rev D user guide | **v2 thermal programs, amplification tiers and 26/8/98 run**, sections 1.5, 2.2, 3.1, 4.2 |
| `upstreamdata_CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf(.txt)` | 10x v3 user guide | **v3 oligos** (appendix "Oligonucleotide Sequences"), thermal programs, cycle numbers, run parameters |
| `upstreamdata_CG000731_ChromiumGEM-X_SingleCell3_ReagentKits_v4_UserGuide_RevA.pdf(.txt)` | 10x v4 (GEM-X) user guide | **v4 bead oligo and final library**, thermal programs, run parameters |
| `10x_CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf(.txt)` | 10x v3.1 single-index guide, appendix pp. 57–58 | v3.1 bead, TSO, cDNA primers, adapter, 8-nt i7 single-index primer and library |
| `10x_CG000315_ChromiumNextGEMSingleCell3-_GeneExpression_v3.1_DualIndex__RevE.pdf(.txt)` | 10x v3.1 dual-index guide, appendix pp. 61–63 | v3.1 bead, TSO, cDNA primers, adapter, 10-nt i5/i7 TT primers and library |
| [10x Dual Index Kit TT Set A index map](https://cdn.10xgenomics.com/raw/upload/v1655151898/support/in-line%20documents/Dual_Index_Kit_TT_Set_A.json) | official vendor JSON, linked by the [10x index-set page](https://www.10xgenomics.com/support/spatial-gene-expression-fresh-frozen/documentation/steps/library-construction/sample-index-sets-for-visium-spatial-gene-expression) | 96 named well-to-i7/i5 pairs and both i5 sequencing workflows; checked SI-TT-A1 below |
| [10x Chromium i7 Multiplex index map](https://cdn.10xgenomics.com/raw/upload/v1655155614/support/in-line%20documents/chromium-shared-sample-indexes-plate.json) and [Single Index Kit T Set A map](https://cdn.10xgenomics.com/raw/upload/v1655155616/support/in-line%20documents/Single_Index_Kit_T_Set_A.json) | vendor JSON, linked by the [single-cell 3′ index-set page](https://www.10xgenomics.com/support/universal-three-prime-gene-expression/documentation/steps/library-prep/sample-index-sets-for-single-cell-3) | v2/v3 and v3.1 single-index well assignments; each well contains four distinct 8-nt i7 oligos |
| `ncomms14049_PMC5241818.xml(.txt)` | the paper, JATS full text (Europe PMC) | v1 workflow and conditions (contrast only); `.txt` made by hand, see §9 |
| `ncomms14049_supp_ncomms14049-s1.pdf(.txt)` (= `…MOESM828_ESM.pdf`) | Supplementary Information | figure legends only; no oligos, no methods |
| `…-s2..s5.xlsx` (= `MOESM829..832`) | Supplementary Data 1–4 | per-barcode tables, gene lists; not chemistry |
| `ncomms14049_PMC5241818.html(.txt)` | PMC article page | **not usable** — the fetch returned a reCAPTCHA page |

Not fetched / not available:

| What | Where | Why it matters |
|---|---|---|
| CG00052 revision contemporaneous with CG000108 Rev A | 10x historical support archive | Rev D is available and checked; earlier differences remain unverified 🟡 |

---

## 1. What it is

Droplet (GEM) single-cell 3' RNA-seq. Each gel bead carries millions of copies of one
oligo: **partial TruSeq Read 1 · 16-nt 10x barcode · UMI · poly(dT)**. In the droplet the
bead dissolves, the cell lyses, and the bead oligo primes RT on poly(A) mRNA; the RT
**template-switches** onto a TSO, so every first-strand cDNA is flanked by Read 1 at one
end and the SMART handle at the other. 🟢 (CG000108 introduction; CG000183/CG000731
introduction.) Emulsion broken, cDNA amplified in bulk with two primers, then **enzymatically
fragmented, end-repaired, A-tailed** and ligated to a **Read 2 adapter**; a sample-index PCR
adds P5/P7. Only the 3'-end fragment (the one carrying Read 1 + barcode + UMI) has the P5
side primer site, so the library is 3'-end counting. 🟢 workflow, 🟡 the "only the 3'
fragment amplifies" reasoning.

| | New relative to v1 (the paper) | Recurs in |
|---|---|---|
| 1 | Barcode + UMI moved **into Read 1** (TruSeq Read 1 on the bead), cDNA read in **Read 2** | every later 10x 3' kit; most bead-based methods (Drop-seq puts them in Read 1 too) |
| 2 | 16-nt barcode (v1: 14), sample index in **i7** (v1: i5) | — |
| 3 | **Enzymatic fragmentation + end repair + A-tail + T/A ligation** of a Y-less, short-bottom-strand Read 2 adapter (paper's v1 used Covaris shearing) | standard TruSeq library prep |
| 4 | TSO = the SMART/ISPCR handle + `ACAT` + `rGrGrG` — the SMART-seq TSO with three ribo-G | [template switching](../ref/concepts/template-switching.md), `rt.SMART_HANDLE` |

Version differences (🟢 per the respective vendor document unless marked):

| | v2 (CG000108) | v3 (CG000183) | v4 GEM-X (CG000731) |
|---|---|---|---|
| UMI | 10 nt | 12 nt | 12 nt |
| poly(dT) on bead, as written | `T30VN` | `T30` (appendix) but `Poly(dT)VN` in the overview figure | `T30VN` |
| second bead oligo (Capture Seq 1, Feature Barcode) | — | yes, "not used for 3' Gene Expression" | yes ("Capture Sequence 1") |
| cDNA reverse primer | 27 nt (handle + `ACAT`) | 22 nt ("Partial TSO") | CG000731 lists PN-2000089 but omits bases; same PN is printed as 22 nt in CG000315 Rev F 🟡 |
| Read 2 adapter top strand | 33 nt, ends `…CAGTCAC` | 32 nt, ends `…CAGTCA` | 🔴 not given |
| sample index | single, i7 8 nt | single, i7 8 nt | **dual**, i7 10 nt + i5 10 nt (Dual Index Kit TT Set A) |
| GEM-RT | 53 °C 45 min, 85 °C 5 min 🟢 (CG00052 Rev D §1.5) | 53 °C 45 min, 85 °C 5 min | **48 °C** 45 min, 85 °C 5 min |
| run (R1 / i7 / i5 / R2) | 26 / 8 / 0 / 98 🟢 (CG00052 Rev D §4.2) | 28 / 8 / 0 / 91 | 28 / 10 / 10 / 90 |

**Disagreement with upstream.** The scg_lib_structs page says the v4 library structure
"is exactly the same as v3 and v3.1". CG000731 shows otherwise: the v4 library is
**dual-indexed** (P5 · i5 N10 · full TruSeq Read 1 …), and the i7 is 10 nt, not 8. 🟢
(CG000731 appendix, "Gene Expression Library Sample Index PCR Product".) Upstream's v3/v4
final library is therefore right only for single-index v3 (and v3.1 single index). The
upstream page also writes the v3/v4 bead as `(T)30VN`; CG000183's appendix writes v3 as
plain `T30` (its own overview figure says `VN`) — an inconsistency inside the vendor
document itself. 🟢 both readings.

## 2. Oligos

All sequences 🟢 **verbatim** from the vendor documents, copied with their dashes and
spacing exactly; `rG` = ribo-G, `N` = barcode/UMI/sample-index position. Part numbers as
printed. Each line was located in the source text with `grep -F`.

### v2 (CG000108, Fig. 4)

```
Gel Bead Oligo Primer (PN-220104)
  5’-CTACACGACGCTCTTCCGATCT-NNNNNNNNNNNNNNNN-NNNNNNNNNN-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3’
RT Primer (TSO) (PN-310354)
  5’-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3’
cDNA Primer Mix (PN-22016)  Forward primer
  5’-CTACACGACGCTCTTCCGATCT-3’
                            Reverse primer
  5’-AAGCAGTGGTATCAACGCAGAGTACAT-3’
Adapter (Read 2) (PN-220026), double-stranded
  5’-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3’
  3’-TCTAGCCTTCTCG-5’
Sample Index PCR, Forward primer: P5 – Partial Read 1 (PN-220111)
  5’-AATGATACGGCGACCACCGA-GATCTACACTCTTTCCCTACACGACGCTC-3’
Sample Index PCR, Reverse primer: P7 – Sample Index – Partial Read 2 (PN-220103)
  5’-CAAGCAGAAGACGGCATACGAGAT-NNNNNNNN-GTGACTGGAGTTCAGACGTGT-3’
```

### v3 (CG000183 appendix, "Oligonucleotide Sequences")

```
Gel Bead Primers (TruSeq Read 1 · 10x Barcode · UMI · Poly(dT))
  5’-CTACACGACGCTCTTCCGATCT-NNNNNNNNNNNNNNNN-NNNNNNNNNNNN-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-3’
Gel Bead, Nextera Read 1 (Read 1N) · Barcode · UMI · Capture Seq 1 — "Primer not used for 3’ Gene Expression"
  5’-GTCAGATGTGTATAAGAGACAG-NNNNNNNNNNNNNNNN-NNNNNNNNNNNN-TTGCTAGGACCGGCCTTAAAGC-3’
Template Switch Oligonucleotide (PN-3000228)
  5’-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3’
cDNA Primers (PN-2000089)  Forward, Partial Read 1
  5’-CTACACGACGCTCTTCCGATCT-3’
                           Reverse, Partial TSO
  5’-AAGCAGTGGTATCAACGCAGAG-3’
Adaptor Oligos (PN-2000094), Partial Read 2
  5’- GATCGGAAGAGCACACGTCTGAACTCCAGTCA-3’
  3’-TCTAGCCTTCTCG-5’
Sample Index PCR Primer (PN-2000095)
  5’-AATGATACGGCGACCACCGAGATCT-ACACTCTTTCCCTACACGACGCTC-3’
Chromium i7 Sample Index (PN-220103)
  5’-CAAGCAGAAGACGGCATACGAGAT-NNNNNNNN GTGACTGGAGTTCAGACGTGT-3’
```

### v4 GEM-X (CG000731 appendix)

```
Gel Bead Primers
  5’-CTACACGACGCTCTTCCGATCT-N16-N12-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3’
```

CG000731 gives no TSO, cDNA-primer, adaptor or full index-primer sequences — only the
bead oligo and final product (§5). Its kit contents (p. 9) identify **cDNA Primers
PN-2000089**. CG000315 Rev F, "Oligonucleotide Sequences" (p. 62), prints the same
part number with forward `CTACACGACGCTCTTCCGATCT` and reverse
`AAGCAGTGGTATCAACGCAGAG`. 🟢 Each document's own claim; 🟡 assigning those bases to
v4 on the assumption that a stable part number means a stable oligo formulation.
The v4 TSO and adaptor sequences remain unresolved. The TT index **read** bases are
available in the manufacturer CSV, but the full index-primer oligos are not established.

### Sequencing primers (upstream only)

Upstream lists the standard Illumina TruSeq Read 1 / Read 2 primers and an i7 "Sample
Index sequencing primer" `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC`. They are
`illumina.TRUSEQ_READ1`, `illumina.TRUSEQ_READ2` and `illumina.INDEX1_PRIMER` exactly. 🟡
(computed); not listed in the vendor documents read — the user guides only name "TruSeq
Read 1 / Read 2" as the primer sites. 🟢 names only.

### Upstream vs vendor, oligo by oligo

Every upstream oligo matches the vendor text base for base (v2 against CG000108; v3
"cDNA reverse primer" and 32-nt adapter against CG000183). Differences are only in
annotation: upstream writes the v3/v4 bead as `(T)30VN` (CG000183 appendix: `T30`), and
upstream does not show the v3/v4 Capture Seq 1 bead oligo on this page. 🟡

## 3. How the oligos interlock — 🟡 (computed)

With `lib/` (`chemdraw.revcomp`, `illumina`, `rt`, `nextera`):

- **Bead handle = `illumina.TRUSEQ_READ1[11:]`** — the 3'-most 22 nt of the 33-nt TruSeq
  Read 1 primer site ("Partial Read 1"). The cDNA forward primer *is* this handle.
- **TSO = `rt.SMART_HANDLE` + `ACAT` + `rGrGrG`** (27 nt DNA + 3 ribo-G). Its DNA part is
  the SMART-seq TSO handle; the SMART-seq2 TSO has the same bases but ends `rGrG+G`
  (LNA). 🟡
- **v2 cDNA reverse primer = the TSO without its rG's** (27 nt; `SMART_HANDLE` + `ACAT`).
  **v3 reverse primer = `SMART_HANDLE[:22]`** — one base shorter than the ISPCR primer
  and a prefix of the v2 one.
- **Read 2 adapter**: the v2 top strand (33 nt) **equals `illumina.INDEX1_PRIMER`**, and
  its reverse complement is the first 33 nt of `illumina.TRUSEQ_READ2`. It starts with
  `illumina.STEM` (`GATCGGAAGAGC`). The v3 top strand is the same minus the final `C`
  (32 nt). The bottom strand, read 5'→3', is `GCTCTTCCGATCT` (13 nt): its reverse
  complement is `AGATCGGAAGAGC` = `A` + the first 12 nt of the top strand. So the bottom
  strand pairs with the top strand's 12-nt stem and puts a **3'-T overhang** opposite the
  insert's A-tail — a short, partly double-stranded adapter, not a Y-adapter.
- **Sample-index forward primer = `illumina.NEBNEXT_UNIVERSAL_PRIMER[:49]`**, i.e.
  `illumina.P5` + `TRUSEQ_READ1[4:24]`. P5 ends `ACAC` and TruSeq Read 1 starts `ACAC`,
  so the two **overlap by 4 nt** in the single-index library. Its 3' end overlaps the bead
  handle (Partial Read 1) by **13 nt** (`CTACACGACGCTC`), which is what it primes on.
- **Sample-index reverse primer = `illumina.P7` + i7 (8 N) + `TRUSEQ_READ2[:21]`**, 53 nt.
  Its reverse complement is the 3' end of the final library; `revcomp(TRUSEQ_READ2[:21])`
  is the top-strand adapter from position 12 on, i.e. it primes on the ligated adapter's
  3' part.
- **v3/v4 Capture Seq 1 bead oligo**: its 22-nt 5' handle = `nextera.ADAPTOR_S5[11:]`
  (the 3' 3 nt of the s5-specific part, `GTC`, + the full 19-nt mosaic end `nextera.ME`) — "Partial Nextera Read 1".
  The 3' Capture Sequence 1 (`TTGCTAGGACCGGCCTTAAAGC`) does not match anything in `lib/`.
  Not used by the gene-expression library. 🟡

## 4. Step by step

Conditions below are revision-specific: v2 from CG00052 Rev D, v3 from CG000183,
and v4 from CG000731. CG000108 Rev A names the v2 steps (1.5 GEM-RT, 2.2 cDNA
amplification, 3.3 adapter ligation, 3.5 sample-index PCR), but does not print
their thermal programs.

1. **GEM generation.** Cells, RT master mix (with TSO), gel beads and partitioning oil
   meet on the chip; the bead dissolves, its oligos are released and the cell lyses. 🟢
2. **GEM-RT.** v2 Rev D and v3: 53 °C 45 min, 85 °C 5 min, 4 °C.
   v4: **48 °C** 45 min, 85 °C 5 min. 🟢
   The bead oligo primes at the poly(A) tail; at the mRNA 5' end the RT adds non-templated
   C's and switches to the TSO. First strand, 5'→3':
   `Partial Read 1 · <cbc16> · <umi> · T30VN · <cDNA> · CCCATGTACTCTGCGTTGATACCACTGCTT`
   🟢 (CG000108 "Reverse Transcript Product"; CG000183 likewise with 12-N UMI and `T30`).
   The final 30 nt are `CCC` + revcomp(TSO DNA part). 🟡 (computed) —
   [template switching](../ref/concepts/template-switching.md),
   [reverse transcription](../ref/concepts/reverse-transcription.md).
3. **Break GEMs, Dynabeads MyOne SILANE cleanup** of the single-stranded cDNA. 🟢
4. **cDNA amplification** with the forward (Partial Read 1) and reverse (TSO-handle)
   primers: v2 Rev D uses 98 °C 3 min; cycles of 98 °C 15 s, **67 °C** 20 s,
   72 °C 1 min; 72 °C 1 min. Its recommended cycles are 14 for <2,000 cells,
   12 for 2,000–6,000, 10 for 6,000–10,000 and 8 for >10,000. 🟢
   v3 uses 98 °C 3 min; cycles of 98 °C 15 s, 63 °C 20 s, 72 °C 1 min; 72 °C 1 min.
   v4: 98 °C 45 s; 98 °C 20 s, 63 °C 30 s, 72 °C 1 min; 72 °C 1 min. Total cycles by
   targeted cell number: 13 (<500), 12 (500–6,000; v4 writes 501–6,000), 11 (>6,000) in
   both. 🟢 Then
   SPRIselect cleanup. Product: double-stranded
   `Partial Read 1 · <cbc16> · <umi> · T30(VN) · <cDNA> · CCC · TSO-handle'` 🟢
5. **Fragmentation, end repair, A-tailing** (one tube): pre-cooled 4 °C block, 32 °C
   5 min, 65 °C 30 min (v3 and v4). 🟢 Then double-sided SPRIselect size selection. Three
   kinds of fragment result 🟡 (as upstream reasons, and as follows from the primers):
   - 3'-end fragment: `Partial Read 1 · <cbc16> · <umi> · T30VN · <cDNA>` + A-tail;
   - internal cDNA fragments, A-tailed at both ends;
   - 5'-end fragment carrying the TSO handle.
6. **Adaptor ligation**: 20 °C 15 min (v3, v4). 🟢 The T-overhang Read 2 adapter is
   ligated to every A-tailed end. 🟡
7. **Post-ligation SPRIselect cleanup.** 🟢
8. **Sample-index PCR**: 98 °C 45 s; cycles of 98 °C 20 s, 54 °C 30 s, 72 °C 20 s;
   72 °C 1 min (v3, v4); cycles by cDNA input (v3: 14–16 for 1–25 ng, 12–14 for
   25–150 ng; v4: 14–16 for 0.25–50 ng, 12–14 for 50–250 ng, 10–12 for 250–600 ng, 8–10
   for 600–1,100 ng, 6–8 for 1,100–1,500 ng, 5 for >1,500 ng). 🟢 Only fragments with the bead's Partial Read 1 on one end and the adapter on
   the other are amplified exponentially: the P5 primer needs the Read 1 handle, the P7
   primer needs the adapter. Internal fragments carry the adapter on both ends and only
   the P7-side primer site; 5'-end fragments carry the TSO handle, for which no primer is
   present. 🟡 Upstream assumes the TSO / Read 1 5' ends are blocked so that the adapter
   ligates only to cDNA ends — not stated in the vendor documents. 🔴
9. **Double-sided SPRIselect** size selection, QC, sequencing. 🟢

## 5. Final library — 5'→3', top strand

v2, single index 🟢 (CG000108 "Sample Index PCR Product", as written):

```
5’-AATGATACGGCGACCACCGA-GATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-NNNNNNNNNNNNNNNN-NNNNNNNNNN-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-cDNA_Insert-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC-NNNNNNNN-ATCTCGTATGCCGTCTTCTGCTTG-3’
```

As segments 🟡 (computed against `lib/`):

```
v2:      P5 · TRUSEQ_READ1[4:] · <cbc 16> · <umi 10> · T30VN · <cDNA> · A · INDEX1_PRIMER · <i7' 8> · P7_RC
v3/v3.1 single index: P5 · TRUSEQ_READ1[4:] · <cbc 16> · <umi 12> · T30(VN) · <cDNA> · A · INDEX1_PRIMER · <i7' 8> · P7_RC
v3.1 dual index/v4: P5 · <i5 10> · TRUSEQ_READ1 · <cbc 16> · <umi 12> · T30VN · <cDNA> · A · INDEX1_PRIMER · <i7' 10> · P7_RC
```

- `P5 · TRUSEQ_READ1[4:]` (58 nt) is exactly `illumina.NEBNEXT_UNIVERSAL_PRIMER`; the
  right-hand end `A · INDEX1_PRIMER · <i7'> · P7_RC` is 66 nt (v2/v3) or 68 nt (v4).
  Fixed length besides the cDNA insert: **182 nt (v2), 184 nt (v3), 200 nt (v4)** —
  58 + 16 + UMI + 32 (T30VN) + 66, and for v4 29 + 10 + 33 + 16 + 12 + 32 + 68. 🟡
  (computed; T30 without VN would make v3 182 nt.)
- `A` is the A-tail added to the insert, not part of the adapter oligo: CG000108 writes
  the junction as `cDNA_Insert-AGATCGG…` (the adapter oligo itself starts `GATCGG…`). 🟡
- `<i7'>` is the top-strand copy of the N8/N10 in the P7 primer, i.e. the reverse
  complement of the primer's index as written. 🟡 (computed: the P7 primer's reverse
  complement is the 3' 53 nt of the v2 right-hand end).
- 🟢 For TT Set A well SI-TT-A1, the 10x index map gives i7
  `GTAACATGCG`, i5 workflow A `AGTGTTACCT` and i5 workflow B
  `AGGTAACACT` (its reverse complement). The map gives all 96 wells;
  these are sequencing/sample-sheet index sequences, whereas CG000315 prints N10
  placeholders in the synthesis-primer bodies.
- 🟢 For v2/v3 single-index wells, the vendor specifies **four** 8-nt i7 oligos
  per well, not one fixed index. Both named index maps list SI-GA-A1 as
  `GGTTTACT`, `CTAAACGG`, `TCGGCGTC`, `AACCGTAA`. The `<i7' 8>` segment above
  represents one member of that well's four-member mixture; the complete per-well
  list is in the linked 10x files.
- CG000315 and CG000731 show the 10-nt i5 **between P5 and the full 33-nt TruSeq Read 1**
  (`…GATCTACAC-N10-ACACTCTTTCCC…`) for v3.1 dual index and v4, so the 4-nt P5/Read 1 overlap of the single-index versions is gone. 🟢
  (verbatim) / 🟡 (interpretation).

v4 final product, verbatim 🟢 (CG000731 appendix):

```
5’-AATGATACGGCGACCACCGAGATCTACAC-N10-ACACTCTTTCCCTACACGACGCTCTTCCGATCT-N16-N12 -T30-VN-cDNA_InsertAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC-N10-ATCTCGTATGCCGTCTTCTGCTTG-3’
```

The cDNA insert on the top strand is **antisense** to the mRNA (the top strand is the
bead-primed first strand). 🟡

## 6. Read layout / sequencing

| Read | Primer | Reads | Cycles |
|---|---|---|---|
| Read 1 | TruSeq Read 1 (`illumina.TRUSEQ_READ1`) on the bottom strand | `<cbc16>` + `<umi>` | v2 26 🟢 (CG00052 Rev D §4.2); v3 28 🟢; v4 28 🟢 |
| i7 | `illumina.INDEX1_PRIMER` = the adapter top strand 🟡 | sample index | v2 8 🟢 (CG00052 Rev D §4.2); v3 8 🟢; v4 10 🟢 |
| i5 | P5-side (`illumina.INDEX2_PRIMER_RC` or flow-cell primed, instrument dependent) | v4 only | 0 (v2/v3) / 10 (v4) 🟢 |
| Read 2 | TruSeq Read 2 (`illumina.TRUSEQ_READ2`) | cDNA, **sense** to the mRNA, starting at the fragmentation end 🟡 | v2 98 🟢 (CG00052 Rev D §4.2); v3 91 🟢; v4 90 🟢 |

- v2 Rev D: paired-end, single indexing, 26 / 8 / 0 / 98. v3: paired-end, single indexing,
  28 / 8 / 0 / 91. v4: paired-end, dual indexing, 28 / 10 / 10 / 90. 🟢
  28 / 10 / 10 / 90, minimum 20,000 read pairs per cell. 🟢
- Read 1 primer's 3' end = the bead handle, so Read 1 base 1 is the first barcode base;
  26 = 16 + 10 (v2), 28 = 16 + 12 (v3/v4). 🟡 (computed)
- `rc(TRUSEQ_READ2)` = `A` + the 33-nt adapter top strand: the Read 2 primer's 3' T sits on
  the A-tail, so Read 2 starts with the last insert base before the A-tail. 🟡 (computed)
- The 10x documents say the 16-bp barcode and the UMI "are encoded in Read 1, while Read 2
  is used to sequence the cDNA fragment" and the sample index is "the i7 index read"
  (v3) / "i7 and i5 … sample index reads" (v4). 🟢 (short quotes).

Barcode inclusion lists: the [10x Cell Ranger algorithm guide](https://www.10xgenomics.com/support/software/cell-ranger/7.2/algorithms-overview/cr-gex-algorithm)
assigns approximately 737,000 allowed 16-base barcodes to v2 and approximately 3
million to v3/v3.1. The original 10x `737K-august-2016.txt` in its
[Supernova tenkit repository](https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/barcodes/737K-august-2016.txt)
has 737,280 16-base rows; its SHA-256 matches the decompressed cached scg copy
exactly (`b0dda4b114d8fc7ea8def94fb472e6205c1831a72e3ac245ed13cc332d9f20d1`).
The cached scg `3M-february-2018.txt.gz` has 6,794,880 16-base rows
(streamed count). That file's line count does not match the
guide's approximate v3 inclusion-list size or [CG000201 Rev A's](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000201_TechNote_Chromium_Single_Cell_3___v3_Reagent__Workflow___Software_Updates_RevA.pdf)
3.5-million *gel-bead diversity* figure; those are different claims and the
mirror's provenance and precise relation to the v3 reagent remain unverified.
The v4 mirror is named `3M-3pgex-may-2023.txt.gz`; its members have not been
traced to an original manufacturer file here. A barcode inclusion list does
not assign a particular barcode to a physical bead or well.

## 7. Open questions

- 🟡 **v2 revision boundary**: CG00052 Rev D gives the programs and cycle tiers
  above, but CG000108 Rev A is earlier. Exact changes between these revisions
  have not been established. The v2 model uses CG000108 for molecular structure.
- 🟢 **v3.1** CG000204 Rev D and CG000315 Rev E were subsequently checked. Both print
  the 16-nt-barcode/12-nt-UMI `T30VN` bead primer, TSO
  `AAGCAGTGGTATCAACGCAGAGTACATrGrGrG`, cDNA primers and adaptor. CG000204 has
  one 8-nt i7 index; CG000315 prints P5–N10 i5 and P7–N10 i7 TT primers and the
  corresponding dual-index final product. The v3.1 dual-index arrangement above is
  therefore source-supported, rather than an inference from v4.
- 🔴 **v4 TSO and adaptor**: not printed in CG000731. Whether the v4 TSO is still the v3
  sequence is not documented in what was read. The guide does list cDNA Primers PN-2000089;
  CG000315 Rev F p. 62 prints both sequences under that same PN, giving a conditional
  cross-guide assignment rather than direct v4 sequence disclosure. 🟡
- 🔴 **v3 bead `T30` vs `T30VN`**: CG000183 appendix omits VN, its own overview figure
  shows `Poly(dT)VN`; v2 and v4 both write VN.
- 🔴 **v3 adapter 32 vs 33 nt**: the v3 adapter top strand ends `…CAGTCA`, but CG000183's
  own ligation product shows `…CAGTCAC`. The final `C` is in any case restored by the
  sample-index reverse primer, whose reverse complement ends `…CAGTCAC` — so the final
  library is unchanged. 🟡 (computed)
- 🔴 Whether the TSO and the Read 1 handle are 5'-blocked (upstream's assumption to
  explain why only cDNA ends take the adapter).
- 🟡 i5 read primer for v4: depends on the instrument workflow (forward vs reverse
  complement); not stated in CG000731 beyond "dual indexing".

## 8. The paper's own chemistry (v1 / GemCode), for contrast

From the paper (`ncomms14049_PMC5241818.xml.txt`), 🟢 paraphrased:

- Bead oligo: sequencing adapters and primers, a **14-bp barcode** from ~750,000 designed
  sequences, a **10-bp UMI** (some early libraries 5 bp), and an **anchored 30-bp
  oligo-dT** (lines 42–44). cDNA ends with the template-switching oligo; PCR uses primers
  complementary to the switch oligo and sequencing adapters (lines 55–57).
- GEM-RT 55 °C 2 h, 85 °C 5 min; Silane + 0.6× SPRI cleanup; cDNA PCR 98 °C 3 min,
  14 × (98 °C 15 s, 67 °C 20 s, 72 °C 1 min), 72 °C 1 min; **Covaris shearing to ~200 bp**;
  end repair + A-tailing, adapter ligation, SPRI, sample-index PCR (lines 300–307).
- NextSeq 500: **98 bp Read 1 (cDNA), 14 bp I7 (barcode), 8 bp I5 (sample index), 10 bp
  Read 2 (UMI)** (lines 58–59, 309–310).

So v2 reorganised the whole read layout (barcode + UMI from I7/Read 2 into Read 1, sample
index from I5 to I7, cDNA from Read 1 to Read 2) and replaced shearing by enzymatic
fragmentation. No oligo sequence is printed in the paper or its supplement. 🟢 (absence
checked by `tools/scrape_primers.py`: only ddPCR probes and barcode lists were found).

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the upstream page, the Europe PMC XML and all supplements
(twice: Europe PMC zip and Springer, identical by md5). The PMC HTML page came back as a
**reCAPTCHA stub** (195-byte `.txt`). `tools/doctext.py` on the JATS XML produced a
`.txt` that is a byte-for-byte copy of the XML on a single line (no tag stripping), so it
was replaced by a tag-stripped, line-folded text made by hand. The paper has no v2–v4
sequences at all; the vendor PDFs that upstream itself cites (CG000108, CG000183,
CG000731) were linked from the upstream page but **not fetched by get_sources** — they
were downloaded by hand from the scg_lib_structs site into the same directory and
converted with `doctext.py`, and `scrape_primers.py` then found every oligo in them.
`scrape_primers.py --find` reports 0 locations for oligos that the vendor writes with
`-` segment separators and N runs (e.g. the bead oligo, the P7 index primer), even though
the as-written strings are there; those were confirmed with `grep -F` on the exact
as-written text instead. The scraper's top hits on the supplements were 14-nt cell
barcodes from per-cell tables (expected false positives).


## Source-led model corrections (2026-10-16)

🟢 Both the cDNA-amplification and sample-index-PCR primer architectures are
printed in the version-specific manufacturer source (v2: CG000108 Rev A Fig. 4;
v3: CG000183 Rev A Appendix pp. 54–55). The schematic now uses these sequences
in its PCR geometry instead of claiming that they were unavailable. Copied TSO
bases are derived from the printed TSO, and the post-template-switch state
includes that handle before recovery and amplification. These are repository
model corrections.

For v3, CG000183 Rev A Appendix p. 54 prints a 32-base adaptor ending `CAGTCA`,
while its ligation-product drawing ends `CAGTCAC`. The public page now shows
this source discrepancy. The model retains the printed 32-base input; the
sample-index reverse primer supplies the terminal C in the completed library.
The guide does not resolve which pre-PCR drawing reflects the supplied oligo.


## Followed manufacturer references (2026-10-16)

🟢 The source review followed CG00052 Rev D to CG000089 Rev B, CG000085 Rev A,
CG000120 Rev A and CG000061 Rev A, plus the original Chromium i7 Multiplex
index map. CG000089 and CG000085 explicitly use 26/8/98. CG000120 Table 2
uses 26/8/91 on the NovaSeq 100-cycle kit, whose available cycle budget is
125 cycles. The shorter insert read is an explicit run variant rather than
a sequence discrepancy. CG000085 further links CG000080 Rev B: its Table 1
prints both 26/8/98 and an alternative 150/8/150 run; p. 4 says Cell Ranger
does not use R1 cycles beyond 26. The displayed default remains CG00052
Rev D's run.

🟢 The original manufacturer `chromium-shared-sample-indexes-plate.json`
contains 96 named wells with four eight-base i7 members per well. SI-GA-A1
contains `GGTTTACT`, `CTAAACGG`, `TCGGCGTC`, `AACCGTAA`, matching the named
well in the separately downloaded Single Index Kit T map. The model's N8
index denotes one member, not the entire well mixture concatenated.

CG000061 studies a **Genome** library made with Genome Reagent Kits v2
(CG00043); it explains general double-sided SPRI selection and does not
replace the RNA-specific size-selection instructions in CG00052. Its linked
CG00043/CG00047/CG00048 concern Genome library preparation and QC and are
outside the RNA construct claim. CG00053 Rev C and the image-only CG000126
Rev A flowchart were also inspected for their input-cell preparation scope;
they concern suspension quality, washing, counting and loading before the
cell input represented by this schematic. They are not additional sources
for the barcode-bearing library's fixed oligo sequences. Sample-specific
thawing, counting-instrument and viability-stain manuals were not imported
into the molecular model. Full source locations and independent-review
status are recorded in the source-review checklist.
