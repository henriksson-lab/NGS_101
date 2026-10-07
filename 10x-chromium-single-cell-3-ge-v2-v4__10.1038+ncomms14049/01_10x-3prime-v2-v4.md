# 10x Chromium Single Cell 3' Gene Expression v2, v3/v3.1 and v4 (GEM-X)

> **Evidence marking.** 🟢 verbatim from a primary source (10x Genomics technical note /
> user guide, or the paper) · 🟡 derived, inferred, or taken only from the secondary
> upstream scg_lib_structs page · 🔴 not published / not in any source read. Relationships
> marked 🟡 *(computed)* were worked out with `lib/` while writing this note; they are not
> yet asserted in a self-test, because this protocol has no `tools/` module yet (status
> `notes` in `catalogue/ours.tsv`).

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
- **CG000183 Rev A** — Chromium Single Cell 3' Reagent Kits **v3** User Guide.
- **CG000731 Rev A** — Chromium **GEM-X** Single Cell 3' Reagent Kits **v4** User Guide.

Sources read (in `$CHEM_DATA/sources/10x-chromium-single-cell-3-ge-v2-v4__10.1038+ncomms14049/`,
never committed):

| File | What | Used for |
|---|---|---|
| `upstream_10xChromium3.html(.txt)` | scg_lib_structs page "10x Chromium Single Cell 3' Gene Expression" | second source: oligos, steps, final library, read layout; checked against the vendor PDFs |
| `upstreamdata_CG000108_AssayConfiguration_SC3v2.pdf(.txt)` | 10x technical note, v2 | **every v2 oligo**, verbatim (fetched by hand from the upstream site's `data/10X-Genomics/`) |
| `upstreamdata_CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf(.txt)` | 10x v3 user guide | **v3 oligos** (appendix "Oligonucleotide Sequences"), thermal programs, cycle numbers, run parameters |
| `upstreamdata_CG000731_ChromiumGEM-X_SingleCell3_ReagentKits_v4_UserGuide_RevA.pdf(.txt)` | 10x v4 (GEM-X) user guide | **v4 bead oligo and final library**, thermal programs, run parameters |
| `ncomms14049_PMC5241818.xml(.txt)` | the paper, JATS full text (Europe PMC) | v1 workflow and conditions (contrast only); `.txt` made by hand, see §9 |
| `ncomms14049_supp_ncomms14049-s1.pdf(.txt)` (= `…MOESM828_ESM.pdf`) | Supplementary Information | figure legends only; no oligos, no methods |
| `…-s2..s5.xlsx` (= `MOESM829..832`) | Supplementary Data 1–4 | per-barcode tables, gene lists; not chemistry |
| `ncomms14049_PMC5241818.html(.txt)` | PMC article page | **not usable** — the fetch returned a reCAPTCHA page |

Not fetched / not available:

| What | Where | Why it matters |
|---|---|---|
| CG00052, v2 user guide (thermal programs, cycle numbers for v2) | 10x support site, `https://www.10xgenomics.com/support` (login/JS) | v2 conditions are only inferred from v3 here 🔴 |
| CG000204 / CG000315, v3.1 single- and dual-index user guides | 10x support site | v3.1 is named by upstream but no v3.1 document was read 🔴 |
| Dual Index Kit TT Set A primer sequences (v4, v3.1 dual) | 10x "Dual Index Kit TT Set A" sample-index sheet | v4 i5/i7 primers are given only as part of the final library 🔴 |

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
| cDNA reverse primer | 27 nt (handle + `ACAT`) | 22 nt ("Partial TSO") | 🔴 not given in CG000731 |
| Read 2 adapter top strand | 33 nt, ends `…CAGTCAC` | 32 nt, ends `…CAGTCA` | 🔴 not given |
| sample index | single, i7 8 nt | single, i7 8 nt | **dual**, i7 10 nt + i5 10 nt (Dual Index Kit TT Set A) |
| GEM-RT | 🔴 (v2 guide not read) | 53 °C 45 min, 85 °C 5 min | **48 °C** 45 min, 85 °C 5 min |
| run (R1 / i7 / i5 / R2) | 26 / 8 / 0 / 98 🟡 (upstream; its 98-cycle Read 2 is not tied to a version) | 28 / 8 / 0 / 91 | 28 / 10 / 10 / 90 |

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

CG000731 gives no TSO, cDNA primer, adaptor or index-primer sequences — only the bead
oligo and the final product (§5). 🔴 for the v4 TSO / cDNA primers / adaptor. Upstream
assumes they equal v3. 🟡

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

Conditions 🟢 from CG000183 (v3) / CG000731 (v4) as noted. The v2 guide (CG00052) was not
read; its steps are named in CG000108 (1.5 GEM-RT, 2.2 cDNA amplification, 3.3 adapter
ligation, 3.5 sample-index PCR) but conditions are 🔴 here.

1. **GEM generation.** Cells, RT master mix (with TSO), gel beads and partitioning oil
   meet on the chip; the bead dissolves, its oligos are released and the cell lyses. 🟢
2. **GEM-RT.** v3: 53 °C 45 min, 85 °C 5 min, 4 °C. v4: **48 °C** 45 min, 85 °C 5 min. 🟢
   The bead oligo primes at the poly(A) tail; at the mRNA 5' end the RT adds non-templated
   C's and switches to the TSO. First strand, 5'→3':
   `Partial Read 1 · <cbc16> · <umi> · T30VN · <cDNA> · CCCATGTACTCTGCGTTGATACCACTGCTT`
   🟢 (CG000108 "Reverse Transcript Product"; CG000183 likewise with 12-N UMI and `T30`).
   The final 30 nt are `CCC` + revcomp(TSO DNA part). 🟡 (computed) —
   [template switching](../ref/concepts/template-switching.md),
   [reverse transcription](../ref/concepts/reverse-transcription.md).
3. **Break GEMs, Dynabeads MyOne SILANE cleanup** of the single-stranded cDNA. 🟢
4. **cDNA amplification** with the forward (Partial Read 1) and reverse (TSO-handle)
   primers: 98 °C 3 min; cycles of 98 °C 15 s, 63 °C 20 s, 72 °C 1 min; 72 °C 1 min (v3).
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
v3/v3.1: P5 · TRUSEQ_READ1[4:] · <cbc 16> · <umi 12> · T30(VN) · <cDNA> · A · INDEX1_PRIMER · <i7' 8> · P7_RC
v4:      P5 · <i5 10> · TRUSEQ_READ1 · <cbc 16> · <umi 12> · T30VN · <cDNA> · A · INDEX1_PRIMER · <i7' 10> · P7_RC
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
- v4: CG000731 shows the 10-nt i5 **between P5 and the full 33-nt TruSeq Read 1**
  (`…GATCTACAC-N10-ACACTCTTTCCC…`), so the 4-nt P5/Read 1 overlap of v2/v3 is gone. 🟢
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
| Read 1 | TruSeq Read 1 (`illumina.TRUSEQ_READ1`) on the bottom strand | `<cbc16>` + `<umi>` | v2 26 🟡 (upstream); v3 28 🟢; v4 28 🟢 |
| i7 | `illumina.INDEX1_PRIMER` = the adapter top strand 🟡 | sample index | v2/v3 8 🟢 (v3); v4 10 🟢 |
| i5 | P5-side (`illumina.INDEX2_PRIMER_RC` or flow-cell primed, instrument dependent) | v4 only | 0 (v3) / 10 (v4) 🟢 |
| Read 2 | TruSeq Read 2 (`illumina.TRUSEQ_READ2`) | cDNA, **sense** to the mRNA, starting at the fragmentation end 🟡 | 98 🟡 (upstream, no version given; v2 guide not read); v3 91 🟢; v4 90 🟢 |

- v3: "Paired-end, single indexing", 28 / 8 / 0 / 91. v4: "Paired-end, dual indexing",
  28 / 10 / 10 / 90, minimum 20,000 read pairs per cell. 🟢
- Read 1 primer's 3' end = the bead handle, so Read 1 base 1 is the first barcode base;
  26 = 16 + 10 (v2), 28 = 16 + 12 (v3/v4). 🟡 (computed)
- `rc(TRUSEQ_READ2)` = `A` + the 33-nt adapter top strand: the Read 2 primer's 3' T sits on
  the A-tail, so Read 2 starts with the last insert base before the A-tail. 🟡 (computed)
- The 10x documents say the 16-bp barcode and the UMI "are encoded in Read 1, while Read 2
  is used to sequence the cDNA fragment" and the sample index is "the i7 index read"
  (v3) / "i7 and i5 … sample index reads" (v4). 🟢 (short quotes).

Barcode whitelists (upstream, 🟡): v2 `737K-august-2016.txt`, v3/v3.1
`3M-february-2018.txt`, v4 `3M-3pgex-may-2023.txt` (all from Cell Ranger).

## 7. Open questions

- 🔴 **v2 conditions**: GEM-RT temperature/time, cDNA cycles and fragmentation for v2
  (CG00052) were not read; §4 uses v3/v4 values.
- 🔴 **v3.1** is named by upstream but no v3.1 guide (CG000204 single index, CG000315 dual
  index) was read. A v3.1 dual-index library would presumably look like v4 (i5 + i7
  10 nt). 🟡 (guess)
- 🔴 **v4 TSO, cDNA primers and adaptor**: not printed in CG000731. Whether the v4 TSO is
  still the v3 sequence is not documented in what was read.
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
