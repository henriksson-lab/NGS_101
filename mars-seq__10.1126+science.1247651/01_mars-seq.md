# MARS-seq / MARS-seq2.0 — FACS into 384-well plates, T7 linear amplification, three barcodes

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only where we could not fetch it. Relationships marked 🟡
> *(computed)* were worked out with `lib/` while writing this note.
>
> **The scg_lib_structs page is inspiration, not evidence.** Anything whose only support
> is that page is 🔴 *(upstream-only)*: it is not cited as a source, no sequence is taken
> from it, and agreement with it is not verification. Nothing in `tools/marsseq.py` or
> `tools/page_parts.py` depends on an upstream-only sequence. Note also that a filename
> prefix is not provenance: `scg_jaitin-sm.pdf` contains *Science*'s own cover sheet and
> Tables S7–S9 (naming DOI 10.1126/science.1247651 and the 14 February 2014 issue),
> reached through the scg_lib_structs repository mirror. It is an original table
> extract, not the complete supplement.

**MARS-seq** — Jaitin DA, Kenigsberg E, Keren-Shaul H, Elefant N, Paul F, Zaretsky I,
Mildner A, Cohen N, Jung S, Tanay A, Amit I. "Massively parallel single-cell RNA-seq for
marker-free decomposition of tissues into cell types." *Science* 343:776–779 (2014).
doi:[10.1126/science.1247651](https://doi.org/10.1126/science.1247651) · PMID 24531970 ·
PMC4412462 (author manuscript).

**MARS-seq2.0** — Keren-Shaul H, Kenigsberg E, Jaitin DA, David E, Paul F, Tanay A, Amit I.
"MARS-seq2.0: an experimental and analytical pipeline for indexed sorting combined with
single-cell RNA sequencing." *Nature Protocols* 14:1841–1862 (2019).
doi:[10.1038/s41596-019-0164-4](https://doi.org/10.1038/s41596-019-0164-4) · PMID 31101904.
Data: GEO GSE123392.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/mars-seq-mars-seq2-0__10.1126+science.1247651/`, never committed):

| File | What | Used for |
|---|---|---|
| `science.1247651_PMC4412462.html.txt` | 2014 author manuscript, main text only | three-level barcoding, 1536 cells/lane, 384-well FACS |
| `scg_jaitin-sm.pdf.txt` | **Jaitin et al. 2014 *Science* Supplementary Material, two-page table extract** — the publisher's cover sheet and original Tables S7–S9, mirrored in scg_lib_structs (`MANIFEST.tsv`: `scg_lib_structs repo file data/MARS-seq/jaitin-sm.pdf`) | **Tables S7–S8**: MARS-seq RT1, the 8 ligation adapters, the second RT primer and both library-PCR primers. The file's cover sheet says it includes Methods and figures, but the local PDF itself has only two pages and does not contain them. |
| [Coauthor-uploaded full 2014 Supplementary Material](https://www.researchgate.net/publication/279842956_Jaitin-SM) | Alexander Mildner's ResearchGate upload of the original Jaitin-SM.pdf; the 5.02 MB PDF download itself returned HTTP 403, but the page displays its full original text | Original Materials and Methods, source page lines 124–390: capture, RT, pooling, Exo I, second strand, IVT, zinc fragmentation, T4 RNA ligase I, RT2 and PCR. The page text and the separate two-page extract agree on the supplement's identity and table content. |
| `s41596-019-0164-4_41596_2019_164_MOESM4_ESM.xlsx.txt` | MARS-seq2.0 Supplementary Table 1, sheet "MARS-seq2.0_RT1 primers" | **all 384 RT1 primers** |
| `s41596-019-0164-4_41596_2019_164_MOESM5_ESM.xlsx.txt` | MARS-seq2.0 Supplementary Table 2, sheet "MARS-seq2.0_Ligation adapter" | **all 32 ligation adapters**, their read-out barcodes |
| `s41596-019-0164-4_41596_2019_164_MOESM1_ESM.pdf.txt` | MARS-seq2.0 Supplementary Information (Supp. Figs 1–10, Manuals 1–3) | step names, QC, read-design syntax, FACS handling |
| `s41596-019-0164-4_41596_2019_164_MOESM2_ESM.pdf.txt` | Reporting summary | batch sizes (192–384 cells) |
| `s41596-019-0164-4_41596_2019_164_MOESM3_ESM.zip` | robot (Bravo `.vzp`) protocol files — listed, not opened | names of the liquid-handling steps |
| [Keren-Shaul et al. MARS-seq2.0 author procedure](https://doi.org/10.17504/protocols.io.7hkhj4w) | original 54-page protocols.io PDF, fetched into the MARS-seq2.0 source bundle | reagent setup p. 29 prints RT2 and both enrichment primers; resolves their sequence identity for 2.0 |
| `upstream_MARS-seq.html.txt` | scg_lib_structs MARS-seq page — **inspiration only, not a source** | orientation in the chemistry while reading the papers; no sequence and no claim on this page rests on it |

Could **not** be fetched (get by hand and drop in the same directory):

| What | URL | Why it matters |
|---|---|---|
| MARS-seq2.0 Nature Protocols main text | https://www.nature.com/articles/s41596-019-0164-4 (paywalled) | version-specific procedure and conditions not covered by the public author procedure; its RT2 and enrichment-primer sequences are available from protocols.io |

So: 🟢 covers the 2014 MARS-seq oligos and methods as well as the MARS-seq2.0 barcoded
oligos and QC/analysis supplement. The unavailable boundary is now the MARS-seq2.0 main
protocol, not the 2014 library structure.

---

## 1. What it is

Single cells are **FACS-sorted** (optionally index-sorted, recording surface-marker
intensities per well) into **384-well plates** with lysis solution and a well-specific
barcoded oligo-dT **RT1** primer. The 2014 Methods pool first-strand products into
24 intermediate wells, treat them with Exo I, and then pool 192-cell batches. The
pooled material undergoes second-strand synthesis and **T7 in-vitro transcription**
(linear amplification), RNA
fragmentation, **ligation of a pool-barcoded adapter** to the RNA 3' end, a second RT
(**RT2**), and PCR. 🟢 for the general flow (2014 main text; MARS-seq2.0 Supp. Fig. 5
legend, which lists RT, second-strand synthesis, IVT, shearing, PCR in that order).

It is a 3'-end counting method in the CEL-seq lineage: the T7 promoter rides on the
oligo-dT primer, so only barcoded first-strand products are amplified. 🟡

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Three levels of barcode**: molecule (UMI), cell/well (on RT1), plate/pool (on the ligated adapter) 🟢 (2014: "molecular, cellular, and plate level tags") | the "pool barcode added after pooling" idea returns in every split-pool and plate-pool method |
| 2 | Pool barcode put on by **single-strand ligation of a 3'-blocked DNA adapter to fragmented aRNA**, not by PCR index | adapter ligation to RNA: [small-RNA ligation](../ref/concepts/small-rna-ligation.md); CEL-seq instead uses Illumina RNA adapters + PCR index |
| 3 | Pool early: only RT1 happens per well; all later steps are pooled and automated | CEL-seq, Drop-seq-era plate methods |
| 4 | (2.0) **Index sorting** + robotics + empty-well background control, longer barcode/UMI, 32 pool barcodes | — |

## 2. Oligos

### 2.1 MARS-seq2.0 RT1 (well barcode) — 🟢

Supplementary Table 1 (`MOESM4`, sheet "MARS-seq2.0_RT1 primers"): 384 primers, columns
WellPosition / Name / Barcode / RT1 primer; four 96-well groups named `v3_gr1` … `v3_gr4`.
No modifications are written. First row, verbatim:

```
A1  v3_gr1  CTATTCG  CGATTGAGGCCGGTAATACGACTCACTATAGGGGCGACGTGTGCTCTTCCGATCTCTATTCGNNNNNNNNTTTTTTTTTTTTTTTTTTTTN
```

General form 🟡 (computed: all 384 rows match it exactly, 91 nt):

```
RT1  CGATTGAGGCCGGTAATACGACTCACTATAGGGGCGACGTGTGCTCTTCCGATCT <7-nt well barcode> NNNNNNNN TTTTTTTTTTTTTTTTTTTT N
```

### 2.2 MARS-seq2.0 ligation adapters (pool barcode) — 🟢

Supplementary Table 2 (`MOESM5`): 32 adapters, `/5Phos/` = 5' phosphate, `/3SpC3/` = 3'
C3 spacer (IDT notation, as given). Columns: name, adapter, "Seq barcode", "Read".

```
lig_N5X4_ix1   /5Phos/GACTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   AGTC  NNNNNAGTC
lig_N5X4_ix2   /5Phos/CATGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CATG  NNNNNCATG
lig_N5X4_ix3   /5Phos/CCAANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TTGG  NNNNNTTGG
lig_N5X4_ix4   /5Phos/CTGTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   ACAG  NNNNNACAG
lig_N5X4_ix5   /5Phos/GTAGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CTAC  NNNNNCTAC
lig_N5X4_ix6   /5Phos/TGATNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   ATCA  NNNNNATCA
lig_N5X4_ix7   /5Phos/ATCANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TGAT  NNNNNTGAT
lig_N5X4_ix8   /5Phos/TAGANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TCTA  NNNNNTCTA
lig_N5X4_ix9   /5Phos/AAGTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   ACTT  NNNNNACTT
lig_N5X4_ix10  /5Phos/GGCGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CGCC  NNNNNCGCC
lig_N5X4_ix11  /5Phos/GTTTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   AAAC  NNNNNAAAC
lig_N5X4_ix12  /5Phos/GCGCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GCGC  NNNNNGCGC
lig_N5X4_ix13  /5Phos/GAAANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TTTC  NNNNNTTTC
lig_N5X4_ix14  /5Phos/TACCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GGTA  NNNNNGGTA
lig_N5X4_ix15  /5Phos/CGGANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TCCG  NNNNNTCCG
lig_N5X4_ix16  /5Phos/CCCTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   AGGG  NNNNNAGGG
lig_N5X4_ix17  /5Phos/TCAGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CTGA  NNNNNCTGA
lig_N5X4_ix18  /5Phos/CTCGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CGAG  NNNNNCGAG
lig_N5X4_ix19  /5Phos/CTACNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GTAG  NNNNNGTAG
lig_N5X4_ix20  /5Phos/CTTANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TAAG  NNNNNTAAG
lig_N5X4_ix21  /5Phos/TGGCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GCCA  NNNNNGCCA
lig_N5X4_ix22  /5Phos/AGCTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   AGCT  NNNNNAGCT
lig_N5X4_ix23  /5Phos/CAGCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GCTG  NNNNNGCTG
lig_N5X4_ix24  /5Phos/ACTTNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   AAGT  NNNNNAAGT
lig_N5X4_ix25  /5Phos/TCTANNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   TAGA  NNNNNTAGA
lig_N5X4_ix26  /5Phos/ACCGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CGGT  NNNNNCGGT
lig_N5X4_ix27  /5Phos/ATGCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GCAT  NNNNNGCAT
lig_N5X4_ix28  /5Phos/GATCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GATC  NNNNNGATC
lig_N5X4_ix29  /5Phos/GGACNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GTCC  NNNNNGTCC
lig_N5X4_ix30  /5Phos/GTCCNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   GGAC  NNNNNGGAC
lig_N5X4_ix31  /5Phos/CGAGNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   CTCG  NNNNNCTCG
lig_N5X4_ix32  /5Phos/GCATNNNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/   ATGC  NNNNNATGC
```

### 2.3 MARS-seq (2014) RT, ligation and library-PCR primers — 🟢

The Jaitin 2014 *Science* Supplementary Material **is on disk** and prints all of these.
Verbatim from `scg_jaitin-sm.pdf.txt`, **Table S7 "Primers used during MARS-seq library
construction"** (lines 23–57; PDF line wraps rejoined, which is why the strings appear
split in the extracted text):

```
barcoded RT primer   CGATTGAGGCCGGTAATACGACTCACTATAGGGGCGACGTGTGCTCTTCCGATCTXXXXXXNNNNTTTTTTTTTTTTTTTTTTTTN
                     (XXXXXX = cell barcode, NNNN = RMT, i.e. the UMI — the table says so)
ligation adapter     AGATCGGAAGAGCGTCGTGTAG, modified with a phosphate group at 5' and a C3 spacer (blocker) at the 3'
Second RT primer     TCTAGCCTTCTCGCAGCACATC
P5_Rd1 PCR forward   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
P7_Rd2 PCR reverse   CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

and **Table S8 "Indexed ligation adapters for 384-well plates"**, all eight rows verbatim:

```
lig_NNNX4_ix1  /5Phos/GACTNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix2  /5Phos/CATGNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix3  /5Phos/CCAANNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix4  /5Phos/CTGTNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix5  /5Phos/GTAGNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix6  /5Phos/TGATNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix7  /5Phos/ATCANNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
lig_NNNX4_ix8  /5Phos/TAGANNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
```

So the 2014 library-PCR pair is **published, not upstream**:

| Primer | Status | Source |
|---|---|---|
| `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` | 🟢 | Jaitin et al. 2014, *Science* 343:776, Supplementary Material, Table S7, row "P5_Rd1 PCR forward" |
| `CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` | 🟢 | the same table, row "P7_Rd2 PCR reverse" |
| `TCTAGCCTTCTCGCAGCACATC` | 🟢 | the same table, row "Second RT primer" |

`tools/marsseq.py` holds the first two as `MARS1_P5_PCR` / `MARS1_P7_PCR` and the third as
`MARS1_RT2_PRINTED`, and the MARS-seq graph consumes them as exact `PCRGeometry` and RT2
inputs. Table S7 is also the source of the sequence `illumina`-equality relations in §3.

#### A version-specific RT2 found after checking the author protocol

```
CTACACGACGCTCTTCCGATCT        MARS-seq2.0 "Second RT primer" (author protocols.io procedure)
```

The 2014 Table S7 instead prints `TCTAGCCTTCTCGCAGCACATC`. The later authors'
[MARS-seq2.0 laboratory procedure](https://doi.org/10.17504/protocols.io.7hkhj4w),
reagent setup p. 29, independently prints `CTACACGACGCTCTTCCGATCT` as **its** RT2.
That sequence equals `revcomp(AGATCGGAAGAGCGTCGTGTAG)`, the reverse complement of the
ligation-adapter constant, and the 3′ 22 nt of `illumina.TRUSEQ_READ1` (🟡 computed,
§3). The 2.0 model holds the author-published RT2 and enrichment primers; the 2014
model retains its distinct printed RT2, whose landing conflicts with the 2014 Methods'
description of a primer complementary to the adapter.

## 3. How the oligos interlock — 🟡 (computed with `lib/`)

**RT1 (55-nt constant part + barcode + UMI + dT):**

- nt 1–13 `CGATTGAGGCCGG` — a leader 5' of the promoter (T7 RNA polymerase binds poorly at
  the very end of a duplex). Function not stated in the fetched sources. 🔴
- nt 14–33 **T7 promoter** `TAATACGACTCACTATAGGG` (the class III consensus with its +1 G),
  followed by one more `G`: the primer reads `…CACTATAGGGGCGACG…`. T7 starts at the first G
  of `GGGG`, so every aRNA begins `5'-GGGGCGACGUGUGCUCUUCCGAUCU…` (25 nt before the well
  barcode).
- nt 36–55 `GACGTGTGCTCTTCCGATCT` — the **3'-terminal 20 nt of `illumina.TRUSEQ_READ2`**
  (34 nt). The preceding `GC` differs from TruSeq (`…TTCA|GACG…`), so only these 20 nt are
  shared.
- then 7-nt well barcode, 8-nt UMI (`N8`), `T20`, and a 3' `N` (not `VN`).
- 384 well barcodes, all distinct, minimum pairwise Hamming distance **2**. 🟡

**P7_Rd2_PCR** (🟢 Table S7) **= `illumina.P7` + `illumina.TRUSEQ_READ2`** exactly (24 + 34 = 58 nt); its 3'
20 nt anneal to the RT1-derived handle, and the remaining 14 nt of Read 2 plus P7 are added
by the primer. No i7 index in it: the library has no sample index read. 🟡

**Ligation adapter = `/5Phos/` + 4-nt pool barcode + N5 + `AGATCGGAAGAGCGTCGTGTAG` +
`/3SpC3/`** (31 nt, 2.0). The 22-nt constant part is the first 22 nt of
`illumina.INDEX2_PRIMER_RC` (= the TruSeq Read 1 site, as it appears on the bottom
strand). Its reverse complement is `CTACACGACGCTCTTCCGATCT`, which is the 3' 22 nt of
`illumina.TRUSEQ_READ1`. 🟡 (computed) So the adapter *presents* a TruSeq Read 1 site that
a primer of that sequence would bind — but **the authors' Second RT primer is
`TCTAGCCTTCTCGCAGCACATC`** (🟢 Table S7), which is not it. Where RT2 binds in the 2014
chemistry is therefore not established by this identity; the sequence
`CTACACGACGCTCTTCCGATCT` is an upstream claim about MARS-seq2.0 only (🔴, §2.3).

- The 5' phosphate makes the adapter the **donor** for ligation onto the 3'-OH of an aRNA
  fragment; the 3' C3 spacer stops it acting as an acceptor, so adapters cannot
  concatenate (the logic in [small-RNA ligation](../ref/concepts/small-rna-ligation.md),
  trick 1, without pre-adenylation). 🟡
- **P5_Rd1_PCR** (🟢 Table S7) **= `illumina.P5` + `TRUSEQ_READ1[4:]`** (29 + 29 = 58 nt;
  it is `illumina.NEBNEXT_UNIVERSAL_PRIMER`). 🟡 (computed) Its 3' 22 nt equal
  `revcomp(AGATCGGAAGAGCGTCGTGTAG)`, the ligation-adapter constant's complement, so it
  binds the adapter end of the RT2 product and carries P5 onto it.
- **The "Seq barcode" column is the reverse complement of the adapter's 4 nt** for all 32
  (GACT → AGTC …; ix2, 12, 22, 28 are palindromes), and "Read" = `NNNNN` + that. This is
  what Read 1 shows: the adapter is read through as its complement, so Read 1 begins with
  the 5 random bases, then the pool barcode, then cDNA. 🟡
- 32 pool barcodes, all distinct, minimum Hamming distance **2**. The 8 MARS-seq1
  barcodes `lig_NNNX4_ix1..8` (upstream) are **identical to `lig_N5X4_ix1..8`**; only the
  random part grew from N3 to N5. 🟡
- Capacity: 384 well barcodes × 32 pool barcodes = **12,288 cells** distinguishable in one
  sequencing pool. 🟡

## 4. Step by step

Reaction order from the [coauthor-uploaded 2014 Supplementary Methods](https://www.researchgate.net/publication/279842956_Jaitin-SM),
cross-checked against the original 2014 article and the two-page primer-table extract.
The full Methods are readable on that page even though its PDF download is blocked.
The later MARS-seq2.0 notes below are identified separately.

1. **Capture plates**: 384-well plates containing 2 µL hypotonic lysis solution, RNase
   inhibitor and one barcoded RT1 primer per well 🟢 (2014 Methods, source page lines
   124–172). The later 2.0 protocol's capture-plate robot settings are documented in
   its own supplement.
2. **Sort** one cell per well, leaving **two wells empty** per 384-well plate as no-cell
   controls; spin, snap-freeze and store at −80 °C 🟢 (2014 Methods, source page lines
   108–122). The later 2.0 protocol uses a different set of control wells.
3. **RT1** in each well: RT1 anneals to poly(A) and is extended into first-strand cDNA,
   writing well barcode + UMI onto every molecule. Superscript III RT, 2 min at 42 °C,
   50 min at 50 °C and 5 min at 85 °C 🟢 (2014 Methods, source page lines 174–252).
   See [reverse transcription](../ref/concepts/reverse-transcription.md).
4. **Pool** barcoded first-strand cDNA from the wells into 24 wells 🟢 (2014 Methods,
   source page lines 254–267).
5. **Exonuclease I** after the first pooling to remove RT primer leftovers, then SPRI
   cleanup and final pooling 🟢 (2014 Methods, source page lines 269–305).
6. **Second strand** on pooled cDNA with the NEB second-strand synthesis kit, 2.5 h at
   16 °C 🟢 (2014 Methods, source page lines 311–315). The original Methods name the kit,
   not its individual enzyme composition; RNase H and DNA polymerase I remain unverified.
   The double-stranded T7 promoter is now complete.
7. **IVT** with the NEB T7 High Yield kit overnight at 37 °C, followed by Turbo DNase I
   template removal 🟢 (2014 Methods, source page lines 317–324). The linear amplification yields antisense aRNA
   `5'-GGGGCGACG…-<well bc>-<UMI>-U20-<antisense cDNA>-3'` 🟡 (computed from RT1; upstream
   draws it starting at `GG`, see §7).
8. **Fragment** the aRNA in Zn²⁺ fragmentation solution for 3 min at 70 °C 🟢 (2014
   Methods, source page lines 326–331). Only the **5' fragment of each aRNA** carries the well barcode / UMI
   handle, so only it can be amplified by the P7_Rd2 primer later. 🟡
9. **Ligate** a `lig_N5X4_ix<k>` adapter (one per pool) to fragmented RNA with **T4 RNA
   ligase I**, 2 h at 22 °C 🟢 (2014 Methods, source page lines 333–368). Pool barcode
   now sits on the other end of the insert.
10. **RT2** on the ligated adapter with AffinityScript RT and the published **Second RT primer**
    `TCTAGCCTTCTCGCAGCACATC` 🟢 (Table S7); the product is
    `5'-<RT2 end>-<N3'>-<pool bc'>-<sense cDNA>-A20-<UMI'>-<well bc'>-AGATCGGAAGAGCACACGTCGCCCC-3'`
    🟡 (computed; the 3' end is the complement of the aRNA 5' end, so it ends in `CCCC`).
    🟢 Enzyme and timing from 2014 Methods, source page lines 370–377. 🔴 Which part of the ligated adapter that primer anneals to is not stated, and it is
    not the adapter's own reverse complement (§3). MARS-seq2.0 uses the distinct RT2
    published in its later author protocols.io procedure (§2.3).
11. **PCR** with P5_Rd1_PCR `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`
    and P7_Rd2_PCR `CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` 🟢
    (both Table S7, verbatim); qPCR QC at three points —
    QC1 pooled cDNA, QC2 after RT2, QC3 library Ct, usually on "Actb 3'" 🟢 (Supp. Manual 2
    field names). Typical final library average 384 bp on TapeStation 🟢 (Supp. Fig. 10).

## 5. Final library — 🟡 (assembled from the oligos above)

MARS-seq2.0, top strand (the strand whose sequence Read 1 reports), 5'→3':

```
5'- P5 (29) · TruSeq Read1 site minus its first 4 nt ACAC, which P5 already ends in (29) · N5 · <pool bc, 4, as "Seq barcode"> · <sense cDNA insert> · A(n) · <UMI', 8> · <well bc', 7> · AGATCGGAAGAGCACACGTC (rc of RT1 nt 36-55 = rc partial Read2) · TGAACTCCAGTCAC · P7' (24) -3'
```

where the last three segments together are `revcomp(P7_Rd2_PCR)` = `TRUSEQ_READ2` reverse
complement + `illumina.P7_RC` (computed from the Table S7 primer).
Note the leader and T7 promoter of RT1 are **not** in the library (IVT starts at the +1 G),
and neither is the aRNA's first 5 nt `GGGGC` (RT1 nt 31–35): the RT2 cDNA ends in
`…CACGTCGCCCC-3'`, but the P7_Rd2 primer anneals only over RT1 nt 36–55 and its 5' part
replaces the `GCCCC` in the copies made from it. 🟡

MARS-seq (2014): same, with N3 instead of N5, a 6-nt well barcode and 4-nt UMI 🟢
(Supplementary Tables S7–S8).

## 6. Read layout / sequencing

- Read 1 (TruSeq Read 1 primer): N5 (2.0) or N3 (MARS-seq1), pool barcode (4), then cDNA
  (sense strand, near the 3' end of the transcript). MARS-seq1 is 🟢 from Tables S7–S8;
  the 2.0 placement is 🟡 derived from its published adapter table.
- Read 2 (TruSeq Read 2 primer): MARS-seq1 reports the 6-nt well barcode, 4-nt UMI, then
  poly(T) 🟢; MARS-seq2.0 analogously uses 7 + 8 nt 🟡.
- No index reads are needed: demultiplexing is by inline pool and well barcodes. 🟢 for
  the 2014 construct.
- The pipeline's read-design syntax 🟢 (Supp. Manual 2): fields `R1_design`, `I5_design`,
  `R2_design`, codes P = pool barcode, M = mRNA, W = well barcode, R = UMI, I = ignore. The
  worked example is **`R1_design = 3I.4P.50M`, `R2_design = 7W.8R.5I`** 🟢 (Supplementary
  Information, read-design fields, verbatim). Read 2 matches the 2.0 RT1 exactly (7 + 8,
  then 5 nt of dT ignored), while Read 1 `3I.4P` matches the **MARS-seq1** N3 adapter.
  🟡 (derived) With the 2.0 N5 adapter Read 1 would instead be `5I.4P` (as the `MOESM5`
  "Read" column `NNNNNAGTC` implies) — that `5I` is derived from the v2 adapter and is
  **not** a published read design. The schematic's read caption says which is which.
- 2014: 1536 cells per sequencing lane, ~22,000 aligned reads per cell 🟢 (main text).
  Instrument and cycle numbers 🔴.

## 7. What the upstream page suggested, and what the authors actually say

The scg_lib_structs page is inspiration only. This table records what it prompted us to
look for and what the **author sources** then established; the right-hand column is the
evidence, and agreement with upstream is never itself a verdict. Where the papers are
silent, the row stays 🔴 and the model asserts nothing.

| Upstream prompted a look at | Author source checked | What the source establishes |
|---|---|---|
| MARS-seq2.0 RT1 = constant + 7-bp cell barcode + 8-bp UMI + T20 + N | all 384 rows of `MOESM4` | **agrees** 🟢 |
| 32 adapters `lig_N5X4_ix1..32`, sequences listed | diffed against `MOESM5` | **identical** 🟢 |
| MARS-seq2.0 plate barcode "9 bp" | `MOESM5` "Read" column: N5 + 4-nt barcode | agrees in length; only 4 nt are the barcode, 5 are random 🟢 |
| MARS-seq1 UMI 4 bp and cell barcode 6 bp | 2014 Supplementary Table S7 | **agrees** 🟢 |
| the second RT primer | Jaitin 2014 Table S7; Keren-Shaul et al. 2019 author protocols.io reagent setup p. 29 | 2014 prints `TCTAGCCTTCTCGCAGCACATC`; 2.0 prints the distinct `CTACACGACGCTCTTCCGATCT` 🟢. The latter is the adapter reverse complement; the 2014 printed oligo is not. Both models use their respective author-published sequence. |
| aRNA drawn as `5'-GGCGACGUGUG…` | T7 +1 is the first G of `GGGG` in RT1 | upstream drops two G's; immaterial to the library 🟡 |
| "the 5'-/acrydite/iSpPC/ is omitted for simplicity" | `MOESM4` RT1 has no modifications | appears to be a leftover from another page (inDrop); MARS-seq RT1 is unmodified as listed 🟡 |
| "not sure what NNN between Partial Rd1 and plate barcode is" | `MOESM5` | in 2.0 it is N5 and the pipeline treats it as ignore (I); a UMI-like role is not stated 🔴 |
| the library-PCR pair | Jaitin 2014 Table S7, rows "P5_Rd1 PCR forward" and "P7_Rd2 PCR reverse" | both printed verbatim 🟢 — `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` and `CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`. They are the page's evidence; the identical strings on the upstream page are irrelevant to it |
| a MARS-seq2.0-specific RT2 and PCR pair | Keren-Shaul et al. 2019 author protocols.io procedure, reagent setup p. 29 | RT2 and both P5-Read1/P7-Read2 enrichment primers are printed verbatim 🟢 and used by the 2.0 page. |

## 8. Open questions

- 🔴 Whether the inaccessible *Nature Protocols* main text differs from the authors'
  public detailed protocols.io procedure in any version-specific condition or correction.
- 🔴 What the 13-nt leader `CGATTGAGGCCGG` is for, and why the 3' anchor is a single `N`.
- 🔴 Purpose of the N5 / N3 random bases next to the pool barcode (diversity for the first
  sequencing cycles? an extra UMI?). The pipeline ignores them.
- 🔴 Why the 2014 Second RT primer differs from the 2.0 primer, and where the 2014 primer
  `TCTAGCCTTCTCGCAGCACATC` anneals on the ligated adapter — the supplement prints the
  oligo, while its Methods describe a primer complementary to the adapter.
- 🟡 "Template switching errors" are a QC category (Supp. Fig. 6G) although the chemistry
  has no TSO: presumably RT1-primed chimeras; not defined in the fetched text.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` got the PMC main text, all five Nature Protocols supplements and
— through the `scg:` handler mirroring `data/MARS-seq/jaitin-sm.pdf` — **a two-page extract
of the 2014 *Science* Supplementary Material containing Tables S7–S9**. The 2014 Methods
and figure pages are absent from that extract. Coauthor Alexander Mildner uploaded the
full original supplement; its PDF download is blocked here, but its complete text is
readable on ResearchGate and supplies the 2014 reaction order. The separate PMC figure
PDF remains gated, as does the *Nature Protocols* main text.
`tools/scrape_primers.py` found the RT1 and adapter tables; it labelled each RT1 hit with
its barcode rather than a name, and reported no `lib/` sequence inside RT1 (the 20-nt
partial TruSeq Read 2 is below its match length). The reaction order came from the
coauthor-uploaded original 2014 Methods. Reading them exposed an ordering error in
the shared graph: the paper pools first-strand cDNA before second-strand synthesis.

One fetch caveat to know about: the bundle was written under the scg_lib_structs slug
`mars-seq-mars-seq2-0__10.1126+science.1247651`, which combines the two protocols, so it
is invisible under this directory's own slug unless
`$CHEM_DATA/sources/mars-seq__10.1126+science.1247651` points at it. Tables S7–S8 were
read there, at `scg_jaitin-sm.pdf.txt` lines 23–57.

The final `--find` check could not confirm the oligos that contain `N`: a query `N`
matches only A/C/G/T, never a literal `N` in the table, so every RT1 and adapter query
returned 0 locations (while `GACTNNNNNAGATCG` "matched" unrelated RT1 rows). The 🟢
sequences were therefore confirmed by exact text match instead: the RT1 A1 row and all 32
adapter rows of §2 are byte-identical to `MOESM4` / `MOESM5`, and their N-free parts are
found by `--find`.
