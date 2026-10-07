# MARS-seq / MARS-seq2.0 — FACS into 384-well plates, T7 linear amplification, three barcodes

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only where we could not fetch it. Relationships marked 🟡
> *(computed)* were worked out with `lib/` while writing this note; they are not yet
> asserted in a self-test, because this protocol has no `tools/` module yet. Claims taken
> only from the upstream scg_lib_structs page are 🟡 *(upstream)*: a careful secondary
> source, but not the authors.

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
| `s41596-019-0164-4_41596_2019_164_MOESM4_ESM.xlsx.txt` | MARS-seq2.0 Supplementary Table 1, sheet "MARS-seq2.0_RT1 primers" | **all 384 RT1 primers** |
| `s41596-019-0164-4_41596_2019_164_MOESM5_ESM.xlsx.txt` | MARS-seq2.0 Supplementary Table 2, sheet "MARS-seq2.0_Ligation adapter" | **all 32 ligation adapters**, their read-out barcodes |
| `s41596-019-0164-4_41596_2019_164_MOESM1_ESM.pdf.txt` | MARS-seq2.0 Supplementary Information (Supp. Figs 1–10, Manuals 1–3) | step names, QC, read-design syntax, FACS handling |
| `s41596-019-0164-4_41596_2019_164_MOESM2_ESM.pdf.txt` | Reporting summary | batch sizes (192–384 cells) |
| `s41596-019-0164-4_41596_2019_164_MOESM3_ESM.zip` | robot (Bravo `.vzp`) protocol files — listed, not opened | names of the liquid-handling steps |
| `upstream_MARS-seq.html.txt` | scg_lib_structs MARS-seq page | full oligo list incl. PCR / RT2 primers, step drawings, read layout |

Could **not** be fetched (get by hand and drop in the same directory):

| What | URL | Why it matters |
|---|---|---|
| 2014 supplementary materials (methods, Tables S7/S8 = MARS-seq1 oligos) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4412462/bin/NIHMS666898-supplement-Sup_Figures.pdf (PMC download gate; Europe PMC and science.org both 403) | the only author source for MARS-seq1 RT1, RT2, PCR primers and every reaction condition |
| MARS-seq2.0 main text (procedure, reagents, Box/Table of primers) | https://www.nature.com/articles/s41596-019-0164-4 (paywalled) | every volume, temperature and enzyme; RT2 and PCR primer sequences; sequencing run parameters |

So: 🟢 covers the MARS-seq2.0 barcoded oligos and the QC/analysis supplement; the reaction
order and all chemistry conditions are 🟡 (upstream, plus Supp. Fig. legends) or 🔴.

---

## 1. What it is

Single cells are **FACS-sorted** (optionally index-sorted, recording surface-marker
intensities per well) into **384-well plates** that already hold a well-specific barcoded
oligo-dT **RT1** primer (and, presumably, lysis buffer 🟡: not stated in the fetched text). After first-strand synthesis everything
from one plate (or half-plate amplification batch) is **pooled** and processed as one
tube: second strand, **T7 in-vitro transcription** (linear amplification), RNA
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

### 2.3 Primers not in the fetched author sources — 🟡 (upstream)

From the scg_lib_structs page only; the MARS-seq2.0 main text and the 2014 supplement
(where these are given) could not be fetched.

```
2nd RT primer (RT2)   CTACACGACGCTCTTCCGATCT
P5_Rd1_PCR primer     AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
P7_Rd2_PCR primer     CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

MARS-seq (2014) per upstream (Tables S7/S8 of the 2014 paper, unchecked by us 🔴):

```
1st RT primer (MARS-seq)  CGATTGAGGCCGGTAATACGACTCACTATAGGGGCGACGTGTGCTCTTCCGATCT [6-bp cell barcode][4-bp UMI] TTTTTTTTTTTTTTTTTTTTN
lig_NNNX4_ix1..ix8        /5Phos/ <4-nt> NNN AGATCGGAAGAGCGTCGTGTAG /3SpC3/     4-nt = GACT CATG CCAA CTGT GTAG TGAT ATCA TAGA
```

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

**P7_Rd2_PCR = `illumina.P7` + `illumina.TRUSEQ_READ2`** exactly (24 + 34 = 58 nt); its 3'
20 nt anneal to the RT1-derived handle, and the remaining 14 nt of Read 2 plus P7 are added
by the primer. No i7 index in it: the library has no sample index read. 🟡

**Ligation adapter = `/5Phos/` + 4-nt pool barcode + N5 + `AGATCGGAAGAGCGTCGTGTAG` +
`/3SpC3/`** (31 nt, 2.0). The 22-nt constant part is the first 22 nt of
`illumina.INDEX2_PRIMER_RC` (= the TruSeq Read 1 site, as it appears on the bottom
strand). Its reverse complement is `CTACACGACGCTCTTCCGATCT` = **RT2**, which is the 3'
22 nt of `illumina.TRUSEQ_READ1`. So RT2 primes on the ligated adapter. 🟡

- The 5' phosphate makes the adapter the **donor** for ligation onto the 3'-OH of an aRNA
  fragment; the 3' C3 spacer stops it acting as an acceptor, so adapters cannot
  concatenate (the logic in [small-RNA ligation](../ref/concepts/small-rna-ligation.md),
  trick 1, without pre-adenylation). 🟡
- **P5_Rd1_PCR = `illumina.P5` + `TRUSEQ_READ1[4:]`** (29 + 29 = 58 nt; it is
  `illumina.NEBNEXT_UNIVERSAL_PRIMER`). Its 3' 22 nt equal RT2, so it extends the RT2 end.
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

Reaction order from the upstream drawing, cross-checked against the step names in the
MARS-seq2.0 supplement; conditions are 🔴 (main text not fetched) unless stated.

1. **Capture plates**: 384-well plates pre-loaded with one barcoded RT1 primer per well,
   stored at −80 °C, thawed on ice and spun down before unsealing 🟢 (Supp. Manual 1); that
   the wells also hold lysis buffer is not stated in the fetched sources 🟡 (inferred). **ERCC spike-in**
   (mix 1) added per well — pipeline defaults 1.25e-05 dilution, 0.01 µL per well 🟢 (Supp.
   Manual 2, amp_batches fields). Robot file "Cell capture plates preparation" 🟢 (zip).
2. **Sort** one cell per well (100 µm nozzle, 384-well calibration every 5 plates, optional
   index sort); keep **four wells empty** as background controls (e.g. O1, O2, P1, P2);
   seal, freeze on dry ice, store −80 °C 🟢 (Supp. Manual 1).
3. **RT1** in each well: RT1 anneals to poly(A) and is extended into first-strand cDNA,
   writing well barcode + UMI onto every molecule. Robot file "RT mix addition" 🟢; a
   "low-volume" RT was one of the 2.0 optimisations that lowered background (Supp. Fig. 2A
   legend) 🟢. See [reverse transcription](../ref/concepts/reverse-transcription.md).
4. **Exonuclease** to remove unused RT1 primer before pooling 🟡 (robot file "Exonuclease
   mix addition" 🟢; the purpose is inferred: free barcoded primers carried into the pool
   would mis-barcode other cells' RNA).
5. **Pool** the wells of an amplification batch (192–384 cells 🟢, reporting summary), clean
   up. Robot file "Pooling" 🟢.
6. **Second strand**: RNase H + DNA polymerase I (Gubler–Hoffman style) 🟡 (upstream);
   2.0 compared second-strand enzyme mixes for background: composition "A", "A" with
   RNase H added, "A" without DNA polymerase, and composition "B" diluted 1:4 or 1:8
   (Supp. Fig. 2B legend 🟢). The double-stranded
   T7 promoter is now complete.
7. **IVT** with T7 RNA polymerase: linear amplification, antisense aRNA
   `5'-GGGGCGACG…-<well bc>-<UMI>-U20-<antisense cDNA>-3'` 🟡 (computed from RT1; upstream
   draws it starting at `GG`, see §7).
8. **Fragment** the aRNA (heat, per upstream 🟡; "sheared at random positions" 🟢 Supp.
   Fig. 5 legend). Only the **5' fragment of each aRNA** carries the well barcode / UMI
   handle, so only it can be amplified by the P7_Rd2 primer later. 🟡
9. **Ligate** a `lig_N5X4_ix<k>` adapter (one per pool) to the 3'-OH of the fragments with
   **T4 RNA ligase 1** 🟡 (upstream). Pool barcode now on the other end of the insert.
10. **RT2** with `CTACACGACGCTCTTCCGATCT` on the adapter: cDNA
    `5'-<RT2>-<N5'>-<pool bc'>-<sense cDNA>-A20-<UMI'>-<well bc'>-AGATCGGAAGAGCACACGTCGCCCC-3'`
    🟡 (computed; the 3' end is the complement of the aRNA 5' end, so it ends in `CCCC`).
11. **PCR** with P5_Rd1_PCR and P7_Rd2_PCR 🟡 (upstream); qPCR QC at three points —
    QC1 pooled cDNA, QC2 after RT2, QC3 library Ct, usually on "Actb 3'" 🟢 (Supp. Manual 2
    field names). Typical final library average 384 bp on TapeStation 🟢 (Supp. Fig. 10).

## 5. Final library — 🟡 (assembled from the oligos above)

MARS-seq2.0, top strand (the strand whose sequence Read 1 reports), 5'→3':

```
5'- P5 (29) · TruSeq Read1 site minus its first 4 nt ACAC, which P5 already ends in (29) · N5 · <pool bc, 4, as "Seq barcode"> · <sense cDNA insert> · A(n) · <UMI', 8> · <well bc', 7> · AGATCGGAAGAGCACACGTC (rc of RT1 nt 36-55 = rc partial Read2) · TGAACTCCAGTCAC · P7' (24) -3'
```

where the last three segments together are `revcomp(P7_Rd2_PCR)` = `TRUSEQ_READ2` reverse
complement + `illumina.P7_RC` (computed: equal to the upstream drawing's top-strand tail).
Note the leader and T7 promoter of RT1 are **not** in the library (IVT starts at the +1 G),
and neither is the aRNA's first 5 nt `GGGGC` (RT1 nt 31–35): the RT2 cDNA ends in
`…CACGTCGCCCC-3'`, but the P7_Rd2 primer anneals only over RT1 nt 36–55 and its 5' part
replaces the `GCCCC` in the copies made from it. 🟡

MARS-seq (2014): same, with N3 instead of N5, a 6-nt well barcode and 4-nt UMI 🟡
(upstream; 🔴 against the paper).

## 6. Read layout / sequencing

- Read 1 (TruSeq Read 1 primer): N5 (2.0) or N3 (MARS-seq1), pool barcode (4), then cDNA
  (sense strand, near the 3' end of the transcript). 🟡
- Read 2 (TruSeq Read 2 primer): well barcode (7), UMI (8), then poly(T). 🟡
- No index reads needed: demultiplexing is by pool barcode in Read 1 and well barcode in
  Read 2. 🟡
- The pipeline's read-design syntax 🟢 (Supp. Manual 2): fields `R1_design`, `I5_design`,
  `R2_design`, codes P = pool barcode, M = mRNA, W = well barcode, R = UMI, I = ignore. The
  worked example is `R1_design = 3I.4P.50M`, `R2_design = 7W.8R.5I`: Read 2 matches the
  2.0 RT1 exactly (7 + 8, then 5 nt of dT ignored), while Read 1 `3I.4P` matches the
  **MARS-seq1** N3 adapter; with the 2.0 N5 adapter it would be `5I.4P` (as the "Read"
  column `NNNNNAGTC` implies). 🟡
- 2014: 1536 cells per sequencing lane, ~22,000 aligned reads per cell 🟢 (main text).
  Instrument and cycle numbers 🔴.

## 7. Upstream page vs. author sources

| Upstream says | Check | Verdict |
|---|---|---|
| MARS-seq2.0 RT1 = constant + 7-bp cell barcode + 8-bp UMI + T20 + N | all 384 rows of `MOESM4` | **agrees** 🟢 |
| 32 adapters `lig_N5X4_ix1..32`, sequences listed | diffed against `MOESM5` | **identical** 🟢 |
| MARS-seq2.0 plate barcode "9 bp" | `MOESM5` "Read" column: N5 + 4-nt barcode | agrees in length; only 4 nt are the barcode, 5 are random 🟢 |
| MARS-seq1 UMI 4 bp vs "4–8 bp" in 2014 methods; plate bc 7 bp vs "6 bp" | 2014 supplement not fetched | unverified 🔴 |
| "2nd RT primer" in Table S7 printed 3'→5' | not fetched | unverified 🔴; the RT2 sequence upstream uses is the correct orientation (= revcomp of adapter constant) 🟡 |
| aRNA drawn as `5'-GGCGACGUGUG…` | T7 +1 is the first G of `GGGG` in RT1 | upstream drops two G's; immaterial to the library 🟡 |
| "the 5'-/acrydite/iSpPC/ is omitted for simplicity" | `MOESM4` RT1 has no modifications | appears to be a leftover from another page (inDrop); MARS-seq RT1 is unmodified as listed 🟡 |
| "not sure what NNN between Partial Rd1 and plate barcode is" | `MOESM5` | in 2.0 it is N5 and the pipeline treats it as ignore (I); a UMI-like role is not stated 🔴 |
| P5_Rd1 / P7_Rd2 / RT2 sequences | not in fetched author files | 🟡 upstream only |

## 8. Open questions

- 🔴 All conditions: lysis buffer, RT enzyme and volume, exonuclease, second-strand mix,
  IVT time, fragmentation buffer/time, ligation, RT2, PCR cycles — in the paywalled
  MARS-seq2.0 main text.
- 🔴 What the 13-nt leader `CGATTGAGGCCGG` is for, and why the 3' anchor is a single `N`.
- 🔴 Purpose of the N5 / N3 random bases next to the pool barcode (diversity for the first
  sequencing cycles? an extra UMI?). The pipeline ignores them.
- 🔴 Whether RT2 / PCR primers in MARS-seq2.0 are exactly the upstream ones.
- 🟡 "Template switching errors" are a QC category (Supp. Fig. 6G) although the chemistry
  has no TSO: presumably RT1-primed chimeras; not defined in the fetched text.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` got the PMC main text, the upstream page and all five Nature
Protocols supplements, but not the 2014 supplement (PMC download gate) or the Nature
Protocols main text (paywall). `tools/scrape_primers.py` found the RT1 and adapter tables;
it labelled each RT1 hit with its barcode rather than a name, and reported no `lib/`
sequence inside RT1 (the 20-nt partial TruSeq Read 2 is below its match length). The
reaction order came from the upstream drawing and the supplementary figure legends.

The final `--find` check could not confirm the oligos that contain `N`: a query `N`
matches only A/C/G/T, never a literal `N` in the table, so every RT1 and adapter query
returned 0 locations (while `GACTNNNNNAGATCG` "matched" unrelated RT1 rows). The 🟢
sequences were therefore confirmed by exact text match instead: the RT1 A1 row and all 32
adapter rows of §2 are byte-identical to `MOESM4` / `MOESM5`, and their N-free parts are
found by `--find`.
