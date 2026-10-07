# STRT-seq family — barcoded 5'-end single-cell RNA-seq (STRT, STRT/C1, STRT-seq-2i)

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available to us. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note; they are not yet asserted in a self-test, because
> this protocol has no `tools/` module yet (status `notes`). Claims that only the upstream
> scg_lib_structs page makes are 🟡 *(upstream only)* — a secondary source.

Three generations from the Linnarsson lab, all **5'-end counting by template switching**,
all reading the transcript's 5' end in Read 1 right after the TSO's `GGG`:

| Variant | Paper | Cell address | UMI | Library chemistry |
|---|---|---|---|---|
| **STRT** | Islam S, Kjällquist U, Moliner A, Zajac P, Fan J-B, Lönnerberg P, Linnarsson S. "Characterization of the single-cell transcriptional landscape by highly multiplex RNA-seq." *Genome Res* 21:1160–1167 (2011). doi:[10.1101/gr.110882.110](https://doi.org/10.1101/gr.110882.110), PMC3129258 | 6-nt barcode in the **TSO**, 96 wells pooled after RT | none | DNase I/Mn²⁺ fragmentation, end repair, A-tail, P2 adapter ligation, SalI, PCR |
| **STRT/C1** ("STRT-seq-C1") | Islam S, Zeisel A, Joost S, La Manno G, Zajac P, Kasper M, Lönnerberg P, Linnarsson S. "Quantitative single-cell RNA-seq with unique molecular identifiers." *Nat Methods* 11:163–166 (2014). doi:[10.1038/nmeth.2772](https://doi.org/10.1038/nmeth.2772) | Fluidigm C1 chamber; 8-nt barcode in the **Tn5 adapter** | 5 nt in the TSO | barcoded Tn5 on amplified cDNA, no library PCR |
| **STRT-seq-2i** | Hochgerner H, Lönnerberg P, Hodge R, Mikes J, Heskol A, Hubschle H, Lin P, Picelli S, La Manno G, Ratz M, Dunne J, Husain S, Lein E, Srinivasan M, Zeisel A, Linnarsson S. *bioRxiv* doi:[10.1101/126268](https://doi.org/10.1101/126268) (v1, 20 Apr 2017, CC-BY 4.0); published as *Sci Rep* 7:16327 (2017), doi:[10.1038/s41598-017-16546-4](https://doi.org/10.1038/s41598-017-16546-4), PMC5703850 | 9600-well array: 5-nt **well** index in PCR primer × 8-nt **subarray** index in Tn5 adapter | 6 nt in the TSO | barcoded Tn5, streptavidin, ss elution, 8-cycle PCR |

Also cited by upstream: Islam S *et al.* "Highly multiplexed and strand-specific
single-cell RNA 5′ end sequencing." *Nat Protoc* 7:813–828 (2012),
doi:10.1038/nprot.2012.022 — the detailed STRT protocol that upstream's STRT drawing is
based on (not available to us); and nanoCAGE (Plessy *et al.* *Nat Methods* 2010,
doi:10.1038/nmeth.1470) for "semi-suppressive PCR". Data: GEO GSE29087 (STRT).

Sources read (in `$CHEM_DATA/sources/strt-seq-family__10.1101+gr.110882.110/`, never
committed):

| File | What | Used for |
|---|---|---|
| `gr.110882.110_PMC3129258.html(.txt)` | STRT, Genome Res 2011, full text from PMC (fetched by hand) | all STRT oligos and conditions |
| `nmeth.2772_..._MOESM266_ESM.pdf` | STRT/C1 Supplementary Figures + Note | list of changes vs STRT, sequencing |
| `nmeth.2772_..._MOESM268_ESM.xlsx` | STRT/C1 oligo table | **every C1 oligo, 96 Tn5 barcodes** |
| `nmeth.2772_..._MOESM267_ESM.xlsx` | gene CV table | not chemistry |
| `126268_v1.full.pdf` | STRT-seq-2i preprint | methods, "Primer sequences" table |
| `s41598-017-16546-4_PMC5703850.html(.txt)` | STRT-seq-2i, Sci Rep 2017 (fetched by hand) | Table 1 (oligos), methods — checked against the preprint |
| `s41598-017-16546-4_MOESM1_ESM.pdf` | Sci Rep Supplementary Information | Fig. S1/S2 legends (read length, optimisations) |
| `upstream_STRT-seq_family.html(.txt)` | scg_lib_structs page | second source, checked below |

Could not be fetched (🔴; list in `MANIFEST.tsv` as `(manual)`):

- STRT Supplemental Table 1 — the 96 **STRT-V2-n barcodes** — https://genome.cshlp.org/content/21/7/1160/suppl/DC1 (cshlp answers 403/404 to scripts; PMC has no copy).
- STRT/C1 **Online Methods** (reaction conditions, read-1 primer, PvuI step if any) — https://doi.org/10.1038/nmeth.2772 (paywalled, no PMC copy).
- STRT Nature Protocols 2012 — https://doi.org/10.1038/nprot.2012.022 (paywalled).
- The 32 well-index (`XXXXX`) and 96 subarray-index (`YYYYYYYY`) sequences of STRT-seq-2i — not in the paper, preprint or supplement.

---

## 1. What it is, and what is new

STRT ("single-cell tagged reverse transcription") was the first **early-pooling**
single-cell RNA-seq: the cell barcode goes on during RT, so 96 cells are pooled and
handled as one sample from the cDNA PCR on. 🟢 The barcode rides in on the
**template-switching oligo**, so every read starts at the cap site — STRT counts 5' ends
(TSS-resolved), not full transcripts. See
[template switching](../ref/concepts/template-switching.md) and
[reverse transcription](../ref/concepts/reverse-transcription.md).

| | New thing | Variant |
|---|---|---|
| 1 | Cell barcode **in the TSO**, pool after RT, one PCR tube for 96 cells | STRT 🟢 |
| 2 | Keep only **5' fragments** by streptavidin capture of the biotinylated cDNA end + a **restriction site in the oligo-dT** handle (SalI in STRT; PvuI in C1 🟡) that releases the 3'-end fragments | STRT 🟢 / C1 🟡 |
| 3 | **UMI in the TSO** (molecule counting); **all-RNA TSO** so it hydrolyses in Mg²⁺ during hot start and cannot prime in PCR; TSO shortened to 29 nt | C1 🟢 (Supplementary Note) |
| 4 | Read-1 handle = **Illumina P5 itself** on TSO and oligo-dT; cell barcode moved into a **Tn5 adapter** that carries P7, so the tagmented library needs **no library PCR** | C1 🟢 |
| 5 | **Dual index**: well index added by a PCR primer that docks on the TSO's TruSeq-Read-1 handle, subarray index by Tn5; single-end reads | 2i 🟢 |

## 2. Oligos

### 2a. STRT (Genome Res 2011) — 🟢 verbatim from the Methods

```
STRT-V3-T30    5′-biotin-AAGCAGTGGTATCAACGCAGAGTCGACT30VN-3′
STRT-V2-n      5′-AAGCAGTGGTATCAACGCAGAGTGCAGTGCTXXXXXXrGrGrG-3′     n = 1..96, XXXXXX = 6-bp barcode, rG = riboguanine
STRT-PCR       5′-biotin-AAGCAGTGGTATCAACGCAGAGT-3′
P2 adapter     5′-CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT-3′
               3′-PHO-GTTCGTCTTCTGCCGTATGCTCGAGAAGGCTAG-PHO-5′
library PCR 1  5′-AATGATACGGCGACCACCGAGATCTAAGCAGTGGTATCAACGCAGAGT-3′
library PCR 2  5′-CAAGCAGAAGACGGCATACGAG-3′
```

The 96 barcodes are in Supplemental Table 1, which we do not have. 🔴 The read-1
sequencing primer is called only "a custom primer" (Fig. 1 legend); its sequence is not
in the paper. 🔴

### 2b. STRT/C1 (Nat Methods 2014) — 🟢 verbatim from the oligo table (`MOESM268`)

```
C1-P1-RNA-TSO  Bio-AAUGAUACGGCGACCACCGAUNNNNNGGG        (all-RNA oligo)
C1-P1-T31      Bio-AATGATACGGCGACCACCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C1-P1-PCR-2    Bio-GAATGATACGGCGACCACCGAT
C1-TN5-U       PHO-CTGTCTCTTATACACATCTGACGC             (mix with C1-TN5-x to make adapter)
C1-TN5-x       CAAGCAGAAGACGGCATACGA <barcode> GCGTCAGATGTGTATAAGAGACAG     x = 1..96, 53 nt
               e.g. C1-TN5-1  CAAGCAGAAGACGGCATACGACGTCTAATGCGTCAGATGTGTATAAGAGACAG  (barcode CGTCTAAT)
```

The table lists all 96 C1-TN5 oligos with their 8-nt barcodes. No sequencing primers are
listed. 🔴

### 2c. STRT-seq-2i (preprint and Sci Rep Table 1) — 🟢 verbatim (Sci Rep spelling)

```
Lysis/RT
STRT-P1-T31            5′Bio-AATGATACGGCGACCACCGATCG-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
P1B-UMI-RNA-TSO        5′Bio-rCrTrArCrArCrGrArCrGrCrTrCrTrTrCrCrGrArTrCrT-rNrNrNrNrNrN-rGrGrG
PCR 1
DI-PCR-P1A             5′Bio-AATGATACGGCGACCACCGA
DI-P1A-idx[1–32]-P1B   5′Bio-AATGATACGGCGACCACCGAGATCTACAC-XXXXX-CTACACGACGCTCTTCCGATC
Tagmentation
STRT-Tn5-Idx[1–96]     CAAGCAGAAGACGGCATACGA-YYYYYYYY-GCGTCAGATGTGTATAAGAGACAG
STRT-TN5-U             5′PHO-CTGTCTCTTATACACATCTGACGC
PCR 2
P1_2nd_PCR             AATGATACGGCGACCACCGAGATC
P2_2nd_PCR             CAAGCAGAAGACGGCATACGAGAT
Sequencing
DI-Read1-Seq           ATGATACGGCGACCACCGAGATCTACAC-NNNNNN-CTACACGACGCTCTTCCGATCT
STRT-Tn5-U             5′PHO-CTGTCTCTTATACACATCTGACGC
DI_idxP1A-Seq          AATGATACGGCGACCACCGAGATCTACAC
```

Preprint and published table agree base for base; they differ only in hyphenation (the
preprint drops the separator hyphen in five oligos). 🟢
Names drift in the methods prose: "STRT-P1T31" / "STRT-P1-T31"; the second-PCR primers
appear as "4K-P1_2ndPCR" and "P2_4K_2ndPCR", the former at "200mM" (surely 200 nM). 🟡
Upstream calls the oligo-dT "C1-P1-T31" — it is the same sequence as C1's (computed). 🟡

### How they interlock — 🟡 (computed with `lib/`)

**STRT**

- The **SMART handle** `rt.SMART_HANDLE` (`AAGCAGTGGTATCAACGCAGAGT`, 23 nt) starts the
  oligo-dT, the TSO and *is* STRT-PCR — the Clontech SMART single-primer design.
- STRT-V3-T30 = handle + `CGAC` + T30 + `VN` (59 nt). The handle's final `GT` + `CGAC`
  spell **`GTCGAC` = SalI**, so the SalI site sits at the junction handle/poly(T).
- STRT-V2-n = handle + `GCAGTGCT` spacer + 6-nt barcode + `rGrGrG` = **40 nt**. It
  contains no SalI site of its own; a barcode could create one only for 1 of the 4096
  possible 6-mers (computed by enumeration) — whether any of the 96 does is unknowable
  without Supplemental Table 1. 🔴
- Library PCR 1 = `illumina.P5[:25]` (`AATGATACGGCGACCACCGAGATCT`) + the full SMART handle
  (48 nt): P5 grafted onto the TSO end.
- Library PCR 2 = `illumina.P7[:22]`.
- P2 adapter top = `illumina.P7[:22]` + `CTCTTCCGATCT` (= the last 12 nt of
  `illumina.TRUSEQ_READ2`), 34 nt with a 3' **T overhang**; the bottom strand is exactly
  the reverse complement of top[:-1] (33 nt). After ligation the insert is followed by
  `AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG` (= revcomp of the top strand, starting with
  `A` + `illumina.STEM`) — the classic single-end "genomic DNA" adapter read-through.

**STRT/C1**

- The handle is now the first **21 nt of the P5 flow-cell sequence + T**:
  `AATGATACGGCGACCACCGAT` = `illumina.P5[:20]` + `T`. **C1-P1-T31 and C1-P1-RNA-TSO share
  exactly these 21 nt**, so the cDNA has the same handle at both ends and one primer
  amplifies it.
- C1-P1-PCR-2 = `G` + `illumina.P5[:20]` + `T` (22 nt): the handle plus an extra
  untemplated 5' `G`.
- C1-P1-T31 = handle + `CG` + T31 (54 nt). `CCGAT`+`CG` spells **`CGATCG` = PvuI** at
  position 17 — the C1 equivalent of STRT's SalI site. The TSO end has no PvuI site
  unless the UMI starts with `CG`: **64 of 1024** 5-mers (6.25 %) create one (computed by
  enumeration over the TSO with flanks). If the PvuI step exists (🟡, upstream only),
  those molecules are cut off the bead and lost — a UMI-composition bias.
- C1-P1-RNA-TSO: 21-nt handle + 5-nt UMI + `GGG` = **29 nt**, all ribo — matches the
  Supplementary Note ("shortened the TSO to 29 bp"). 🟢
- C1-TN5-x = `illumina.P7[:21]` + 8-nt barcode + `GCGTCAGATGTGTATAAGAGACAG`. That
  24-nt tail is `nextera.ADAPTOR_S5[-24:]` = last 5 nt of s5 (`GCGTC`) + `nextera.ME`.
  **All 96** table rows have exactly this structure, with the barcode written **as named**
  in the "Barcode" column; the 96 barcodes are unique, minimum pairwise Hamming distance 4.
- C1-TN5-U is the **exact reverse complement of that 24-nt tail** (`ME_RC` + `GACGC`), and
  equals `nextera.INDEX2_PRIMER[:24]`. The duplex is 24 bp, not just the 19-bp ME.

**STRT-seq-2i**

- STRT-P1-T31 is identical to C1-P1-T31 (PvuI site and all).
- The TSO handle changes to **`CTACACGACGCTCTTCCGATCT` = `illumina.TRUSEQ_READ1[-22:]`**
  ("P1B"); TSO = 22 + 6-nt UMI + 3 = **31 nt**, all RNA.
- DI-P1A-idx-P1B = `illumina.P5` (29) + 5-nt well index + `TRUSEQ_READ1[-22:-1]` (21 nt,
  the TSO handle minus its last `T`) = 55 nt. It primes only on the TSO end.
- DI-PCR-P1A = `illumina.P5[:20]`: it primes the oligo-dT end (which starts with the same
  20 nt) and, once the index primer has run, the TSO end too. **Both** PCR primers are
  biotinylated, so both ends of the amplified cDNA carry biotin. 🟡
- STRT-Tn5-Idx and STRT-TN5-U have the same layout as C1-TN5-x / C1-TN5-U.
- P1_2nd_PCR = `illumina.P5[:24]`, P2_2nd_PCR = `illumina.P7` (24 nt).
  - P1_2nd_PCR matches the 5' (TSO-end) fragments over all 24 nt, but on the 3' (oligo-dT)
    fragments only the first 20 nt match: its last 4 nt `GATC` face `TCGT` — a 4-nt 3'
    mismatch. This, plus the absence of the P1B read-1 site, is what excludes 3'
    fragments in 2i; the methods mention **no PvuI step** (see §6). 🟡
  - P2_2nd_PCR's last 3 nt (`GAT`) sit where the Tn5 adapter has the first 3 barcode bases
    (adapter = `P7[:21]` + `YYYYYYYY`). Only barcodes beginning `GAT` match fully; in the
    C1 set (if reused — unknown 🔴) 1/96 does, 14 have 1, 42 have 2 and 39 have 3
    mismatches there. Consequences in §4. 🟡
- DI-Read1-Seq = `illumina.P5[1:]` (28 nt, P5 without its first `A`) + **6 N** + the 22-nt
  TSO handle = 56 nt. The template has only a **5-nt** well index there, so a 6-N primer
  is one base too long — upstream flags this too. Unchanged in the published Table 1. 🟡
- DI_idxP1A-Seq = `illumina.P5` exactly; it reads the 5-nt well index. STRT-TN5-U is
  reused as the Index 1 primer and reads the 8-nt subarray index.

## 3. Step by step

### STRT (2011) — 🟢 conditions from the Methods

1. **Capture plate**: 96-well, 5 µL STRT buffer (20 mM Tris pH 8.0, 75 mM KCl, 6 mM MgCl₂,
   0.02 % Tween-20) with 400 nM STRT-V3-T30 and 400 nM STRT-V2-n — a different barcoded
   TSO per well. One picked cell per well, frozen on dry ice.
2. **RT + template switch**: add 5 µL RT mix (4 mM DTT, 2 mM dNTP, 5 U/µL SuperScript II,
   **6 mM MnCl₂**, 2500 molecules of 8 synthetic control mRNAs); 10 °C 10 min, 42 °C 45 min.
   No heat lysis — it fragmented RNA and caused template-switching hotspots.
   First strand: `5'-bio-handle-CGAC-T30-VN-<cDNA>-CCC`; the TSO's `rGrGrG` pairs the `CCC`
   and RT copies the TSO, so the cDNA's other end gets handle-`GCAGTGCT`-<barcode>-`GGG`
   (as complement). 🟡 (standard [template switching](../ref/concepts/template-switching.md))
3. **Pool and clean**: MyOne carboxylate beads in 14 % PEG / 0.9 M NaCl added per well;
   all wells pooled, ethanol washed, eluted in 37 µL.
4. **Single-primer PCR, 20 cycles**: 200 nM biotin STRT-PCR, Advantage2; 95 °C 1 min;
   20 × (95 °C 15 s, 65 °C 30 s, 68 °C 4 min). Both ends carry the handle, so short
   products fold back on themselves (suppression PCR). 🟡 Both strands now 5'-biotinylated.
5. **Immobilise** on MyOne C1 streptavidin beads (2 M NaCl BWT), 10 min.
6. **Fragment** on beads: DNase I 0.0003 U/µL + 10 mM MnCl₂ (Mn²⁺ favours double-strand
   breaks), exactly 8 min RT, 120 µL. Each molecule leaves two bead-bound ends: the TSO
   end and the oligo-dT end. 🟡
7. **End repair** (NEBNext, 30 min RT), **A-tail** (Klenow exo⁻, 37 °C 30 min).
8. **Ligate P2 adapter + SalI-HF in one step**: 2 µM adapter, T4 ligase, 2 U/µL SalI-HF,
   NEBuffer 4, 37 °C 30 min. SalI cuts at the T30 handle and **releases the 3'-end
   fragments** from the beads 🟢; the TSO end has no site (barring an unlucky barcode). 🟡
9. **Library PCR, 12 cycles** on the beads: library PCR 1 + 2 (400 nM), Phusion; 98 °C
   30 s; 12 × (98 °C 10 s, 65 °C 30 s, 72 °C 30 s); 72 °C 5 min. AMPure; 1–3 ng/µL.
10. **Size select** 200–400 bp on 2 % E-gel.

### STRT/C1 — 🟢 only what the Supplementary Note says; the rest 🟡 (upstream only)

Changes from STRT, per the Supplementary Note 🟢: no well barcode in the TSO and TSO
shortened (36 → 29 nt; but see §6 — STRT's TSO is 40 nt by the 2011 sequence);
cDNA not purified, just diluted 1:10 before PCR; all-RNA TSO; tagmentation replaces
fragmentation/end-repair/A-tailing/ligation, and **library PCR is eliminated** (total PCR
33 → 21 cycles); TSO carries the Illumina read-1 adapter, Tn5 brings read 2. Capture,
imaging, wash, lysis, RT and PCR run on the C1 AutoPrep; tagmentation on a Biomek 2000.
The source of the Tn5 enzyme is not stated in our sources. 🔴

1. C1 capture, lysis, RT with C1-P1-T31 + C1-P1-RNA-TSO (UMI). 🟢 order; conditions 🔴.
2. cDNA PCR **21 cycles** with C1-P1-PCR-2. 🟢 cycle count; conditions 🔴.
3. **Tagment** each cell's cDNA with its own barcoded Tn5 (C1-TN5-x/U), 96 per run, pool. 🟢
4. Streptavidin capture + **PvuI** to drop 3'-end fragments; denature (NaOH, per upstream)
   to release the non-biotinylated strand. 🟡 (upstream only; plausible because
   C1-P1-T31 carries the PvuI site — computed). The released strand is the one the Tn5
   adapter was joined to, so it is contiguous from P7 to P5′ and needs no gap fill. 🟡
5. Pooled library sequenced **directly, without amplification** (samples 187ss, 187ss-2);
   one aliquot was amplified for comparison (187ds). 🟢 Upstream draws a P1/P2 library
   PCR here — true only for 187ds. 🟡

### STRT-seq-2i — 🟢 conditions from the methods (preprint = Sci Rep)

1. **Dispense or FACS** cells into the 9600-well array (96 subarrays × 100 wells);
   image (CellTracker Green), select single-cell wells; all later dispenses go only to
   those wells.
2. **Lysis** (Alt A): 50 nL with 500 nM STRT-P1-T31, dNTP, 2 % Triton X-100, DTT, RNase
   inhibitor; 72 °C 3 min. (Supp. Fig. S2c reports better yield at 95 °C; methods keep 72 °C.)
3. **RT** (Alt A): 85 nL with SuperScript II, 12.6 mM MgCl₂, 1.79 M betaine, 10.5 µM
   P1B-UMI-RNA-TSO; 42 °C 90 min. (Alt B, frozen plates: one 70 nL lysis-RT mix.)
4. **Indexed PCR**: one of 32 DI-P1A-idx-P1B primers per well (200 nM) + 100 nM DI-PCR-P1A,
   KAPA HiFi; 95 °C 3 min; 5 × (98 °C 30 s, 67 °C 1 min, 72 °C 6 min); 15 × (98 °C 30 s,
   68 °C 30 s, 72 °C 6 min); 72 °C 5 min. The first cycles graft P5 + well index onto the
   TSO end; after that DI-PCR-P1A amplifies both ends. 🟡
5. **Extract** by centrifuging the inverted array into a 96-well plate — one well per
   subarray, i.e. pooled by subarray.
6. **Tagment** each subarray pool with its own STRT-Tn5-Idx transposome (6.25 µM adapter +
   6.25 µM Tn5, 37 °C 1 h): 3 µL + 2 µL cDNA in 20 µL CutSmart, 55 °C 20 min.
7. **Streptavidin** (MyOne C1 in 0.5 % SDS binding buffer) 15 min, pool all 96,
   wash; **ExoSAP-IT** 37 °C 15 min removes leftover adapters; elute the **single-stranded
   library** in water, 70 °C 10 min. No PvuI step is described. 🟢
8. **PCR 2, 8 cycles**: P1_2nd_PCR + P2_2nd_PCR, KAPA HiFi; 98 °C 30 s, 65 °C 10 s, 72 °C
   20 s. AMPure 0.7×, then 0.5×/1× double size selection.

## 4. Final libraries — 🟡 (assembled from the oligos above; strand shown = Read-1 sense)

**STRT**

```
5'- P5[:25] · SMART handle (23) · GCAGTGCT · <bc 6> · GGG(…) · <mRNA 5' end, sense> · A · GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -3'
```

The last segment is revcomp(P2 adapter top) minus its first base; with the preceding `A`
from A-tailing it is the full 34 nt. The P7 end is only `P7[:22]` — the 2011 Solexa layout.

**STRT/C1**

```
5'- G · P5[:20] · T · <UMI 5> · GGG · <mRNA 5' end, sense> · ME_RC · GACGC · <revcomp bc 8> · revcomp(P7[:21]) -3'
```

(The leading `G` comes from C1-P1-PCR-2; upstream's final drawing omits it.) No full P5
or P7 — only the 20/21-nt Solexa P1/P2 graft sequences.

**STRT-seq-2i**

```
5'- P5 (29) · <well idx 5> · CTACACGACGCTCTTCCGATCT · <UMI 6> · GGG · <mRNA 5' end, sense> · ME_RC · GACGC · <revcomp subarray idx 8> · ATCTCGTATGCCGTCTTCTGCTTG -3'
```

Drawn as if P2_2nd_PCR is incorporated whole, so the last segment is `illumina.P7_RC`.
But the primer's last 3 nt (`GAT`) overlap the first 3 barcode bases (§2, interlock), so
two outcomes are possible: either the mismatched primer extends as is and its `GAT`
overwrites barcode bases 1–3 in the copies (the 8-nt index read would then end in `ATC`
for every subarray), or the polymerase's proofreading trims the mismatches and the
library keeps the full barcode followed by `revcomp(P7[:21])`. Which happens is not
stated. 🟡

## 5. Read layout / sequencing

| | Read 1 | Index read(s) | Instrument |
|---|---|---|---|
| STRT | custom primer from the P1 side 🟢; sequence 🔴 (upstream gives `GATCTAAGCAGTGGTATCAACGCAGAGTGCAGTGCT` = `P5[20:25]` + TSO body, computed 🟡). Reads `<bc6>`, then G's, then the 5' end of the mRNA. ~55 bp (cost estimate). | none — the barcode is in Read 1 | GA IIx 🟢 |
| STRT/C1 | **50 bp**: 5-nt UMI, then into the mRNA 5' end 🟢. Primer sequence 🔴 — must end at the handle's `…CCGAT`. | separate index read for the 8-nt barcode 🟢; primer presumably C1-TN5-U (as in 2i) 🟡; reads **revcomp** of the table barcode 🟡 (computed) | HiSeq 2000 🟢 |
| STRT-seq-2i | single-end, DI-Read1-Seq: 6-nt UMI, `GGG`, transcript 🟢. **45 bp** (methods) vs **48 bp** (Supp. Fig. S1) — inconsistent within the paper 🟡 | Index 1 = STRT-Tn5-U, 8 nt subarray (revcomp of `YYYYYYYY` 🟡); Index 2 = DI_idxP1A-Seq, 5 nt well (as written in the primer 🟡) 🟢 | HiSeq 2000/2500, SE 50-cycle kit 🟢 |

The reads begin at the TSO junction, so the first transcript base after the `GGG` is the
cap-proximal 5' end. 2i analysis trims up to nine G's and discards reads without three. 🟢
The 2011 Fig. 1 legend says each read has the barcode "followed by three to six Cs";
reading the TSO-sense strand they are G's (as 2i says) — a wording slip. 🟡

## 6. Upstream (scg_lib_structs) checked against the papers

Agreements 🟢/🟡: STRT TSO, STRT-PCR (= "ISPCR"), library PCR 1/2, P2 adapter top and
bottom sequences; all three C1 oligos and the Tn5 top/bottom; every 2i oligo including
the 6-N DI-Read1-Seq (upstream flags the 6 N as a likely error, which our computation
supports).

Disagreements:

- **STRT oligo-dT**: upstream `5'-Bio-TTAAGCAGTGGTATCAACGCAGAGTCGACT…T(30)VN` has an extra
  5' `TT` relative to the 2011 STRT-V3-T30 (computed: upstream[2:] equals the 2011 oligo).
  Probably the Nat Protoc 2012 version, which we could not read. 🟡
- **P2 adapter modifications**: 2011 gives the bottom strand phosphorylated at **both**
  ends (`3′-PHO-…-PHO-5′`) and no phosphorothioate; upstream writes the top as
  `…GATC*T` (`*` read as the usual phosphorothioate mark; the page does not define it) and no phosphates. 🟡 (protocol
  version difference, unverified)
- **STRT TSO length**: 2011 sequence gives 40 nt (31 + 6 + 3, computed); the C1
  Supplementary Note says the old TSO was 36 bp. 🔴 unexplained.
- **C1 read-1 primer**: upstream draws it as `5'- TTACTATGCCGCTGGTGGCTA` — that is the
  plain **complement** (not reverse complement) of `AATGATACGGCGACCACCGAT`, i.e. the
  handle written in the wrong orientation. 🟡 (computed)
- **C1 Tn5 bottom**: upstream draws only 19 nt (ME) annealed; C1-TN5-U is 24 nt and pairs
  over all 24. 🟡 (computed)
- **C1 library PCR**: upstream shows a P1/P2 amplification; the paper's main data were
  sequenced unamplified. 🟢
- **C1 PvuI / NaOH** and **2i PvuI**: neither is in our sources. For C1 we could not read
  the Online Methods (🔴); for 2i the methods describe streptavidin + ExoSAP + heat
  elution and no PvuI — upstream's 2i PvuI step looks carried over from C1. 🟡
- **2i oligo-dT name**: upstream "C1-P1-T31", paper "STRT-P1-T31" — same sequence. 🟡

## 7. Open questions

- 🔴 STRT barcodes (Supplemental Table 1) and the STRT read-1 primer.
- 🔴 STRT/C1 Online Methods: RT/PCR/tagmentation conditions, read-1 and index primers,
  whether PvuI is used. If it is, 6.25 % of UMIs (those starting `CG`) are lost. 🟡
- 🔴 2i well-index (32 × 5 nt) and subarray-index (96 × 8 nt) sequences; whether the
  subarray set is the C1 barcode set.
- 🟡 2i: how cleanly a 4-nt 3' mismatch on P1_2nd_PCR, with a proofreading polymerase,
  excludes the oligo-dT-end fragments (they cannot be read by DI-Read1-Seq anyway).
- 🟡 2i DI-Read1-Seq has 6 N over a 5-nt index — a typo in both versions, or a deliberate
  bulge?
- 🔴 Why DI-Read1-Seq omits the first `A` of P5.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the upstream page, the 2i preprint and the C1 supplements,
but **not the defining Genome Research paper**: its DOI prefix `10.1101/` is Cold Spring
Harbor's, shared by bioRxiv and Genome Research, and the tool routes every `10.1101/`
DOI to bioRxiv, which fails silently (no `(manual)` row either). The PMC full text
(PMC3129258) and the published 2i paper (PMC5703850) + supplement were fetched by hand and
converted with `tools/doctext.py`. `tools/scrape_primers.py` found the C1 table, the 2i
tables and STRT-PCR / adapter / library primers, but missed STRT-V3-T30 (written with the
`T30VN` shorthand) and STRT-V2-n (`XXXXXX` placeholder + `rGrGrG`), split
DI-P1A-idx-P1B at its `XXXXX` placeholder, and did not flag the `U`-containing
C1-P1-RNA-TSO as RNA. Reaction order had to be read from the methods.
