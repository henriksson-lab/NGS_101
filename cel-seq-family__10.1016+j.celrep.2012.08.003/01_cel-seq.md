# CEL-seq and CEL-seq2: barcoded 3' RNA-seq by T7 linear amplification

> **Evidence marking.** 🟢 verbatim from a primary source (paper, supplement, detailed
> protocol) · 🟡 derived, inferred, or only in a secondary source (the upstream
> scg_lib_structs page) · 🔴 not published / not available to us. Relationships marked
> 🟡 *(computed)* were worked out with `lib/` while writing this note. They are not yet
> asserted in a self-test, because this protocol has no `tools/` module yet (status
> `notes`).

**CEL-Seq**: Hashimshony T, Wagner F, Sher N, Yanai I. "CEL-Seq: single-cell RNA-Seq by
multiplexed linear amplification." *Cell Reports* 2 (2012) 666–673.
doi:[10.1016/j.celrep.2012.08.003](https://doi.org/10.1016/j.celrep.2012.08.003), PMID 22939981.
**We could not fetch this paper** (see below), so everything about the original protocol
comes from the CEL-Seq2 paper and the upstream page.

**CEL-Seq2**: Hashimshony T, Senderovich N, Avital G, Klochendler A, de Leeuw Y, Anavy L,
Gennert D, Li S, Livak KJ, Rozenblatt-Rosen O, Dor Y, Regev A, Yanai I. "CEL-Seq2:
sensitive highly-multiplexed single-cell RNA-Seq." *Genome Biology* 17 (2016) 77.
doi:[10.1186/s13059-016-0938-8](https://doi.org/10.1186/s13059-016-0938-8), PMC4848782
(open access). Data: GEO GSE78779. Pipeline: github.com/yanailab/CEL-Seq-pipeline.

This one note covers the whole family: original CEL-Seq, the intermediate "CEL-Seq + UMI"
primer that the CEL-Seq2 paper uses as its baseline, CEL-Seq2 in tubes or plates, and
CEL-Seq2 on the Fluidigm C1.

## Sources

Fetched by `tools/get_sources.py` into
`_data/sources/cel-seq-family__10.1016+j.celrep.2012.08.003/`, never committed:

| File | What | Used for |
|---|---|---|
| `upstream_CEL-seq_family.html(.txt)` | scg_lib_structs "CEL-seq / CEL-seq2" page | the only source for the Illumina TruSeq Small RNA oligos (RA5, RA3, RTP, RP1, RPI), the read primers, and how the original CEL-Seq library was built |
| `s13059-016-0938-8_PMC4848782.html(.txt)` | CEL-Seq2 paper, full text (PMC) | what changed from CEL-Seq, barcode design, C1 implementation, sequencing |
| `…_MOESM4_ESM.docx(.txt)` | Additional file 4: "CEL-Seq2 Protocol" (detailed bench protocol) | **all 96 CEL-Seq2 RT primers**, randomhexRT primer, every volume and temperature |
| `…_MOESM5_ESM.pdf(.txt)` | Additional file 5, Table S2 | **CEL-Seq, CEL-Seq + UMI, CEL-Seq2 primer designs**, library RT primer |
| `…_MOESM6_ESM.docx(.txt)` | Additional file 6: C1 mixes and steps | CEL-Seq2 on the C1 |
| `…_MOESM3_ESM.pdf(.txt)` | Table S1: C1 capture and barcode per chamber | not chemistry (barcode-to-chamber map) |
| `…_MOESM1/2_ESM.pdf(.txt)` | Figures S1, S2 | optimisation data only (RT enzyme, second-strand kit, clean-up) |
| `…_PMC4848782.xml` | full-text JATS XML | the Methods (CEL-Seq2 changes, barcode design, Sequencing, GEO accession); these sections are missing from the PMC `.html.txt` |
| `…_supplementary.zip`, `…Fig*_HTML.*` | duplicates / figures | not used |

Not fetched (🔴):

| What | URL | Why it matters |
|---|---|---|
| CEL-Seq paper (Cell Reports 2012) and its Supplemental Information | https://doi.org/10.1016/j.celrep.2012.08.003 · https://www.cell.com/cell-reports/fulltext/S2211-1247(12)00228-8 | the original methods: the 8-nt barcode list ("as previously published", Table S2 of CEL-Seq2), which Illumina kit and adapters were used, IVT kit details, read lengths. The publisher page returned HTTP 403 and Europe PMC has no full text. |

---

## 1. What it is

A **3'-tag, early-barcoded** single-cell RNA-seq method whose amplification step is
**in vitro transcription (IVT) with T7 RNA polymerase**, not PCR. 🟢 (CEL-Seq2 paper,
Background.) Each cell's mRNA is reverse transcribed with its own oligo-dT primer that
carries, 5' to 3': a T7 promoter, the Illumina small-RNA 5' adapter, (in CEL-Seq2 a UMI,) a
cell barcode, and an anchored poly(T). After second-strand synthesis, every cDNA carries a
double-stranded T7 promoter. Barcoded cells can then be **pooled** and amplified together
by one IVT. The amplified antisense RNA (aRNA) is fragmented and turned into an Illumina
library that keeps only the fragment at the barcoded 3' end. 🟢 (CEL-Seq2 Fig. 1a legend,
Methods; mechanism details 🟡)

| | New thing | Builds on / recurs in |
|---|---|---|
| 1 | **Linear amplification by IVT** after early barcoding, so cells are pooled before amplification | Eberwine aRNA amplification; recurs in MARS-seq, inDrop, DR-seq, sci-RNA-seq-like IVT variants |
| 2 | **T7 promoter on the RT primer**, upstream of the sequencing adapter | the 5' end of each aRNA is then the read-1 adapter |
| 3 | (CEL-Seq) 3' adapter added to fragmented aRNA by **RNA ligation** (TruSeq Small RNA RA3) | [small-RNA ligation](../ref/concepts/small-rna-ligation.md) |
| 4 | (CEL-Seq2) ligation replaced by **random-hexamer RT with a 5' tail** carrying the 3' adapter | [reverse transcription](../ref/concepts/reverse-transcription.md), random priming |
| 5 | (CEL-Seq2) **UMI** between adapter and barcode; shorter primer (92 → 82 nt) | UMIs (Kivioja 2012), Grün 2014 |

CEL-Seq2 lists the changes as: shorter primers with a UMI upstream of the barcode;
SuperScript II and the SuperScript II double-stranded cDNA kit reagents in place of the
Ambion MessageAmp II kit; bead clean-ups in place of columns; random priming with a tailed
hexamer in place of adapter ligation. 🟢 (Methods, items 1–4; the MessageAmp II and column details are from Results.) The paper reports
that removing the ligation raised the fraction of reads with barcodes that map from
60.9 % to 93.8 %. 🟢

## 2. Oligos

### 2.1 RT primers (one per cell)

🟢 Verbatim from Table S2 (`MOESM5_ESM.pdf.txt`, lines 1–21). PDF line wraps rejoined.
Unmodified, desalted DNA. No modifications are given in any source.

```
CEL-Seq        CGATTGAGGCCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATC[8 base barcode]TTTTTTTTTTTTTTTTTTTTTTTTV
CEL-Seq + UMI  CGATTGAGGCCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATCNNNNN[6 base barcode]TTTTTTTTTTTTTTTTTTTTTTTTV
CEL-Seq2       GCCGGTAATACGACTCACTATAGGGAGTTCTACAGTCCGACGATCNNNNNN[6 base barcode]TTTTTTTTTTTTTTTTTTTTTTTTV
```

The 96 CEL-Seq2 primers, 🟢 verbatim from the detailed protocol (`MOESM4_ESM.docx.txt`,
lines 51–146). Each primer is the CEL-Seq2 line above with the 6-nt barcode filled in,
for example:

```
1s   GCCGGTAATACGACTCACTATAGGGAGTTCTACAGTCCGACGATCNNNNNNAGACTCTTTTTTTTTTTTTTTTTTTTTTTTV
96s  GCCGGTAATACGACTCACTATAGGGAGTTCTACAGTCCGACGATCNNNNNNGAGTGATTTTTTTTTTTTTTTTTTTTTTTTV
```

Barcodes, as they appear in the primers (5'→3' on the RT primer), 🟢:

```
 1s AGACTC   2s AGCTAG   3s AGCTCA   4s AGCTTC   5s CATGAG   6s CATGCA   7s CATGTC   8s CACTAG
 9s CAGATC  10s TCACAG  11s AGGATC  12s AGTGCA  13s AGTGTC  14s TCCTAG  15s TCTGAG  16s TCTGCA
17s TCGAAG  18s TCGACA  19s TCGATC  20s GTACAG  21s GTACCA  22s GTACTC  23s GTCTAG  24s GTCTCA
25s GTTGCA  26s GTGACA  27s GTGATC  28s ACAGTG  29s ACCATG  30s ACTCTG  31s ACTCGA  32s ACGTAC
33s ACGTTG  34s ACGTGA  35s CTAGAC  36s CTAGTG  37s CTAGGA  38s CTCATG  39s CTCAGA  40s CTTCGA
41s CTGTAC  42s CTGTGA  43s TGAGAC  44s TGCAAC  45s TGCATG  46s TGCAGA  47s TGTCAC  48s TGTCGA
49s TGGTAC  50s GACATG  51s GATCAC  52s GATCTG  53s GATCGA  54s GAGTAC  55s AGACAG  56s AGACCA
57s AGTGAG  58s AGGAAG  59s AGGACA  60s CAACAG  61s CAACCA  62s CAACTC  63s CACTCA  64s CACTTC
65s CAGAAG  66s CAGACA  67s TCACCA  68s TCACTC  69s TCCTCA  70s TCCTTC  71s TCTGTC  72s GTCTTC
73s GTTGAG  74s GTTGTC  75s GTGAAG  76s ACAGAC  77s ACAGGA  78s ACCAAC  79s ACCAGA  80s ACTCAC
81s CTCAAC  82s CTTCAC  83s CTTCTG  84s CTGTTG  85s TGAGTG  86s TGAGGA  87s TGTCTG  88s TGGTTG
89s TGGTGA  90s GAAGAC  91s GAAGTG  92s GAAGGA  93s GACAAC  94s GACAGA  95s GAGTTG  96s GAGTGA
```

The protocol recommends barcodes 1, 4, 5, 9, 10, 23, 25, 26, 31, 46 when pooling as few as
10 cells. 🟢 The 8-nt CEL-Seq barcodes are 🔴: Table S2 says only that they are "as
previously published", which means in the CEL-Seq paper we could not fetch.

### 2.2 Library-construction oligos

CEL-Seq2, 🟢 (protocol line 33; Table S2 "Library RT primer"):

```
randomhexRT primer   GCCTTGGCACCCGAGAATTCCANNNNNN
```

The PCR primers are given only as "RNA PCR Primer (RP1, from Illumina kit)" and "uniquely
indexed RNA PCR Primer (RPIX)", "sequences available from Illumina". 🟢 The sequences
below come **only from the upstream page** and are therefore 🟡. They are the Illumina
TruSeq Small RNA kit oligos. The upstream page writes `[6-bp RPI]` for the index; we write
`<rpi>`.

```
RA5 (RNA 5' adapter)    GUUCAGAGUUCUACAGUCCGACGAUC                      CEL-Seq only, as a sequence reference
RA3 (RNA 3' adapter)    TGGAATTCTCGGGTGCCAAGG                           CEL-Seq only: ligated to aRNA
RTP (RT primer)         GCCTTGGCACCCGAGAATTCCA                          CEL-Seq only
RP1                     AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA
RPI1..48                CAAGCAGAAGACGGCATACGAGAT<rpi>GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
Read 1 seq. primer      GTTCAGAGTTCTACAGTCCGACGATC
Index read primer       TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
Read 2 seq. primer      GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
```

The upstream page lists RPI1–RPI48 as the bases written in the primer, for example
RPI1 `CGTGAT` and RPI2 `ACATCG` (see the upstream file, lines 58–106).

The CEL-Seq2 paper does confirm 🟢 that the original CEL-Seq **ligated** the second adapter
and that CEL-Seq2 "alleviates the need for Illumina's TruSeq Small-RNA kit". It says
nothing about the RA3 or RTP sequence, and it never names RA5. The upstream author says
that the CEL-Seq paper was "not entirely clear" and that the kit identity was inferred from
oligo names, then confirmed by CEL-Seq2.

### 2.3 How the oligos interlock: 🟡 (computed with `lib/`, `chemdraw.revcomp`, `illumina`)

Segments of the RT primers, 5'→3':

| Segment | CEL-Seq (92 nt) | CEL-Seq + UMI (95 nt) | CEL-Seq2 (82 nt) |
|---|---|---|---|
| 5' spacer | `CGATTGAGGCCGG` (13) | same | `GCCGG` (5) |
| T7 promoter `TAATACGACTCACTATAGGG` | 20 | 20 | 20 |
| Illumina small-RNA 5' adapter | full RA5 `GTTCAGAGTTCTACAGTCCGACGATC` (26) | same | **3'-truncated RA5** `AGTTCTACAGTCCGACGATC` (20) |
| UMI | none | `NNNNN` (5) | `NNNNNN` (6) |
| cell barcode | 8 | 6 | 6 |
| anchored oligo-dT | `T24V` (25) | same | same |

- The lengths add up to the paper's own figures: 92 nt for CEL-Seq and 82 nt for
  CEL-Seq2. 🟢 (Results: "from 92 to 82 nucleotides"). The intermediate primer comes to
  95 nt. 🟡
- The three G's at the end of the T7 promoter are the first transcribed bases, so each
  **aRNA begins** `GGGGUUCAGAG…` (CEL-Seq) or `GGGAGUUCUAC…` (CEL-Seq2). 🟡 (standard T7
  class III promoter, +1 = first G after `TATA`)
- **RP1 = `illumina.P5` + the first 21 nt of RA5** (`GTTCAGAGTTCTACAGTCCGA`); 50 nt.
  - On CEL-Seq aRNA-cDNA, all 21 of those bases pair.
  - On CEL-Seq2, the RT primer kept only the last 20 nt of RA5. RP1's 3'-terminal
    **16 nt** (`GAGTTCTACAGTCCGA`) pair, because the T7-derived G stands in for RA5's G.
    The next base along RP1 is a mismatch.
  - After the first PCR cycle, RP1's tail restores the full 26-nt RA5. **The P5 end of
    a CEL-Seq2 library is therefore identical to a CEL-Seq library**, and the standard
    small-RNA Read 1 primer (= RA5) still works.
- **randomhexRT tail = RTP**: the first 22 nt of `GCCTTGGCACCCGAGAATTCCANNNNNN` are exactly
  RTP. RTP is `G` + revcomp(RA3). CEL-Seq2 therefore puts the same sequence onto the
  cDNA at RT that CEL-Seq put there by RA3 ligation followed by RTP priming.
- **RPI 3' tail = `GTGACTGGAGTT` + revcomp(RA3)**. RPI has the same sense as RTP, so it
  primes on the RP1-extended copy of the randomhexRT-primed cDNA, which ends in
  revcomp(RTP). There RPI's 3'-terminal **21 nt** pair. The 22nd base mismatches, because the RTP-derived `C`
  meets the RPI's `T`. After PCR the library carries the RPI version.
- The **index read primer** is exactly revcomp of the RPI 3' tail and begins with RA3.
  The **Read 2 primer** is exactly the RPI 3' tail.
- **RPI starts with `illumina.P7`**. The index read reads the reverse complement of the
  bases written in the primer: RPI1 `CGTGAT` → read as `ATCACG`, RPI2 `ACATCG` →
  `CGATGT`, and so on. (Those read-out forms match how Illumina names TruSeq Small RNA
  indices, from general knowledge, unchecked against an Illumina document. 🟡)
- **CEL-Seq2 barcodes** (96, checked against the paper's design rules, 🟡 computed): all
  96 are distinct. The minimum pairwise Hamming distance is 2, as the paper says 🟢. The GC
  content is 50 % for every barcode, inside the stated 33–67 %. No barcode ends in T, as
  stated. The last base is A, C or G, 32 times each. This last-base rule keeps the
  barcode from running into the oligo-dT. 🟡 (reason inferred)
- The oligo-dT anchor is **`V` alone, not `VN`**. One anchoring base fixes the primer at
  the poly(A) junction. The first non-A base of the transcript is the only one forced.
  🟡 (compare `rt.oligo_dt`, default anchor `VN`)

## 3. Step by step

### 3.1 CEL-Seq2, tubes or plates: 🟢 from the detailed protocol (`MOESM4`) unless marked

1. **Cell into primer.** Either transfer the cell with a micro-pipette into a 0.5 µL drop on a tube cap, remove
   the excess liquid and freeze in liquid nitrogen, or sort it into a plate holding 1.2 µL primer mix. The primer mix is 1 µL
   primer at 25 ng/µL, ERCC spike-in, and 0.5 µL 10 mM dNTP, made up to 6 µL with water.
   That works out to about 5 ng of RT primer per cell. 🟡 (computed). The ERCC spike-in is
   the equivalent of 1 µL of a 1:1,000,000 dilution per cell. Freezing lyses the cell.
2. **Anneal.** 65 °C for 5 min (on a cap: 2 × 2.5 min with a spin in between), then onto ice.
3. **RT.** Add 0.8 µL containing First Strand buffer, DTT, RNase inhibitor and SuperScript
   II. Incubate 42 °C 1 h, then heat-inactivate at 70 °C for 10 min. First strand:
   `RT primer · cDNA (antisense)`.
4. **Second strand.** Add 10 µL Second Strand buffer, dNTP, E. coli DNA ligase, E. coli
   DNA Pol I and RNase H. Incubate 16 °C 2 h. This is RNase H nicking plus Pol I nick
   translation plus ligase, i.e. Gubler–Hoffman. 🟡 (named by us). The T7 promoter is now
   double-stranded. The cells are barcoded from this point.
5. **Pool** the cells that go to one IVT. Clean up with AMPure XP plus PEG/NaCl bead buffer,
   or for one or a few reactions 1.2× beads. Elute in 6.4 µL. The clean-up may be
   replaced by 65 °C for 20 min.
6. **IVT** (MEGAscript T7): A, G, C and U, 10× buffer, T7 enzyme, 37 °C for 13 h. The
   product is aRNA: antisense, starting with the T7 G's and the 5' adapter.
7. **ExoSAP-IT**, 37 °C 15 min, "to remove primers".
8. **Fragment** the aRNA: 22 µL aRNA + 5.5 µL buffer (200 mM Tris-acetate pH 8.1, 500 mM
   KOAc, 150 mM MgOAc), 94 °C 3 min. Stop with 2.75 µL 0.5 M EDTA. Clean up with 1.8×
   RNAClean XP. Expected yield is 500–1000 pg/µL from ~0.1 ng RNA, peaking at ~500 bp.
9. **Library RT**: 5 µL aRNA + 1 µL randomhexRT + 0.5 µL dNTP, 65 °C 5 min, onto ice. Add
   SuperScript II mix, then 25 °C 10 min and 42 °C 1 h. Only the 5'-most fragment of each
   aRNA carries the adapter, UMI and barcode. The hexamer adds the RTP sequence at the
   other end. 🟡
10. **PCR**: Phusion HF master mix with RP1 and one RPI per library. 98 °C 30 s; 11 cycles of
    (98 °C 10 s, 60 °C 30 s, 72 °C 30 s); 72 °C 10 min. Up to 15 cycles if the aRNA was low.
    Clean up twice with 1× AMPure XP. Expect at least ~1 ng/µL, peaking at 200–400 bp.
    Only fragments that carry RA5 (from the RT primer) at one end and the hexamer tail at
    the other are amplified exponentially. Internal aRNA fragments, which lack RA5, are
    not. 🟡

### 3.2 CEL-Seq2 on the Fluidigm C1: 🟢 from `MOESM6` and the Methods

- After capture, a different diluted CEL-Seq2 primer is loaded into each of the 96 outlets
  and goes to chamber 1. The primer stock is 25 ng/µL, mixed 2 µL with 8 µL Lysis Mix.
- **Chemical lysis**: NP-40, Tris pH 8.4, EDTA, SUPERase-In and ERCC. 65 °C 5 min. No
  freezing.
- **RT**: here MessageAmp kit reagents are used (First Strand buffer and **ArrayScript**),
  not SuperScript II. 42 °C 120 min.
- **Second strand**: MessageAmp DNA polymerase plus RNase H, 16 °C 120 min, then
  **heat-inactivation at 65 °C 20 min in place of the bead clean-up**.
- **IVT on chip**: MessageAmp T7, 37 °C 12 h. The IVT is per cell, not pooled.
- Harvest, then pool 5 µL from each outlet. Clean up with 1.8× RNAClean XP, then ExoSAP,
  fragmentation and library preparation "as in regular CEL-Seq2".
- Table S1 maps C1 chamber to barcode. The barcodes on the chip run 1–48, each used
  twice. 🟢 (table)

### 3.3 Original CEL-Seq (2012): 🟡/🔴, reconstructed without the primary paper

From the CEL-Seq2 paper 🟢: the MessageAmp II kit for RT, second strand and IVT; column
clean-ups; a 5'-adapter-carrying RT primer with an 8-nt barcode and no UMI. The aRNA was
made into a library by ligating the second adapter, then RT, then "a few PCR cycles".

From the upstream page only 🟡: fragmented aRNA is ligated at its 3' end to RA3. The upstream
author guesses the RA3 3' end is blocked. RTP primes RT on the ligated RA3, and RP1 + RPI
PCR follows. Ligase, conditions, fragmentation and cycle numbers are 🔴, because they are
in the unfetched paper. See
[small-RNA ligation](../ref/concepts/small-rna-ligation.md) for how the TruSeq Small RNA
adapters behave in their native kit.

## 4. Final libraries: 🟡 (assembled from the oligos above; agrees with upstream)

Written 5'→3' on the strand that begins with P5. `<rpi>'` is the reverse complement of the
6 bases written in the RPI primer.

**CEL-Seq2:**

```
5'- P5 · RA5 (26) · UMI (6) · cell barcode (6) · T24 · V' · <aRNA insert, antisense of mRNA> · N6' · TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC · <rpi>' · P7' -3'
```

**CEL-Seq:**

```
5'- P5 · RA5 (26) · cell barcode (8) · T24 · V' · <aRNA insert> · TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC · <rpi>' · P7' -3'
```

- P5 = `illumina.P5` (29 nt). P7' = `illumina.P7_RC` (24 nt). The 33-nt 3' adapter is
  RA3 + `AACTCCAGTCAC`. Fixed adapter bases total 130 nt for CEL-Seq2 and 126 nt for
  CEL-Seq, counting the UMI/barcode but not T24V, the insert or N6. 🟡 (computed)
- In CEL-Seq2, `N6'` stands for the random-hexamer-primed bases, which belong to the
  insert. In CEL-Seq, the ligation joins the insert directly to RA3.
- **The aRNA's 5' G's** (`GGGG`/`GGG` from the T7 promoter) are **not** in the final
  library. RP1 primes inside RA5, and the cDNA 3' end that copies those G's has no
  partner on P5. 🟡 (computed)
- **Disagreement with upstream (minor, 🟡):** the upstream page draws the CEL-Seq aRNA as
  starting `GGUUCAG…` and the CEL-Seq2 aRNA as `GAGUUC…`. With T7 transcription starting
  at the first G of `GGG`, we get `GGGGUUCAG…` and `GGGAGUUC…`. The difference does not
  reach the final library.
- In its step diagrams, the upstream page writes the CEL-Seq2 second strand as
  `[barcode][UMI]` in top-strand order. That is the reverse complement position of
  `UMI · barcode`, so it is consistent with our layout.

## 5. Sequencing

🟢 (CEL-Seq2 Methods and protocol): paired-end, **Read 1 = 15 nt, index read = 7 nt,
Read 2 = 36 nt**, run on a HiSeq 2500 in rapid mode. Libraries "should be considered
Small-RNA libraries". They were stated to be incompatible with HiSeq v4 high-throughput
reagents. Loading was 8 pM on HiSeq v3 and 12 pM on MiSeq or rapid mode.

| Read | Primer | Reads | |
|---|---|---|---|
| Read 1 | small-RNA Read 1 primer = RA5 | CEL-Seq2: UMI (6) + barcode (6) + `TTT`; CEL-Seq: barcode (8) + T's | 🟡 primer from upstream; content computed |
| Index 1 (i7) | `TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC` | 6-nt `<rpi>'` (+1 cycle) | 🟡 |
| Read 2 | `GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` | starts at the random-primed end of the fragment and reads the cDNA (mRNA sense) toward the poly(A) | 🟡 |

The pipeline demultiplexes on the Read 1 barcode, takes the UMI from Read 1 as well, and
maps Read 2 with Bowtie2. 🟢 Read 1 is too short to reach the transcript. All biology is in
Read 2. 🟡

The read primers are **not** the TruSeq/Nextera ones in `lib/seqprimers.py`. They are the
small-RNA set, hence the paper's warning about kit compatibility. 🟡

## 6. Open questions

- 🔴 The **original CEL-Seq paper and supplement**: the 8-nt barcode list, the ligation
  chemistry (ligase, the RA3 3' block, fragmentation), PCR cycles and read lengths. Fetch
  by hand from the URLs above.
- 🔴 The **5' spacers** (`CGATTGAGGCCGG`, shortened to `GCCGG` in CEL-Seq2) are not
  explained. They are presumably there to let T7 RNA polymerase bind a promoter near a
  duplex end. 🟡
- 🔴 What ExoSAP-IT removes after the IVT. The text says "primers": presumably leftover
  single-stranded RT primer, since exonuclease I acts on ssDNA. 🟡
- 🟡 Every TruSeq Small RNA sequence (RA3, RTP, RP1, RPI, read primers) rests on the upstream
  page only. It should be checked against Illumina's "Illumina Adapter Sequences" document
  before going into a module.
- 🟡 Whether the C1 version really kept ArrayScript and MessageAmp, against the paper's
  statement that CEL-Seq2 is "kit-free". `MOESM6` lists MessageAmp reagents for the chip
  RT, second strand and IVT.

## 7. How this note was made (tool evaluation)

`get_sources.py` got the full CEL-Seq2 package (PMC text plus the Europe PMC supplementary
zip, unpacked, with `.txt` twins). It could not get the CEL-Seq Cell Reports paper (publisher
only, 403). `scrape_primers.py` found all 96 CEL-Seq2 RT primers and the Table S2 designs.
The bracketed placeholders (`[8 base barcode]`) split each Table S2 design into two hits.
None of the small-RNA oligos are in `lib/`, so the scraper recognised only P5/P7 in them.
`--find` **cannot verify oligos containing degenerate bases**: it strips every non-ACGT
character from the query, so `…GATCNNNNNNAGACTC…` is searched as `…GATCAGACTC…` and
reported absent. The full CEL-Seq2 primers were therefore verified with a literal `grep -F`
on the `.txt` (lines 51 and 146 of `MOESM4`).
