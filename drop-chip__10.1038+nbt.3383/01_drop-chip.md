# Drop-ChIP — single-cell ChIP-seq by barcoding nucleosomes in droplets

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet.

**Drop-ChIP** — Rotem A, Ram O, Shoresh N, Sperling RA, Goren A, Weitz DA, Bernstein BE.
"Single-cell ChIP-seq reveals cell subpopulations defined by chromatin state."
*Nature Biotechnology* 33, 1165–1172 (2015). doi:[10.1038/nbt.3383](https://doi.org/10.1038/nbt.3383)
(PMID 26458175, PMC4636926). The paper points to a web portal with an interactive flow
chart (pubs.broadinstitute.org/drop-chip), which was not fetched. No other papers for
this method are listed in `catalogue/scg_lib_structs.tsv`.

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/drop-chip__10.1038+nbt.3383/`, never committed):

| File | What | Used for |
|---|---|---|
| `nbt.3383_PMC4636926.html(.txt)` | author manuscript on PMC, with Online Methods | **all reaction steps**, SC-PCR primers, read layout |
| `nbt.3383_41587_2015_BFnbt3383_MOESM11_ESM.pdf(.txt)` | Supplementary Information (figures, Suppl. Tables 1, 3, Notes) | buffer recipes (Suppl. Table 1), SC-PCR program, legend of Suppl. Fig. 2 |
| `nbt.3383_41587_2015_BFnbt3383_MOESM12_ESM.xlsx(.txt)` | Supplementary Table 2 | **all 1,152 barcode adaptors** (sheet `IDT_COA_07_26_12_RAM_convert`) |
| `nbt.3383_41587_2015_BFnbt3383_MOESM10_ESM.xlsx`, `MOESM13_ESM.xlsx` | Suppl. Tables 4, 5 (signatures, clustering) | not chemistry |
| `nbt.3383_41587_2015_BFnbt3383_MOESM14_ESM.zip` | Suppl. File 1, chip designs (`.dwg`) | not opened (CAD) |
| `MOESM15/16_ESM.avi` | Suppl. Movies 1, 2 | not chemistry |
| `upstream_Drop-ChIP.html(.txt)` | Teichmann-lab scg_lib_structs page | step-by-step drawing; checked below |

Not obtained: the Europe PMC JATS XML (HTTP 500) and the four NIHMS supplement copies
(download gate) — these duplicate the Springer files above. 🔴 **The sequences printed
in Supplementary Fig. 2A/B** (the adaptor in both orientations and the four final
configurations) are an image in `MOESM11_ESM.pdf`; only the legend came through as text.
🔴 **The Illumina adaptor and primers used for the library step are never given as
sequences** in any fetched file.

---

## 1. What it is

Single cells are encapsulated in ~50 µm drops together with detergent and **MNase**,
which cuts chromatin into mainly mono-nucleosomes inside each drop. A second emulsion
holds one double-stranded **barcode adaptor** per drop (1,152 different ones). A
three-inlet merger fuses one cell drop, one barcode drop and a squirt of labeling buffer
(end repair + **blunt ligation**), so every nucleosomal DNA end in the drop gets the same
barcode. Labeled chromatin from ~100 cells is pooled with **unlabeled carrier chromatin
from another species** (human K562) and only then immunoprecipitated, so the IP works at
bulk scale. 🟢 (main text, Online Methods)

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Barcode before IP**: label chromatin per cell, pool, then do one ChIP with carrier chromatin | later single-cell ChIP / CUT&Tag schemes that index first and pull down in bulk |
| 2 | A **symmetric** adaptor ligated blunt to *both* ends, with the barcode on both sides, so a read pair checks itself | — |
| 3 | **PacI to destroy adaptor concatemers** (two adaptors ligated back to back recreate TTAATTAA) and **BciVI** (Type IIS) to cut the adaptor down to a 1-nt A overhang for TA-style Illumina ligation | — |

Fragmentation is by MNase, not by Tn5 (contrast [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md),
where fragmentation and adaptor addition are one step). The Discussion lists replacing
MNase, better ligation and barcoded beads as the improvements needed. 🟢

## 2. Oligos

### Barcode adaptors (BA) — 🟢 Supplementary Table 2

1,152 rows, column "forward sequence", each 60 nt. First and last as written:

```
1      TTAAGGGCTTTCGTATCCGGGGGACCTTAATTAAGGTGGGGGGGATACCTTTCGGGTTAA
1152   TTAAACTTATAGGTATCCGGGGGACCTTAATTAAGGTGGGGGGGATACGATATTCATTAA
```

Ordered from IDT (sheet name "IDT_COA…"), suspended in Quick Ligase buffer at 500 µM in
384-well plates. 🟢 Only the forward strand is listed; how the double-stranded adaptor
was made (separate complementary oligo, annealed) is not stated. 🔴

### SC-PCR primers — 🟢 Online Methods ("Barcode and primer design")

```
SC-PCR1   TAAGGTGGGGGGGATAC    Tm 59.6
SC-PCR2   TAAGGTCCCCCGGATAC    Tm 59.6
```

(The paper gives the two sequences without names; "SC-PCR1 / SC-PCR2" is upstream's
labeling. No modifications are given, so they are taken as unmodified DNA. 🟡)

### Illumina adaptor and library primers — 🔴 not published

The methods say only "Illumina adapters diluted 1:150" and "Illumina Primers at 25 µM".
Upstream draws TruSeq forked adaptors (P5 + TruSeq Read 1 on one arm, TruSeq Read 2 arm +
8-nt i7 + P7 on the other, 3'-T overhang) and P5/P7 as primers, and **marks this itself as
an educated guess**, supported by finding `GATCGGAAGAGCACA…` (`illumina.STEM` + read-2 arm)
next to adaptor sequence in the public reads. 🟡 (upstream; plausible, because BciVI leaves
a single 3'-A that a T-tailed TruSeq adaptor would ligate to)

### How the adaptor is built — 🟡 (computed over all 1,152 rows)

Every row matches one pattern:

```
TTAA · <bc 8> · GTATCCGGGGGACCTTAATTAAGGTGGGGGGGATAC · <reverse of bc> · TTAA
```

- The 36-nt core is identical in all 1,152; only the 8-nt barcodes differ, and the
  1,152 left barcodes are all distinct (minimum Hamming distance 2).
- The **right barcode is the left one reversed** (not reverse-complemented) in all
  1,152. So the complementary strand reads, 5'→3', `TTAA · complement(bc) · … ·
  revcomp(bc) · TTAA`. The adaptor is "symmetric" in the sense that either end, on
  either strand, carries the same 8-nt identity — read as reverse(bc) or revcomp(bc)
  depending on which strand is read (see §5).
- **PacI** (`TTAATTAA`) sits in the middle of the core: `…GACC TTAATTAA GGTG…`. A
  second PacI site arises whenever two adaptors ligate end to end (`…TTAA`+`TTAA…`),
  which is how concatemers are destroyed. An adaptor ligated to genomic DNA does not
  make a site unless the genomic end happens to start with TTAA.
- **BciVI** (`GTATCC`, cutting downstream) sits once on each strand: `GTATCC` at
  position 13 of the forward strand and its reverse complement `GGATAC` at position 43.
  Each cuts on the side away from its own barcode, toward the central PacI site
  (`GTATCC` at 13 cuts after position 24 on the top strand). On a labeled fragment,
  where the adaptor half faces the insert, that means the cut removes the outer
  PacI/primer stub and leaves `GGGGGGGATAC · bc · TTAA` (or its counterpart) on each
  end of the insert.
- Without its G runs the core would be a perfect palindrome
  (`GTATCC·ACC·TTAATTAA·GGT·GGATAC`). The two G runs break the symmetry: the left arm has
  `GGGGG` where the mirror image has `CCCCC`. This is the paper's "5 Guanine nucleotides
  on each side of the barcode [that are] not complementary and form a loop" 🟢 (Online
  Methods), designed to stop the two ends of a labeled fragment from pairing with each
  other (a panhandle) during PCR. 🟡 (interpretation)
- Consequence for the primers: **SC-PCR1 occurs in the forward strand, SC-PCR2 in its
  reverse complement** (each exactly once). A fragment with adaptor halves at both ends
  presents SC-PCR1's sequence at one end and SC-PCR2's complement at the other in the
  most common case, so both primers are needed. 🟡

Barcode-set quirks 🟡 (computed):

- **13 barcodes** (rows 35, 76, 84, 104, 192, 245, 386, 470, 519, 614, 760, 847, 1048)
  create an extra PacI site right at the barcode–insert boundary: in 7 the left barcode
  starts with TTAA (`TTAA·TTAA…`), in 6 it starts with AATT, so the reversed right
  barcode ends with TTAA (`…TTAA·TTAA`). PacI digestion (done twice, see §3) would cut off the barcode
  and these labels should be lost. Whether these 13 are under-represented in the data is
  not reported. 🔴
- **6 barcodes** (rows 101, 237, 249, 529, 1017, 1100) contain an extra BciVI site.
- Read forms: across all 1,152, the 2,304 possible read forms (reverse(bc) and
  revcomp(bc)) are all distinct, minimum Hamming distance 2 — enough to detect, not to
  correct, a single error.
- Size: the paper's Fig. 3A legend calls them "64 bp" adaptors; every listed sequence
  is 60 nt. 🟡 disagreement, unexplained.

## 3. Step by step

Volumes and conditions 🟢 from the Online Methods and Supplementary Table 1 unless marked.

**In drops**

1. **Barcode emulsion**: 96-nozzle parallel drop-maker in a pressure chamber; each nozzle
   draws from one well of a 384-well plate, ~35 µm drops, ~10⁹ copies of one adaptor per
   drop. All 1,152 emulsions are pooled and mixed by rolling.
2. **Cell drops**: cells at 5 M/mL in PBS co-flowed 1:1 with 2× digestion buffer
   (Triton X-100, sodium deoxycholate, CaCl₂, sodium butyrate, MNase) into ~50 µm drops,
   loading tuned so ~1 in 6 drops has a cell. 4 °C 10 min (lysis), 37 °C 15 min (MNase).
3. **Fusion and labeling**: three-inlet merger; electrodes (100 V AC, 25 kHz) fuse one
   cell drop, one barcode drop and 2× labeling buffer (PEG 400, Quick Ligase buffer, ATP,
   End-It end-repair enzyme, a fast ligase, dNTPs, EGTA to stop MNase, sodium butyrate).
   End repair and blunt ligation in one mix; room temperature 2 h. Barcode adaptors go on
   both ends of each nucleosomal fragment, but adaptor concatemers also form (adaptor
   excess ~10⁹ vs ~10⁷ fragments). Fused drops are collected into a bed of carrier drops.

**Off chip**

4. **Split and stop**: 20 µL aliquots (~100 labeled cells) into wells, + 50 µL stopping
   buffer (EGTA, EDTA, Triton, deoxycholate) + 25 µL unlabeled K562 chromatin (~1 M cells,
   prepared by MNase in bulk) as carrier. Demulsifier (perfluoro-octanol), spin, take the
   aqueous phase.
5. **ChIP**: antibody (H3K4me3, H3K4me2 or H3.3) overnight 4 °C, protein-A magnetic beads;
   washes low-salt ×2, high-salt ×2, LiCl ×1, TE ×2; beads left in 21.5 µL TE.
6. **PacI #1, on the beads**: 37 °C 2 h, 65 °C 20 min. Cuts concatemers between adaptors
   and in the middle of every adaptor, leaving ~30-bp pieces that are removed by size
   selection. 🟢 What remains on a labeled fragment is half an adaptor at each end:
   `TAAGG · TGGGGGGG · ATAC · bc · TTAA · insert …` 🟡 (computed from the PacI site)
7. **Elute and clean**: 2× elution buffer (NaCl, SDS, EDTA), RNase 37 °C 20 min,
   proteinase K 37 °C 2 h, 65 °C 30 min; 1.5× AMPure XP.
8. **SC-PCR, 14 cycles** with both SC-PCR primers (Hercules Taq, 4 % DMSO, 3 µM each
   primer): 95 °C 3 min; 14 × (95 °C 30 s, 55 °C 10 s, 72 °C 1 min); 72 °C 10 min.
   1.1× AMPure. 🟢 The 17-nt primers anneal entirely inside the remaining adaptor half
   (`TAAGG` from the PacI half-site + core + `ATAC`). 🟡
9. **Dephosphorylate** (Antarctic Phosphatase, 37 °C 30 min) so that nothing but the
   BciVI-cut ends can take an Illumina adaptor ("to reduce unspecific adapter
   ligation"). 1.1× AMPure.
10. **BciVI**, NEB Buffer 4, 37 °C 1 h. Leaves "an A overhang" on every adaptor-labeled
    end (the paper says "at the 5′ end"). 🟢 The cut removes the primer-derived `TAAGGT`
    stub; the insert is now flanked by `GGGGGGGATAC · bc · TTAA` (or `CCCCCGGATAC · …`
    for the other strand) with a single unpaired A that is a 3' overhang, not 5'. 🟡
    (BciVI cuts GTATCC(N)6/5, giving a 1-nt 3' overhang; the paper's "5′ end" wording
    does not fit the enzyme)
11. **Illumina adaptor ligation**: 1.1× AMPure, evaporate to 4 µL, Quick Ligase + 1.5 µL
    adaptor (1:150), room temperature 15 min. Adaptor sequence 🔴.
12. **PacI #2** (same conditions) to cut concatemers formed in step 11; 0.7× AMPure.
13. **Library PCR, 14 cycles**: PfuUltra II, Illumina primers 0.5 µL of 25 µM; 95 °C
    3 min; 14 × (95 °C 30 s, 55 °C 30 s, 72 °C 1 min); 72 °C 10 min. 0.7× AMPure. 🟢
    The paper says a "sample" index is read to tell the ~100-cell pools apart 🟢; that it
    is added in this PCR is inferred 🟡. Primer sequences 🔴.

**Agreement with upstream**: steps 3–13 match upstream's drawing in order (ligation →
IP → PacI → SC-PCR → BciVI → Illumina ligation → PacI → P5/P7 PCR). Upstream omits the
dephosphorylation step and the carrier chromatin; it places PacI #1 after IP as the
paper does. 🟡

## 4. Final library — 🟡

Assembled from the BA table and the BciVI cut, with upstream's guessed TruSeq adaptor
(🟡 upstream / 🔴 unconfirmed) for the outer parts. The `lib/` check confirms upstream's
orientation-1 product starts with `illumina.NEBNEXT_UNIVERSAL_PRIMER` (P5 + TruSeq Read 1)
and contains `illumina.TRUSEQ_READ1`, `illumina.INDEX1_PRIMER` and `illumina.P7_RC`.
Top strand, orientation 1 (both adaptors in the forward orientation):

```
5'- P5 · TruSeq Read 1 · GGGGGGGATAC · reverse(bc) · TTAA · <nucleosomal DNA> · TTAA · bc · GTATCCGGGGGA · TruSeq Read 2' · i7' · P7' -3'
```

The four configurations (upstream's orientations 1–4, the paper's Suppl. Fig. 2B) come
from each of the two adaptors being ligated in either orientation. They differ only in
the 11-nt stub next to each Illumina arm (`GGGGGGGATAC` on one adaptor strand,
`CCCCCGGATAC` on the other, as read toward the insert) and in whether the barcode
appears as reverse(bc) or revcomp(bc). 🟡 (computed)

Upstream notes that configurations 2 and 3 (same stub at both ends) may amplify poorly by
suppression PCR; the paper does not discuss this. 🟡

## 5. Sequencing

🟢 HiSeq 2500, **2 × 60** paired end, ~320 M reads per lane. "The first 11 sequencing
cycles are dark" to avoid low-complexity failure on the constant adaptor bases.
Barcodes are expected at bases 1–8 of read 1 and 12–19 of read 2, followed by the TTAA
half of the PacI site, then genome. Alignment: `bowtie2 -X 1000 --trim5 23` (mm9); only
pairs whose two barcodes match are kept (with a rule merging two barcodes that co-label
≥10 % of each other's reads, for drops that fused with two barcode drops). 🟢

🟡 (computed with `lib/`, on the upstream-style library): after `illumina.TRUSEQ_READ1`,
exactly **11 nt** (`GGGGGGGATAC` or `CCCCCGGATAC`) precede the barcode; the same is true
after `illumina.TRUSEQ_READ2` (`CCCCCGGATAC` in orientation 1; the 3'-A from BciVI is
paired by the read-2 primer's terminal T). 11 + 8 + 4 (TTAA) = **23**, the `--trim5` value.
So read 1 with 11 dark cycles starts at the barcode (bases 1–8, as stated), and read 2
without dark cycles has the barcode at bases 12–19 (as stated). That would mean the
dark cycles were applied to read 1 only — but then trimming 23 nt from read 1 removes 11
genomic bases. Upstream draws 11 dark cycles on **both** reads. Which is right is not
resolvable from the text. 🔴

Index read: upstream's guess is the standard TruSeq i7 read with `illumina.INDEX1_PRIMER`,
8 nt. 🟡

## 6. Upstream page, checked

| Claim on upstream page | Paper / supplement | Verdict |
|---|---|---|
| BA = `TTAA·N8·GTATCCGGGGGACCTTAATTAAGGTGGGGGGGATAC·N8·TTAA`, 1,152 of them | matches all 1,152 rows of Suppl. Table 2 | agree 🟢 |
| Both barcodes drawn as the same `NNNNNNNN` | the right barcode is the left one **reversed** | upstream simplifies 🟡 |
| SC-PCR1 / SC-PCR2 sequences | identical to Online Methods | agree 🟢 |
| PacI after IP, before SC-PCR; again before final PCR | same | agree 🟢 |
| BciVI leaves an A overhang; TruSeq forked adaptors ligated | A overhang stated; adaptor identity not given | upstream's own "educated guess" 🟡 |
| Dephosphorylation before BciVI | in paper, not on upstream | upstream omits |
| 11 dark cycles on read 1 and read 2 | "first 11 cycles are dark"; barcode at read-2 bases 12–19 | conflict, see §5 🔴 |
| Adaptor length (drawn 60) | paper says 64 bp | paper's number inconsistent with its own table 🟡 |

## 7. Open questions

- 🔴 Exact Illumina adaptor, library primers and sample-index set.
- 🔴 How the BA duplex was formed (complementary oligo not listed).
- 🔴 Dark cycles on read 1 only or on both reads (§5).
- 🔴 Whether the 13 barcodes with a built-in PacI site at the insert boundary (and the 6
  with an extra BciVI site) were usable.
- 🔴 Sequences in Supplementary Fig. 2 (image only) — would settle the four
  configurations and the read-2 layout.
