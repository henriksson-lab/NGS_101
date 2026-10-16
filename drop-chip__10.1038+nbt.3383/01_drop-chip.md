# Drop-ChIP — single-cell ChIP-seq by barcoding nucleosomes in droplets

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; the current protocol model and self-test live in `tools/`.

**Drop-ChIP** — Rotem A, Ram O, Shoresh N, Sperling RA, Goren A, Weitz DA, Bernstein BE.
"Single-cell ChIP-seq reveals cell subpopulations defined by chromatin state."
*Nature Biotechnology* 33, 1165–1172 (2015). doi:[10.1038/nbt.3383](https://doi.org/10.1038/nbt.3383)
(PMID 26458175, PMC4636926). The paper points to a web portal with an interactive flow
chart (pubs.broadinstitute.org/drop-chip), which was not fetched. No other papers for
this method are listed in `catalogue/scg_lib_structs.tsv`.

Sources read (fetched by `tools/get_sources.py`, into
`/husky/henriksson/NGS101_papers/sources/drop-chip__10.1038+nbt.3383/`, never committed):

| File | What | Used for |
|---|---|---|
| `nbt.3383_PMC4636926.xml.txt` and `.html.txt` | original PMC article with Online Methods | **all reaction steps**, SC-PCR primers, read layout |
| `nbt.3383_41587_2015_BFnbt3383_MOESM11_ESM.pdf(.txt)` | Supplementary Information (figures, Suppl. Tables 1, 3, Notes) | buffer recipes (Suppl. Table 1), SC-PCR program, visually inspected Supplementary Fig. 2A/B |
| `nbt.3383_41587_2015_BFnbt3383_MOESM12_ESM.xlsx(.txt)` | Supplementary Table 2 | **all 1,152 barcode adaptors** (sheet `IDT_COA_07_26_12_RAM_convert`) |
| `nbt.3383_41587_2015_BFnbt3383_MOESM10_ESM.xlsx`, `MOESM13_ESM.xlsx` | Suppl. Tables 4, 5 (signatures, clustering) | not chemistry |
| `nbt.3383_41587_2015_BFnbt3383_MOESM14_ESM.zip` | Suppl. File 1, three chip-design `.dwg` files | member inventory checked; device designs, not chemistry |
| `MOESM15/16_ESM.avi` | Suppl. Movies 1, 2 | frames checked; droplet devices, not chemistry |
| `upstream_Drop-ChIP.html(.txt)` | Teichmann-lab scg_lib_structs page | orientation only; never primary evidence |

The JATS XML and PMC Cloud NIHMS supplements are now cached; the latter duplicate or
overlap the publisher supplements above. 🟢 **Supplementary Fig. 2A/B was inspected
visually in the original PDF, page 2.** The extracted text contains only its legend;
the bases in the image required separate transcription, recorded below.
🔴 **The exact final Illumina adaptor and library-primer identities remain unresolved:**
Online Methods names diluted Illumina adapters and Illumina primers but no product or
revision that identifies a matching vendor sequence source.

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

(The paper gives the two sequences without names; SC-PCR1 and SC-PCR2 are labels in
this note. No modifications are given, so unmodified DNA is inferred. 🟡)

### Illumina adaptor and library primers — 🔴 unresolved from checked sources

The original Online Methods names "Illumina adapters diluted 1:150" and "Illumina
Primers at 25 µM" without their product identity, revision, bases or index set
(JATS line 185). A matching vendor guide cannot yet be selected. The outer arms in
the model therefore remain inferred placeholders. No specific TruSeq construct is
established by these sources.

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

## 4. Final library — 🟡

The barcode-adaptor sequences in Supplementary Table 2 and the BciVI step support
the inner constant/barcode/TTAA region on both sides of the nucleosomal insert.
Supplementary Fig. 2B illustrates four orientations. The exact sequence outside
those inner regions depends on the unselected Illumina adapter and primer reagents;
the model leaves both outer arms unresolved rather than assigning P5, P7 or a
particular index sequence.

🟢 The original Supplementary Fig. 2A (PDF page 2) prints this illustrative top
strand in orientation 1, grouped by its colour-coded blocks:

```
5′-TTAA TTTCAATG GTATCC GGGGGACC TTAATTAA GGTGGGGG GGATAC GTAACTTT TTAA-3′  printed top
3′-AATT AAAGTTAC CATAGG GGGGGTGG AATTAATT CCAGGGGG CCTATG CATTGAAA AATT-5′  printed bottom
3′-AATT AAAGTTAC CATAGG CCCCCTGG AATTAATT CCACCCCC CCTATG CATTGAAA AATT-5′  computed complement
```

Orientation 2 swaps the two 8-nt barcode blocks between strands. Fig. 2B shows the
derived post-BciVI constants `GGGGGGGATAC` on the left and `GTATCCGGGGGA` on the
right, the two 8-nt barcode alternatives in each of four ligation orientations,
and `TTAA` at both insert junctions. It labels the outer regions only **“Illumina
adapter”**, without their bases. The pictured `TTTCAATG`/`GTAACTTT` pair is an
illustrative figure example and was not found among the 1,152 listed forward
oligos; the model continues to use listed row 1 rather than substituting the
figure barcode for a published production oligo. 🟢 for image text, 🟡 for its
relation to the production list.

🟡 **Published Fig. 2A pairing discrepancy:** both bottom green priming-site blocks
are printed as `GGGGGTGG` and `CCAGGGGG` in each orientation, whereas their top
partners `GGGGGACC` and `GGTGGGGG` require `CCCCCTGG` and `CCACCCCC`
respectively. The published Supplementary Table 2 forward strands share the
top-strand constant, so the model uses one listed top oligo and computes its
Watson–Crick complement for the duplex. The intended reason for the lower
figure motifs is unresolved; they cannot be used as a physically paired duplex.

## 5. Sequencing

🟢 HiSeq 2500, **2 × 60** paired end, ~320 M reads per lane. "The first 11 sequencing
cycles are dark" to avoid low-complexity failure on the constant adaptor bases.
Barcodes are expected at bases 1–8 of read 1 and 12–19 of read 2, followed by the TTAA
half of the PacI site, then genome. Alignment: `bowtie2 -X 1000 --trim5 23` (mm9); only
pairs whose two barcodes match are kept (with a rule merging two barcodes that co-label
≥10 % of each other's reads, for drops that fused with two barcode drops). 🟢

🟢 Fig. 2B places the blue Read 1 interval at the barcode after the 11-nt outer
constant, while the blue Read 2 interval includes its 11-nt outer constant before
the barcode. This matches the Methods coordinates and the 23-nt Read 2 prefix
(11 constant + 8 barcode + `TTAA`) trimmed during alignment. 🟡 The figure does not
label the instrument's dark-cycle instruction or sequencing primer, so allocating
the stated 11 initial dark cycles specifically to Read 1 remains an inference.

## 6. Source discrepancies and boundaries

The main-paper Fig. 3A legend calls the barcode adapters **64 bp** (JATS line 547),
whereas all 1,152 Supplementary Table 2 forward sequences are **60 nt**. The model
uses the actual listed 60-nt strands; an unlisted four-base extension cannot be
inferred. The Methods describes an A overhang at the **5′ end** (JATS line 185),
while the BciVI cut geometry from the printed recognition motif yields a **3′-A**
overhang in the model. The published orientation wording remains unresolved.
The separate Fig. 2A lower-strand priming-site mismatch is documented in §4.

## 7. Open questions

- 🔴 Exact Illumina adaptor, library primers and sample-index set.
- 🔴 How the BA duplex was formed (complementary oligo not listed).
- 🔴 Dark cycles on read 1 only or on both reads (§5).
- 🔴 Whether the 13 barcodes with a built-in PacI site at the insert boundary (and the 6
  with an extra BciVI site) were usable.
- 🔴 Sequencing-primer bases and the exact instrument configuration for the 11 dark
  cycles; Fig. 2 resolves the molecular read boundaries, not those run settings.
