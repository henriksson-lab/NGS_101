# sci-CAR-seq — research notes

## Sources read

- 🟢 Cao et al. 2018, DOI [10.1126/science.aau0730](https://doi.org/10.1126/science.aau0730), PMC6571013, main article workflow and figure legend.
- 🟢 Cao 2018 published Supplementary Materials (`science.aau0730_supp_NIHMS1033869-supplement-1.pdf`,
  PMC author-manuscript supplement, hand-downloaded 2026-10-17), Materials and Methods pp. 16–18
  (text twin lines 62–310; numbers lost in the twin's justified lines checked on the rendered
  pages 2–3). This is now the authority for reagents, conditions, primers and read recipes.
- 🟢 Cao 2018 supplementary table workbook (`…supplement-Supplimentary_table.xlsx`): sheets
  Table S1–S13 are throughput, differential-expression/accessibility, TF and linkage results.
  It contains no oligo or barcode table (searched for P5, P7, S5, ME and the RT-primer arm:
  no hits).
- 🟢 GEO GSE117089 family metadata (`geo-GSE117089-all-full.txt`, author-deposited, the
  article's data accession), sample extraction protocol at line 344488 (GSM3271040; the
  same text is repeated verbatim for GSM3271041–45). GEO truncates the field after
  "mRNA Second Strand S". Its primer, enzyme and cycle text agrees with the published
  Supplementary Materials (below); it is kept as history, not as the cited authority.
- 🟢 Amini et al. 2014, DOI [10.1038/ng.3119](https://doi.org/10.1038/ng.3119), Online Methods lines 116–122 and original publisher Supplementary Table 4 (`ng.3119_41588_2014_BFng3119_MOESM2_ESM.xlsx.txt`), correction notice lines 1–6 and oligo rows 9–13. This is the cited indexed-Tn5 predecessor; it does not identify sci-CAR's own oligo members.

## Molecular order

🟢 Intact nuclei first receive a well-indexed poly(dT) RT primer carrying a UMI, then a
well-indexed Tn5 transposome. Nuclei are pooled and redistributed by FACS. After cDNA
second-strand synthesis, each well is lysed and split: the RNA aliquot receives an
unindexed Tn5 site and a well-indexed PCR; the ATAC aliquot is amplified from its original
barcoded Tn5 ends with well-indexed PCR primers. RNA and ATAC libraries are pooled and
sequenced separately.

## Identifier semantics

The first RT and Tn5 indexes have different oligo sequences but encode the same first-well
coordinate. The RNA and ATAC PCR indexes likewise encode the same second/FACS well. The
two-coordinate lookup links both library types to one cell; no individual barcode segment
is physically shared between modalities.

## Oligo source review

🟢 The accessible main article (line 78) identifies four structural oligo families:
well-indexed poly(dT)/UMI RT primers, first-well barcoded ATAC transposomes,
unindexed cDNA transposomes, and well-indexed PCR primers for each modality.
The RNA and ATAC oligos in a well encode the same coordinate but have distinct
sequences. 🔴 The article does not print their ordered bases or barcode lists.
The published Supplementary Materials (read 2026-10-17, last section) print the RT
primer and the P5/P7 PCR primers with bracketed index placeholders, and the
supplementary workbook holds no oligo table; no barcode lists were found in either.
The schematic's barcode and PCR-index tokens are therefore placeholders (UMI 8 N and
RT index 10 bp are printed). The Tn5 adapter interiors remain inferred (see "Read
layout and inferred architecture" below).

### Published predecessor correction

🟢 Amini 2014 Supplementary Table 4's **Correction Notice** (extracted lines 1–6)
states that its original online version omitted the first five rows; the corrected
table was posted 20 May 2015. The corrected rows 9–13 are `P5_i5_1` through
`P5_i5_5`, with respective 8-base index values `TATAGCCT`, `ATAGAGGC`,
`CCTATCCT`, `GGCTCTGA`, `AGGCGAAG`. The original values were absent, not
conflicting bases. 🟡 The sci-CAR model uses generic placeholders for its first
well Tn5 index because the Cao 2018 supplementary oligo table remains gated;
it does not assign these five Amini members to sci-CAR.

## Read layout and inferred architecture (bounded pass 2026-10-17)

Sources re-read: Cao 2018 PMC6571013 JATS lines 78 (eight-stage workflow), 140
(transposase complexes "described in (7)"), 475 (Fig. 1 legend); Amini 2014
PMC4409979 JATS Online Methods lines 116–126; Amini Supplementary Table 4
(publisher MOESM2, `sci-atac-seq-family` bundle) correction notice and oligo rows.
The Cao supplementary Methods PDF and oligo workbook remain behind the PMC gate.

- 🟢 Fig. 1 legend (line 475): "RNA-seq: index2 and read1 cover the i5 index, UMI and
  RT barcode; index1 and read2 cover the i7 index and cDNA fragment. ATAC-seq: read1
  and read2 cover genomic DNA sequence. Index 1 and index 2 cover the Tn5 and PCR
  barcodes."
- 🟢 Line 78: both RNA PCR primers (Tn5-adaptor and RT-primer side) "bear a
  well-specific barcode that is the second RNA-seq index"; the ATAC primers likewise
  carry the second ATAC index.
- 🟡 Model correction: the earlier schematic placed the PCR-well and Tn5 barcodes inside
  Read 1/Read 2 and treated i5/i7 as unknown sample indices, and drew the RNA cDNA end
  as a TruSeq dA-ligation adapter. The model now puts the PCR-well index in i5 and i7
  (both features of the second-round cell barcode), the RT barcode and UMI in Read 1,
  and the unindexed-Tn5 end (ME/S7) on the Read 2 / Index 1 side.
- 🟡 Inferred sequences: the RT-primer Read 1 arm is drawn as TruSeq Read 1; the
  unindexed Tn5 end as Nextera ME + S7. The barcoded ATAC transposons use the Amini
  Table 4 architecture (S5–TCCACGC–8-nt index–GCGATCGAGGACGGC–ME and
  S7–CTGTCCCTGTCC–8-nt index–CACCGTCTCCGCCTC–ME, custom extended-mosaic-end read
  primers), justified only by line 140; the selftest checks that architecture against
  Amini rows P5_i5_1, P7_i7_1 and the three printed sequencing primers. Index 2 for ATAC
  is the derived reverse complement of the Read 1 primer.
- Open, not stated in the inspected article (gated supplement unread): sequencing-cycle numbers (the earlier 25/50 and 50/50
  read lengths had no inspected source and were removed), barcode/UMI/PCR-index lengths,
  the UMI versus RT-barcode order, RT primer and second-strand reagents.

## GEO protocol text (deferred follow-up 2026-10-17)

Routes for the gated supplement, retried: PMC instance URLs (download gate), Europe PMC
`supplementaryFiles` ("not open access"), science.org supplement (Cloudflare 403),
OpenAlex open locations (only the PMC author manuscript). The authors' GEO GSE117089
sample records were read instead (line 344488):

- 🟢 RT primer: `5′-ACGACGCTCTTCCGATCTNNNNNNNN[10bp index]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′`
  (25 µM, SuperScript IV, gradient 4–55 °C). This settles the UMI (8 N) before the
  10-bp RT index, and the 18-nt partial TruSeq Read 1 arm.
- 🟢 First-round ATAC: Nextera TD buffer plus "1 to 2 μL indexed TDE1 enzyme (Illumina)",
  55 °C 30 min, stopped with 40 mM EDTA / 1 mM spermidine. Pooled, DAPI-stained nuclei
  are sorted at 25 nuclei per well into 5 µL EB.
- 🟢 Second strand: NEB mRNA Second Strand Synthesis buffer and enzyme, 16 °C 180 min. No
  added primer is listed, so the model's earlier "second-strand primer" input was removed.
  The reaction is stopped with Zymo DNA binding buffer, cleaned with AMPure, and the eluate
  split into two plates.
- 🟢 RNA half: "i7 only TDE1 enzyme (25nM, Illumina)" at 55 °C for 5 min. Then SDS, BSA and
  P7 primer `CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG` at 55 °C for 15 min, then P5 primer
  `AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT` with NEBNext HiFi
  (72 °C 5 min; 18–22 cycles). Sequencing: NextSeq V2 75-cycle kit, Read 1 18, Read 2 52,
  Index 1 10, Index 2 10 cycles.
- 🟢 ATAC half: P5 primer `AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC` with the same
  P7 primer (72 °C 3 min; 18–22 cycles). Sequenced paired-end 51 bp with "custom primers as
  described before" and index reads covering the Tn5 and PCR barcodes.
- 🟡 Model corrections from this text. The Read 1 arm, UMI (8) and RT index (10) are no
  longer inferred. S5 and S7 are primer-established. Read lengths are RNA 18/52/10/10 and
  ATAC 51/51. Still inferred: the ME of the vendor "i7 only TDE1" adaptor, the PCR-index
  lengths (drawn as 10 nt; only RNA index cycles are given), the ATAC transposon interior
  (Amini architecture) and ATAC index cycles.
- Open at that time (resolved 2026-10-17, next section): the supplementary workbook, the
  RNA half after GEO's truncation, and whether the GEO text equals the published supplement
  (it does for every primer, enzyme and cycle number compared).

## Published Supplementary Materials (round deferred-followup-03, 2026-10-17)

Identity: the PDF has no title page (author-manuscript supplement, 34 pp., created
11 Jun 2019); its contents (A549 DEX, HEK293T/NIH3T3 mixtures, mouse kidney, Figs S1–S14
"sci-CAR …", legends of Tables S1–S13) and its PMC instance URL (PMC6571013) match
Cao 2018. Read: "Sample processing for cultured cells", "Mouse kidney nuclei extraction
and fixation" and "Sample processing with fixed nuclei" (text twin lines 62–310;
rendered pp. 16–18). `scrape_primers` finds four oligos, all transcribed below.

- 🟢 RT (p. 16–17): 5,000 nuclei (2 µL) + 1 µL 25 µM anchored oligo-dT primer
  `5′-ACGACGCTCTTCCGATCTNNNNNNNN[10bp index]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′` (IDT) +
  0.25 µL 10 mM dNTP, 55 °C 5 min, ice; 1.75 µL mix (5X SSIV buffer, DTT, SuperScript IV,
  RNaseOUT); 4, 10, 20, 30, 40, 50 °C 2 min each, 55 °C 10 min.
- 🟢 Tagmentation in the same well: 3 µL lysis buffer without IGEPAL, 10 µL Nextera TD
  buffer, 1–2 µL "indexed TDE1 enzyme (Illumina)", 55 °C 30 min; stopped with 20 µL 2X stop
  solution (40 mM EDTA, 1 mM spermidine). Transposon bases are not printed; the article's
  data-availability line attributes the transposase complexes to Amini 2014 (ref 7).
- 🟢 Pool, DAPI 3 µM, singlet-gated FACS of 25 nuclei per well into 5 µL EB. Second strand:
  0.66 µL NEB mRNA Second Strand Synthesis buffer + 0.34 µL enzyme, 16 °C 180 min. Stop with
  6 µL Zymo DNA binding buffer (5 min RT), purify with 24 µL AMPure XP, elute 12.5 µL EB,
  split 6 µL into each of two plates.
- 🟢 RNA plate: 5 µL Nextera TD buffer + 1 µL "i7 only TDE1 enzyme (25 nM, Illumina)",
  55 °C 5 min; 0.4 µL 1% SDS, 0.4 µL BSA, 2 µL 10 µM P7 primer
  `5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′`, 55 °C 15 min; 2 µL 10% Tween-20,
  1.2 µL water, 2 µL 10 µM indexed P5 primer
  `5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`, 20 µL NEBNext
  HiFi 2X; 72 °C 5 min, 98 °C 30 s, 18–22 cycles (98 °C 10 s, 66 °C 30 s, 72 °C 1 min),
  72 °C 5 min. Pool, 0.8 vol AMPure XP, Qubit, 6% TBE-PAGE. NextSeq 500 V2 75-cycle kit:
  Read 1 18, Read 2 52, Index 1 10, Index 2 10 cycles.
- 🟢 ATAC plate: 1 µL 10 µM indexed P5 primer
  `5′-AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC-3′`, 1 µL 10 µM P7 primer (as above),
  2 µL water, 10 µL NEBNext HiFi 2X; 72 °C 3 min, 98 °C 30 s, 18–22 cycles (98 °C 10 s,
  63 °C 30 s, 72 °C 1 min), 72 °C 5 min. Pool, 1 vol AMPure XP, gel extraction to remove
  amplicons < 200 bp, 6% TBE-PAGE. NextSeq 500 mid-output 300-cycle kit, "custom recipe
  (paired end 51 bp reads with index reads that covered both the Tn5 barcode and the library
  amplification barcode) and custom primers as previously described in (29)"; ref 29 is
  Cusanovich et al. 2018 *Nature* 555:538 (PMC5866720), not read in this pass.
- 🟢 Fixed kidney nuclei (p. 18–19): 4% PFA 15 min on ice; thawed nuclei permeabilised with
  0.2% Triton X-100 3 min; after sorting, 6 µL EB + 0.5 µL 1% SDS + 0.5 µL proteinase K,
  65 °C 16 h; 2 µL 10% Tween-20; 6 µL to an ATAC PCR plate (2 µL each primer, Tween-20,
  BSA, 20 µL NEBNext); the RNA lysate is purified with 2 vol AMPure XP, eluted in 5 µL, and
  only then gets second-strand synthesis. The model draws the cultured-cell route; this
  order difference is stated on the page as an undrawn variant.
- 🟡 Model corrections: graph sources moved from GEO to the published Supplementary
  Materials; RT, tagmentation, FACS, second-strand and cDNA-tagmentation inputs now carry
  the printed conditions; the misleading "lyse and split" step is now the printed DNA
  binding buffer + AMPure cleanup (with its supernatant as a physically removed waste) and
  6 µL + 6 µL split; both PCR steps carry their printed primers, polymerase, programme and
  cleanup (0.8× AMPure for RNA; 1× AMPure plus > 200 bp gel extraction for ATAC). Three
  discard reasons that implied physical removal of tagmentation/collision products now say
  that those products are not physically separated and that their fate is computed; the
  DAPI singlet gate is the only represented physical exclusion at FACS.
- Still open: barcode lists for the RT, Tn5 and PCR plates (not in the inspected SM or
  workbook); the "i7 only TDE1" adaptor bases (ME assumed); ATAC transposon interior
  (Amini architecture, inferred); ATAC custom sequencing primers and index-read cycles
  (ref 29, Cusanovich 2018, not read); PCR-index lengths (drawn 10 nt, consistent with the
  RNA run's 10 + 10 index cycles).


## Critical re-review, 2026-10-17

DNA and RNA tagmentation end-class candidates remain unresolved. Missing intended library ends does not establish failure to amplify or a physical separation.

Inspected scope: Local tagmentation-to-PCR declarations and terminal-fate audit; no new reaction-source assertion. This bounded correction does not complete the linked-source or two-reader absence audit.
