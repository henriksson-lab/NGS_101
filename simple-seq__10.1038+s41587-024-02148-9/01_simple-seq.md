# SIMPLE-seq — source intake

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. This is a notes-only entry: there is no model or schematic page yet.
> Relationships marked 🟡 *(computed)* were worked out with `lib/chemdraw.revcomp`.
> The scg_lib_structs drawn page is not a source here.

Defining paper: Bai D, Zhang X, Xiang H, Guo Z, Zhu C, Yi C. "Simultaneous
single-cell analysis of 5mC and 5hmC with SIMPLE-seq." *Nature Biotechnology*
43(1):85–96 (2025; published online 9 Feb 2024), DOI
[10.1038/s41587-024-02148-9](https://doi.org/10.1038/s41587-024-02148-9).
Not in PMC; subscription-only. The publisher PDF (printed DOI, title and authors
checked) was hand-downloaded by the maintainer on 2026-10-17 and is now the authority
for this note; the bench protocol and GEO text below are authors' documents used only
where the Methods are silent.

## Sources read

- 🟢 **Defining article, publisher PDF** (`s41587-024-02148-9.pdf`; MANIFEST row). Read in
  the `.txt`: design and overview (lines 120–185), Methods "Oligonucleotide and model DNA
  synthesis", "Oxidation and malononitrile-based selective labeling of 5hmC", "Oxidation
  and pic-borane-based selective reduction of 5mC" (lines 2582–2650), "Barcoded Tn5
  assembly" (2671–2685), "SIMPLE-seq procedures" through "Sequencing" (2692–2795) and
  "Pre-processing" / "Reads splitting" (2796–2830). Not read: Results beyond the
  overview, Extended Data figures, data-analysis sections after read splitting.
- 🟢 The publisher's Supplementary Tables workbook (MOESM2), sheet
  `Supplementary_Table1_1`, prints the library-building oligos with their chemical
  modifications. Sheets `Supplementary_Table1_2` through `_1_4` print the first,
  second and third rounds of cell-barcoding oligos. `_1_5` holds model oligos with
  5mC/5hmC controls and `_1_6` DNA-integrity qPCR primers, separate from the library.
- 🟢 The publisher's supplementary PDF (MOESM1) is a three-page reporting summary with
  no chemistry protocol.
- 🟢 The authors' code repository linked from the article,
  [github.com/cxzhu/SIMPLE-seq](https://github.com/cxzhu/SIMPLE-seq), contains a
  bench protocol `Protocol/Protocol_SIMPLE-seq_Feb2024.pdf` (8 pages, dated Feb 2024,
  uploaded 8 Jan 2025, commit `13e187f`). Its running header reads "Protocol for
  Paired-Damage-seq" although the title is "Protocol for SIMPLE-seq"; it is an authors'
  document, not the peer-reviewed Methods, and is used below only with that label.
  Its README says Read 2 carries the cell barcode and UMI, Read 1 the genomic
  sequence, and that reads are split into 5mC and 5hmC files "based on the indicator
  sequences".

- 🟢 GEO series GSE197740 (cited by PubMed ID 38336903; contributors Bai and Zhu) and its
  samples GSM5929312 (`mESC_SIMPLE`), GSM5929313 (`PBMC_SIMPLE`) and GSM7291623
  (`mouse_brain_SIMPLE`), saved as `geo_acc_<accession>_self_*.txt`. Each sample's
  `!Sample_extract_protocol_ch1` gives the authors' short protocol: Tn5 tagmentation
  with "identical barcoded assembled Tn5", two rounds of ligation-based combinatorial
  barcoding, lysis and 1× SPRI, "Second adaptor tagging and 5hmC transition and
  hmC-to-T recording", primer extension, "5mC specific labelling", final amplification
  with "Kapa HiFi uracil plus DNA polymerase" and 1.8× AMPure XP, and sequencing "using
  a MGI-2000 platform (MGI) with the following read lengths: PE 150 + 8 + 150 (Read 1 +
  Index 1 + Read 2)". This is a repository record submitted by the authors, not the
  peer-reviewed Methods.

## Oligos (Supplementary Table 1, verbatim)

`Linker-R01` `CTCGAATGCCTCAAGCACGTGGAT`; `Linker-R02`
`GGTCTGAGTTCGCACCGAAACATCGGCCAC`; `Blocker-R01`
`ATCCACGTGCTTGAGAGGCCAGAGCATTCG`; `Blocker-R02`
`GTGGCCGATGTTTCGGTGCGAACTCAGACC`. Modified oligos as printed (the workbook writes
`5Phos/` without a leading slash): `pMENTs` `5Phos/CTGTCTCTTATACACATCT/ddC/`;
`Bridge Adaptor-bottom` `AGATGTGTATAAGAGACAGCTGACGNNNNNN/-NH2`;
`Bridge Adaptor-top` `5Phos/CGTCAGCTGTCTCTTATACACATCTGACGCTGCCGACGA`;
`P7-indicator primer` `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG/caC/TGA/caC/G`
(`caC` = 5-carboxylcytosine). Final PCR: `P7 index primer`
`CAAGCAGAAGACGGCATACGAGAT[6nt i7 index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC`;
`N5 Universal primer` `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTC`.

First-round (Tn5) representative `DNA_#01` `GCATTCGAGACGCAAGATGTGTATAAGAGACAG`
(12 rows); second-round `R01_#01` `/5Phos/GTGCGAACTCAGACCAAACCGGATCCACGTGCTTGAG`
(96 rows); third-round `R02_#01`
`CAGACGTGTGCTCTTCCGATCTNNNNNNNNNNAAACCGGGTGGCCGATGTTTCG` (96 rows). One example
does not establish the oligo in any particular well; the ten random bases of `R02_#`
are a molecular identifier distinct from the three well-barcode rounds.

## Published Methods (defining article)

- 🟢 Barcoded Tn5 assembly: barcoded DNA oligos of 33 nt ("9-bp overhang for ligation
  and 5-nt barcode region and 19-nt mosaic"), printed as
  `5Phos/GCATTCGAGACGCAAGATGTGTATAAGAGACAG`, annealed to
  `5Phos/CTGTCTCTTATACACATCT/ddC` (95 °C 5 min, 1 °C/min to 4 °C); 9 µl 10 µM
  transposon + 10 µl unloaded Vazyme Tn5 (0.5 mg/ml, S601) + 17 µl coupling buffer,
  30 °C 30 min.
- 🟢 Fixation (1% formaldehyde 10 min, glycine), NIB nuclei isolation, 0.3% SDS in
  NEBuffer 2.1 at 60 °C 10 min, Triton X-100 quench at 37 °C 60 min.
- 🟢 Barcoded tagmentation: 48 wells × 10,000 nuclei with 12 well-specific barcoded
  Tn5 (four wells per barcode), 55 °C 30 min, 15 mM EDTA. The first (Tn5) barcode
  round is used as a sample-multiplexing barcode (overview, lines 181–183).
- 🟢 Ligation barcoding: BC Plate 01 with T4 DNA ligase (NEB M0202L) 37 °C 30 min,
  then R01-Blocking-Solution (Blocker-R01) 30 min; pool; BC Plate 02 with Blocker-R02;
  Termination Solution (Blocker-R02 + EDTA). ~2,000 nuclei per sub-library; NaCl/SDS/
  proteinase K 55 °C 2 h; 1× AMPure XP.
- 🟢 Second adaptor tagging: 72 °C 5 min, ice; Exonuclease I 37 °C 60 min "to remove
  excess Linker_R02"; 1.8× AMPure; 95 °C 1 min, ice; Bridge Adapter (10 µM), 7.5%
  PEG8000, 30 U T4 DNA ligase (Thermo EL0013), 20 °C 16 h; 1× AMPure.
- 🟢 5hmC: spike-in T1M/T2H/T3MH (10 pg); 1× AMPure + P-6 column; 0.05 M NaOH 37 °C
  30 min; K2RuO4 on ice 1 h; P-6; 150 mM malononitrile in 10 mM Tris pH 7.0, 37 °C
  20 h. Primer extension with the P7-indicator primer and **KAPA2G Robust HS** DNA
  polymerase in KAPA buffer GC: 95 °C 5 min, 16 × (65 °C 60 s, 68 °C 8 min), 72 °C
  10 min; Exonuclease I 37 °C 60 min; 1.8× AMPure.
- 🟢 5mC: mTET1CD oxidation (8 µM, Fe2+, α-KG, ascorbate, ATP, DTT; 37 °C 80 min);
  Qiagen Protease 50 °C 1 h; 1.8× AMPure; pic-borane (12.5 µl of 3.26 M + 2.5 µl 3 M
  sodium acetate pH 5.2 to 10 µl DNA) 70 °C 4 h; Zymo Oligo Clean & Concentrator or
  P-30 column.
- 🟢 Index PCR: 2× KAPA U+ HiFi with 2.5 µl 10 µM Universal and Index primers;
  98 °C 10 s, then 4–6 × (98 °C 10 s, 68 °C 15 s, 72 °C 1 min), 72 °C 10 min; 1.8×
  AMPure, 20 µl water.
- 🟢 Sequencing: MGISEQ2000 PE 200 + 6 + 100 (Read 1 + Index 1 + Read 2), or NovaSeq
  PE 150 + 8 + 8 + 150 (Read 1 + Index 1 + Index 2 + Read 2).
- 🟢 Read processing: UMI = first 10 bp of Read 2; barcodes at Read 2 positions 11–15
  (BC1), 48–52 (BC2) and 79–83 (BC3); reads are split into 5hmC (both indicator 5caC
  sites converted) and 5mC (unconverted) files.
- 🟢 Design: primer extension records 5hmC-to-T on the newly synthesized strand
  complementary to the original template; TET/borane then converts 5mC on the original
  template, and also turns the primer's pre-deposited 5caC into a T signal (lines
  131–135, 172–175).

### Differences between authors' documents and the published Methods

The published Methods govern; these are not treated as publication errors.

- `DNA_#` 5′ phosphate: the Methods print the barcoded Tn5 oligo with `5Phos/`; the
  workbook `DNA_#` rows print no modification. The Methods value resolves the earlier
  open question about splint ligation of `R01_#` to the `DNA_#` 5′ end; the workbook
  omission is recorded as a source inconsistency.
- Read layout: Methods MGI PE 200 + 6 + 100 or NovaSeq 150 + 8 + 8 + 150; GEO sample
  protocols state PE 150 + 8 + 150 on MGI-2000. The GEO text is an authors' deposit;
  which layout applies to each deposited run is unresolved.
- Primer extension polymerase: Methods KAPA2G Robust HS; bench protocol KAPA HiFi.
  Tagmentation plate: Methods 48 wells, 12 barcodes; bench protocol 96 wells.
  Sub-library size: Methods ~2,000 nuclei; bench protocol 5,000–10,000. Ligation naming:
  the Methods use Blocker-R01, Blocker-R02 and Linker_R02 only, matching the workbook;
  the bench protocol's Linker-R03/Quencher-R03 names are not used in the Methods.

## How they interlock (🟡 computed)

- Read 2 from the `R02_#` arm: N10 (positions 1–10), then the remainder of `R02_#`,
  `R01_#` and `DNA_#`, places the variable 5-nt windows at Read 2 positions 11–15,
  48–52 and 79–83 and the mosaic end at 84–102. That reproduces the published
  BC1/BC2/BC3 positions, and confirms the 5′→3′ order `R02_#`–`R01_#`–`DNA_#`. In
  read order BC1 is the `R02_#` (third-round) barcode and BC3 the Tn5 barcode.
- 🟡 Assuming a Read 1 primer ending at the mosaic end (the Methods do not name the
  sequencing primers), Read 1 starts with the indicator `CTGACG`, whose two `caC` are
  read-1 positions 1 and 5; converted (T) on 5hmC-strand copies, C on copies of the
  original strand carrying the unmodified Bridge Adaptor-top complement. The N6 of the
  3′-amine-blocked bottom strand is not ligated, so it is absent from the library.
- 🟡 On the MGI layout (Read 2 = 100 cycles) Read 2 ends inside the mosaic end, so
  genomic sequence comes from Read 1 there.
- `DNA_#` ends in the mosaic end `AGATGTGTATAAGAGACAG`, complementary to `pMENTs`:
  the barcoded transposon is `DNA_#` annealed to `pMENTs` (🟢 bench protocol A-3
  step 10 anneals "barcoded DNA oligos" to pMENTs and loads Tn5).
- revcomp(`Linker-R01`) = `ATCCACGTGCTTGAG`+`GCATTCGAG`: the first 15 nt pair with the
  3′ end of `R01_#`, the last 9 nt with the 5′ end of `DNA_#`.
- revcomp(`Linker-R02`) = `GTGGCCGATGTTTCG`+`GTGCGAACTCAGACC`: the first 15 nt pair
  with the 3′ end of `R02_#`, the last 15 nt with the 5′ end of `R01_#`.
- `Blocker-R02` is exactly revcomp(`Linker-R02`). `Blocker-R01` begins with the same
  15 nt as revcomp(`Linker-R01`), so it pairs with the `R01_#`-binding half of
  `Linker-R01`.
- revcomp(`Bridge Adaptor-top`) = `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTGACG`: the
  bottom strand's fixed 25 nt pair with the top strand's first 25 nt, leaving the
  bottom N6 (3′-amine-blocked) as the overhang for capturing a single-stranded 3′ end.
  The indicator primer equals this revcomp with the two C of the 3′ `CTGACG` printed as
  `caC`.
- The `P7 index primer` 3′ part `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC` covers the
  `CAGACGTGTGCTCTTCCGATCT` 5′ arm of `R02_#`; the N5 primer ends in
  `TCGTCGGCAGCGTC`, the 5′ arm of the indicator primer and bridge-top tail.

## Step order (authors' bench protocol; superseded where the Methods above differ)

🟢 Bench protocol pages 4–8 (A-3 to A-5): formaldehyde fixation, nuclei isolation,
SDS/Triton treatment; tagmentation in 96 wells with column-specific barcoded Tn5
(55 °C 30 min, EDTA quench); pooled split ligation with annealed barcode plates and
T4 DNA ligase (37 °C 30 min), Blocking solution then a second round quenched with
EDTA "Termination-Solution"; sub-libraries of 5,000–10,000 nuclei; proteinase K lysis
and 1× SPRI; 72 °C 5 min then Exonuclease I "to remove excess Linker_R02", 1.8×
AMPure; 95 °C denaturation and 3′-end ligation of the Bridge Adapter (T4 DNA ligase,
7.5% PEG8000, 20 °C 16 h); P-6 cleanup, 0.05 M NaOH denaturation, K2RuO4 oxidation
(ice, 1 h) and malononitrile labelling of 5fC (37 °C 20 h); primer extension with the
P7-indicator primer and KAPA HiFi (16 cycles of 65 °C/68 °C 8 min) to record 5hmC,
Exo I, 1.8× AMPure; mTET1CD oxidation (37 °C 80 min, optional twice), protease,
pic-borane reduction (70 °C 4 h), Zymo cleanup; final PCR with KAPA HiFi U+, Universal
and Index primers (4–6 cycles), 1.8× AMPure, qPCR quantification.

🟡 The 5′-to-3′ cell-barcode strand is `R02_#`–`R01_#`–`DNA_#` (mosaic end)–genomic
DNA–bridge-adapter copy, with Read 2 from the `R02_#` arm (N10 + barcodes) and Read 1
genomic from the Nextera-side arm; the published Read 2 barcode positions now support
this order (computed above).

## Open questions

- Which read layout (Methods MGI 200 + 6 + 100 / NovaSeq 150 + 8 + 8 + 150, or GEO
  150 + 8 + 150) applies to each deposited run; GSM7291623 also lists "Illumina HiSeq
  2500" as instrument while its protocol text says MGI-2000.
- On the NovaSeq layout, what Index 2 (8 cycles) reads: the printed `N5 Universal
  primer` has no i5 index window, and the P7 index primer carries a 6-nt i7 index read
  in 8 Index 1 cycles. Not explained in the Methods read.
- Whether all `DNA_#` rows carry the 5′ phosphate printed in the Methods (workbook
  omits it).
- The hmC-CATCH and TAPS/TET-pic-borane predecessor methods (refs. 43, 44) were not read.
- The bench protocol uses an Expi293 plasmid named both "pcDNA4-Flag-Tet1CD" and
  "pcDNA3-Flag-mTET1CD"; the Methods name pCDNA3-Flag-mTET1CD. Irrelevant to the library.
