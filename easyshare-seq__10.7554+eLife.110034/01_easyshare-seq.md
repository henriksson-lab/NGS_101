# easySHARE-seq — source intake

Defining source: Soltys et al., *eLife* (2026, version of record),
[10.7554/eLife.110034](https://doi.org/10.7554/eLife.110034), PMC13384510, and its
Supplementary file 1. The directory is named for this published DOI; the bioRxiv preprint
[10.1101/2024.02.26.581705](https://doi.org/10.1101/2024.02.26.581705) (first posting) is
superseded and retained in the source bundle only as history. The earlier "Ma et al."
attribution in this note was wrong. Supplementary
file 1 is byte-identical (md5) to the bioRxiv `media-6` workbook, so its sheet citations below
apply to the published supplement.

🟢 RNA and accessibility outputs share two 7-nt ligation barcodes and a third PCR index.

## Published Methods read (bounded pass 2026-10-17; PMC XML text lines)

- 🟢 Key resources table (269–303); barcoding oligos (404–412); Tn5 preparation (414–418);
  fixation (420–422); tagmentation (424–430); RT (432–436); single-cell barcoding (438–442);
  reverse crosslinking (444–448); streptavidin pull-down (450–452); scATAC, cDNA and scRNA
  library preparation (454–466); sequencing (468–470).
- 🟢 Order: fixed nuclei are **tagmented first** (8 × 10,000 cells, 37 °C 30 min) with Tn5
  loaded from pooled Tn5-A/Tn5-reverse and Tn5-B2S/Tn5-reverse duplexes, then three
  tagmentation reactions are combined per in-nucleus RT (biotin-dT RT primer, Maxima H, PEG).
  The earlier model had RT before tagmentation and a B2S-only transposome.
- 🟢 Round 1 (25 °C 30 min) and round 2 (25 °C 45 min) split ligations with T4 DNA ligase,
  pooled in between; 3,500-cell sub-libraries; SDS/proteinase K reversal; M280 streptavidin
  separates biotinylated cDNA (beads) from transposed DNA (supernatant).
- 🟢 ATAC: MinElute, then Q5 PCR with i7-Tru-Seq-long and one N5XX (72 °C 6 min, 15–17
  cycles), 0.55×/1.4× bead selection. RNA: on-bead template switching (TSO, 25 °C 30 min,
  42 °C 90 min), cDNA PCR (13–15 cycles), Tn5-A-only tagmentation of 25 ng (55 °C 10 min,
  SDS stop), 5-cycle PCR with the same N5XX and i7-Tru-Seq-long.
- 🟢 Appendix 1 (PMC XML 2162–2164) states that the barcoding scheme breaks the Read 2
  primer site between the barcoding oligos and the Tn5/RT primer, and that the initial
  overhang is 8 bp long.
- 🟡 Derived from the oligo sequences (the source describes linkers but draws no duplex):
  ME-reverse covers only the ME, leaving the 5′-phosphorylated `GGGCTCGG` of Tn5-B2S or the RT
  primer single-stranded; the round-1 scaffold's 5′ `CCGAGCCC` pairs with it and its
  `ACGAGACGCTAACC` with the barcode oligo, so the barcode oligo 3′ end ligates to that 5′
  phosphate. The scaffold's 3′ `TCT` pairs with the round-2 oligo's 3′ `AGA`. The 72 °C
  starts are taken as Tn5 gap fill. Only ATAC fragments with one Tn5-A and one barcoded
  B2S end have both PCR sites; A/A and B2S/B2S fragments are candidate off-products and are
  not drawn as physically removed.
- 🟢 The RT primer's `T/iBiodT/T24` is 26 thymidines including the biotin-dT, so the RNA
  library carries A26 next to the UMI.

## Published discrepancy (shown on the page)

🟢 The eLife Key resources table (PMC XML 286–289; eLife v1 XML 257) gives Tn5-B2S as
`CAGGGCATTTGCACCCCATGCC`, whereas Supplementary file 1 (sheet "PCR oligos") gives
`/5Phos/GGGCTCGGAGATGTGTATAAGAGACAG`. The table value has no mosaic end and does not pair
with the round-1 scaffold; the model uses the workbook value. Its origin is unresolved.

## Oligo workbook

🟢 Supplementary file 1 (on disk as the md5-identical `media-6` workbook), sheets “BC oligos R1” and “BC oligos R2”,
prints each barcode together with its scaffold and explicitly says that annealing leaves
single-stranded ligation overhangs. Representative Round1_001 is
`/5Phos/GGTTAGCGTCTCGT` with scaffold
`CCGAGCCCACGAGACGCTAACCTCT`; the terminal `TCT` is retained as bases 8–10 of Index 1.

🟢 The same workbook, sheet “PCR oligos”, gives the exact branch-finishing pair.
`i7-Tru-Seq-long` is `CAAGCAGAAGACGGCATACGAGAT`; the N5XX family supplies the
opposite primer and the eight-base third cell index. Representative N501 is
`AATGATACGGCGACCACCGAGATCTACACAGCGCTAGTCGTCGGCAGCGTC`.

🟢 Both the scATAC-seq and scRNA-seq methods use that pair. The published Sequencing section (PMC XML 468–470) gives
Read 1 150 cycles, Index 1 17 cycles, Index 2 8 cycles and Read 2 150 cycles on a
NovaSeq 6000. Index 1 is the two seven-base ligation barcodes separated by their
three-base linker; Index 2 is the N5XX PCR-well index.

🟢 The published final library is Nextera-ended, not TruSeq-ended. Its Read 2
sequencing-primer site is assembled from the `GTCTCGT` tail of the round-1 barcode
and the `GGGCTCGGAGATGTGTATAAGAGACAG` handle introduced by Tn5-B2S or the RT
primer. Together these form the canonical Nextera Read 2 primer sequence.

## Native graph migration — 2026-10-18 batch02

The source is the published eLife2026 version of record, not a preprint recipe. Supplementary file1 was read through the cached media-6 workbook: MANIFEST.tsv identifies published eLife/PMC downloads with the same md5 ec0bcfbd863f68a46d70cfac79d0e4c6. N501(RNA) and N502(ATAC) are distinct source-printed illustrative choices: the published Sequencing section explicitly uses different N5XX indices for the modalities, without assigning this example pair to a particular run. Both libraries retain the published150/17/8/150 cycle layout. Exact workbook source tests cover the physical RT biotin-dT residue,26T tail, TSO RNA residues, transposome strands, barcode duplexes and both N5 primers.

The local graph owns the named reagents, source order, physical branch decisions, primers, endpoint and read declarations. Independent molecular populations occupy separate Scene rows; subsequent RNA-template fate and complete native end geometry remain explicit model boundaries. Curator-written review targets are local and are not derived from graph contents. This is a capped migration pass, not a completed missing-information audit; no missing-information claim is accepted.

Existing barcode identities are attached to the actual native barcode/primer segments as well as their endpoint complements; orientation transforms preserve identity without treating the reported sequencing base as the physical copied base.
