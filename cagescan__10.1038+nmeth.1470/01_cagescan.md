# CAGEscan — research notes

Plessy *et al.*, *Nature Methods* (2010),
[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470),
[original author manuscript PMC2906222](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906222/).
The earlier note's PMC2900199 citation was a repository error.

## Originals inspected

- 🟢 Actual main Results lines 111–115 and Figure 1 distinguish CAGEscan from nanoCAGE.
- 🟢 Actual Online Methods, `nmeth.1470_supp_NIHMS207200-supplement-supplement_1.pdf`:
  inherited nanoCAGE RT conditions plus CAGEscan-specific Methods, text 1–118, and paired-end
  processing 225–241. This is a separate original from the main article's Results.
- 🟢 Original Supplementary Table 8, MOESM315 PDF p.26, text 5134–5170: all eight CAGEscan
  oligos, including the custom paired-end sequencing primers. Both source text and the
  rendered original table were inspected independently twice.
- 🟢 Rendered original Online Methods printed p.17 (PDF p.3) was checked for primer units.

The oligo scraper was run before reading the full source bundle. Remaining supplements,
linked references and vendor instructions are not all dispositioned; this is a bounded audit.

## Source-owned route

Random N6 priming and preferential template switching at capped RNA use the published
TSO, with an inline three- or six-base DNA sample barcode and terminal RNA `rGrGrG`.
CAGEscan substitutes SuperScript III or PrimeScript for the nanoCAGE SSII enzyme and
cleans first-strand material using RNAClean XP.

The source then specifies two PCR stages, each with a pilot quantitative PCR:

1. A 1.2-µL pilot aliquot of the named 30-µL first-strand pool informs cycle selection;
   a 3-µL aliquot goes to Ex Taq HS semi-suppressive amplification. The model computes
   the residual pool and leaves its later handling unresolved.
2. AMPure XP cleanup is followed by an 8-ng adapter-PCR pilot and the remaining 20-ng
   production reaction. This PCR introduces the paired-end sequencing ends and replaces
   the interfering semi-suppressive reverse end.
3. A second AMPure XP cleanup precedes the source 5-of-30-µL agarose-gel check and
   concentration adjustment for Genome Analyzer II sequencing. The gel check is QC,
   not preparative gel selection.

The graph keeps cap-linked cDNA, other first-strand templates and primer artifacts as
**separate molecular classes**. Semi-suppressive PCR disfavors short/homoadapter templates;
it does not physically discard every artifact. Removed bead-cleanup fractions and QC
aliquots have separate fates. CAGEscan does not inherit nanoCAGE's EcoP15I cleavage/ligation.

## Eight published oligos and unfinished assembly

All eight source oligos remain displayed on the schematic. Six structural inputs have
local molecular states: RT handle plus N6, TSO handle plus variable inline barcode and RNA
tail, both semi-suppressive primers and both final PCR primers. The barcode has structured
sample-identifier meaning; there is no invented eight-base instrument index.

Table 8 prints `TAGTCGAACGAAG...` in the RT primer and `TAGTCGAACTGAAG...` in the
semi-suppressive reverse primer. The difference is preserved. A shorter exact 3′ match
exists, so this is not classified as a proven ordered-oligo error.

Both custom sequencing primers are source-positive. Their complete assembled landing
geometry, exact study barcode members and instrument cycle counts remain model questions.
The unassembled initial PCR-template geometry is explicitly separate from the amplified
library endpoint. This limitation is not evidence that primer sequences were unpublished.

## Apparent source concentration-unit issue

The original Online Methods printed p.17 visibly gives **100 mM** semi-suppressive primers
and **200 mM** final-PCR primers, while nearby dNTPs are printed as **200 µM**. Both readers
checked the original PDF, so the primer units are not merely an extraction artifact.

The public page displays these source locations and values as an apparent unit issue.
The model retains the published oligos and operation order without assigning a corrected
primer concentration. The intended concentration remains unresolved.

## Readout and remaining sources

One custom paired-end read enters at the cap/TSO-derived end, including the inline barcode;
the other enters from the random-primer-defined downstream end. Genome Analyzer II is
explicitly named. Paired-end alignment, rDNA alignment exclusion and collapsing identical
mapping pairs are computational operations, not physical library purification.

The review remains in progress for actual linked materials, including reference 16
(Bentley *et al.* 2008, DOI 10.1038/nature07517), cited nanoCAGE/RT methodological provenance,
manufacturer instructions for RNAClean XP/AMPure XP, and remaining supplement scope.
No unavailable “instrument-context document” is invented from the unfinished model.
No missing-information claim is accepted at this checkpoint.
