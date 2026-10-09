# Nextera Mate Pair

## 1. Protocol boundary

This page follows Illumina's Nextera Mate Pair Library Prep kit, document 15035209 v02.
The defining topology is a long tagmented fragment that is circularized, re-sheared and selected
for the marked old circle junction before a conventional single-index TruSeq library is made.

## 2. Sources read

🟢 Illumina, *Nextera Mate Pair Library Prep Reference Guide*, document 15035209 v02 (2016).

🟢 Illumina, *Data Processing of Nextera Mate Pair Reads on Illumina Sequencing Platforms*.
This technical note prints the junction and external adapter sequences and explains RF mapping.

## 3. Molecular state changes

🟢 Mate Pair Tagment Enzyme fragments and tags long genomic DNA. Strand-displacement fill repairs
the staggered transposition products, and the guide then uses blunt intramolecular ligation to
make circles. Exonuclease removes remaining linear DNA.

🟢 Covaris shearing converts each large circle into short fragments of approximately 300–1,000 bp.
Streptavidin beads retain the biotin-marked fragments that span the former circle-closing bond.

🟢 The common duplicated junction printed by Illumina is
`CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG`, derived in the model as `ME_RC + ME`.
Single junctions in either orientation are also documented and can occur less frequently.

🟢 Captured fragments remain bead-bound through end repair, dA-tailing, ligation of one of twelve
TruSeq DNA LT adapters and PCR. Each adapter carries a six-base index, and the guide explicitly
requires a single-index TruSeq sequencing workflow.

## 4. Read layout

Read 1 and Read 2 enter genomic sequence from opposite ends of the short captured fragment. When
the fragment spans the junction, the genomic pieces on either side were the distant ends of the
original long molecule and align in outward-facing RF orientation. A sufficiently long read may
cross the internal junction; that sequence is trimmed during analysis.
