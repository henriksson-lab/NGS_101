# DOGMA-seq — research notes

Defining source: Mimitou et al., 2021, [doi:10.1038/s41587-021-00927-2](https://doi.org/10.1038/s41587-021-00927-2).

## Sources read

- 🟢 Main article and Methods: joint RNA, chromatin-accessibility and protein workflow.
- 🟢 PMC full-text Methods, especially “Transposition and barcoding for DOGMA-seq”; publisher supplementary protocol MOESM1 and antibody-barcode workbook MOESM3.

## Boundary

🟢 DOGMA-seq is not an alias for ASAP-seq: it adds a gene-expression branch and produces three modality-specific libraries. The pages therefore remain separate even though one paper defines both.

🟡 Proprietary bead regions are shown by role and length rather than assigned unsupported bases.

## Source-printed oligos

🟢 The DOGMA-seq Methods specify the preamplification ADT additive `CCTTGGCACCCGAGAATT*C*C` (two terminal phosphorothioate linkages), the tag-library SI-PCR forward primer `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC`, and the RPI-x reverse/index primer `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`, where `x` is a user-selected sample index. The authors use the commercial 10x Multiome ATAC + Gene Expression kit, revision CG000338 Rev A, for the other branches.

🔴 The exact gel-bead barcode oligo members and complete vendor primer architecture are not specified in the paper or downloaded supplement. MOESM3 does publish antibody-specific barcode sequences, but that list does not establish the complete conjugated oligos.
