# DOGMA-seq — research notes

Defining source: Mimitou et al., 2021, [doi:10.1038/s41587-021-00927-2](https://doi.org/10.1038/s41587-021-00927-2).

## Sources read

- 🟢 Main article and Methods: joint RNA, chromatin-accessibility and protein workflow.
- 🟢 PMC full-text Methods, especially “Transposition and barcoding for DOGMA-seq”; publisher supplementary protocol MOESM1 and antibody-barcode workbook MOESM3.
- 🟢 Official 10x CG000338 Rev A Multiome guide (September 2020), cited by the paper, Appendix p. 80.

## Boundary

🟢 DOGMA-seq is not an alias for ASAP-seq: it adds a gene-expression branch and produces three modality-specific libraries. The pages therefore remain separate even though one paper defines both.

🟢 The cited guide prints both bead-primer constant arms and the final ATAC/RNA
library architectures. 🟡 Variable bead barcode and UMI members are shown by
role and length rather than assigned unsupported values.

## Source-printed oligos

🟢 The DOGMA-seq Methods specify the preamplification ADT additive `CCTTGGCACCCGAGAATT*C*C` (two terminal phosphorothioate linkages), the tag-library SI-PCR forward primer `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC`, and the RPI-x reverse/index primer `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`, where `x` is a user-selected sample index. The authors use the commercial 10x Multiome ATAC + Gene Expression kit, revision CG000338 Rev A, for the other branches.

🟢 The official CG000338 Rev A Appendix (PDF p. 80, lines 2468–2478)
prints PN-2000261 RNA bead primer `CTACACGACGCTCTTCCGATCT-N16-N12-T30VN`,
the ATAC bead primer `AATGATACGGCGACCACCGAGATCTACAC-N16-CGCGTCTG`, and
both final library architectures. The shorthand `T30` here denotes the
poly(dT) tract printed in the guide; it does not specify a particular
barcode/UMI member.

🔴 Exact gel-bead barcode members and complete antibody conjugates remain
unverified. MOESM3 publishes antibody-specific barcode sequences, but that
list alone does not establish the complete conjugated oligos. BioLegend
vendor source checks remain open before an absence claim can be accepted.
