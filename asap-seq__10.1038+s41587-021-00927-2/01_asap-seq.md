# ASAP-seq — research notes

Defining source: Mimitou et al., 2021, [doi:10.1038/s41587-021-00927-2](https://doi.org/10.1038/s41587-021-00927-2).

## Sources read

- 🟢 Main article and Methods: antibody staining, bridge extension, 10x scATAC barcoding and separate ATAC/ADT amplification.
- 🟢 PMC full-text Methods, especially “Transposition and barcoding for ASAP-seq”; publisher supplementary protocol MOESM1 and TotalSeq antibody-barcode workbook MOESM3.
- 🟢 Official 10x CG000209 Rev A v1.1 guide (May 2019), Appendix p. 49; the paper specifically cites Rev D and also CG000168 Rev D for v1, which remain to compare.

## Molecular interpretation

🟢 The key ASAP-seq addition is extension of an antibody-derived tag with a sequence that can be captured in the scATAC GEM. Chromatin fragments and tags consequently inherit the same bead barcode but are recovered in separate libraries.

🟡 The public schematic uses length-preserving role tokens for reagent-generation-dependent tag and bridge sequence. It does not promote one commercial oligo layout to a universal sequence.

## Source-printed oligos

🟢 The paper prints both capture bridges, including their **3′ inverted dT block**:

| Role | 5′→3′ oligo as printed |
|---|---|
| BOA, TotalSeq-A | `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT/3InvdT/` |
| BOB, TotalSeq-B | `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTTGCTAGGACCGGCCTTAAAGC/3InvdT/` |
| RPI-x, TotalSeq-A tag index | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` |
| D7xx_s, TotalSeq-B tag index | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCAGACGTGTGC` |

🟢 `x` in the P7 primers is the user-chosen sample index; the bridge is added to the 10x ATAC barcoding mix. The paper used scATAC kit versions v1 and v1.1. MOESM3 lists antibody-specific tag barcodes.

🟢 The official CG000209 Rev A Appendix (PDF lines 1577–1607) prints the v1.1
PN-2000210 bead primer as `AATGATACGGCGACCACCGAGATCTACAC-N16-TCGTCGGCAGCGTC`,
plus the transposed-DNA and final ATAC-library architectures. It prints a
16-base barcode *position*, not every member. 🟡 The author Methods cite the
later CG000209 Rev D, so version equivalence is still to check.

🔴 The exact 10x gel-bead barcode members, complete antibody conjugates, cited
v1 CG000168 Rev D and v1.1 CG000209 Rev D guides have not yet been established
in this audit. The paper and its available supplements print the bridge bases
above; missing-information claims remain blocked.
