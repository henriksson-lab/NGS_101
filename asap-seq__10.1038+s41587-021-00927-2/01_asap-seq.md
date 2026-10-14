# ASAP-seq — research notes

Defining source: Mimitou et al., 2021, [doi:10.1038/s41587-021-00927-2](https://doi.org/10.1038/s41587-021-00927-2).

## Sources read

- 🟢 Main article and Methods: antibody staining, bridge extension, 10x scATAC barcoding and separate ATAC/ADT amplification.
- 🟢 PMC full-text Methods, especially “Transposition and barcoding for ASAP-seq”; publisher supplementary protocol MOESM1 and TotalSeq antibody-barcode workbook MOESM3.

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

🟢 `x` in the P7 primers is the user-chosen sample index; the bridge is added to the 10x ATAC barcoding mix. The paper used scATAC kit versions v1 and v1.1. MOESM3 lists antibody-specific tag barcodes, while the complete commercial bead oligo and every conjugated antibody-oligo structure are not disclosed.

🔴 Exact 10x gel-bead barcode members and kit-specific oligo architecture remain undisclosed in these source files. These are distinct from the explicitly published bridge bases above.
