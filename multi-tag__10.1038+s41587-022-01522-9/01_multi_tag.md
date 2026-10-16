# MulTI-Tag

**Defining source:** doi:10.1038/s41587-022-01522-9.

## Molecular boundary

🟢 Sequential antibodies recruit pA–Tn5 complexes whose loaded adapters encode protein-target identity; single-cell indexing is added downstream.

## Evidence boundary

🟢 The target barcode is structured separately from the cell barcode. The public protocols.io method is DOI 10.17504/protocols.io.6qpvr4bezgmk/v2.

## Published oligos

🟢 The published article points to Supplementary Table 1 (PMC10188359, MOESM2),
and the linked preprint's `media-2.xlsx` gives the P5_i5 and P7_i7 pA–Tn5 adapter
families and the i5 and i7 amplification-primer families, with their eight-base
barcodes and experiment assignments. The source manifest records the published
MOESM2 workbook and the preprint `media-2.xlsx` as byte-identical (MD5
`5cd0f4d7b6de26836d8c04dac29ecb93`); the local `media-2.xlsx` therefore
also supplies the published table.
The prose of the preprint specifies the common
reverse mosaic-end oligo as 5′-[phos]CTGTCTCTTATACACATCT-3′ (methods,
“Conjugates for MulTI-Tag”). This is the phosphorylated annealing partner for the
loaded adapters.

🟢 For example, table row `P5_i5_1` assigns H3K27me3 the adapter
5′-TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG-3′,
with target barcode `TATAGCCT`. Table row `sci_Ad1.1_TAGATCGC` gives an i5
amplification primer
5′-AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCTCCACGC-3′.
The full table is necessary to represent the 72-by-72 single-cell index choices;
these examples do not identify a universal adapter or primer.

🟡 The current molecular page names these oligo roles but does not consume the
published table sequences. Its generic sequencing-primer limitation therefore
describes a modeling gap, not an absence of the adapter and amplification-primer
sequences from the source. The published paired-end run is 25 + 25 bases (preprint
methods, “Sequencing and data preprocessing”); exact primer landing still needs
reconstruction from a specific table-row choice and the final library.
