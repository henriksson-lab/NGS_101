# snapTotal-seq — source intake

Defining paper: Salmen et al., *Nature Communications* (2024), DOI
[10.1038/s41467-024-50291-3](https://doi.org/10.1038/s41467-024-50291-3).

🟢 Low-temperature MALBAC priming and template switching amplify total RNA. An 8-nt
cell barcode is introduced after cDNA PCR by double-strand conversion, followed by a
fragmentation, enrichment, dA-tailing and Y-adapter workflow.

## Published oligo inventory

🟢 Niu et al. Supplementary Data 12 prints every named library oligo. The three
RT/MALBAC primers share `GTGAGTGATGGTTGAGGATGTGTGGAGNNNNN` and end in
`GGG`, `TTT` or `TTTTTTTTTTTT`; the five random bases are the UMI. `GAT27`
(`GTGAGTGATGGTTGAGGATGTGTGGAG`) amplifies the cDNA. The cell-barcode
conversion primer is
`GTGTGCTCTTCCGATCTNNNNNNNNAGGAGAGTGTGAGTGATGGTTGAGGATGTGTGGAG`.
`Adapter-DSC` is `GCACGACATCTGCTAACGCAGTAGTGTGCTCTTCCGATCT`.

🟢 The same table prints the indexed P5 primer
`AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT`
and the Y-adapter complement as
`P*GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG`;
`P*` denotes the 5′ phosphate as written in the table. Its final custom-primer
variant uses `Adapter-P7`
`CAAGCAGAAGACGGCATACGAGATGCACGACATCTGCTAACGCAGTA`, whereas the
standard-primer alternative uses a TruSeq-P7 primer with an eight-base i7
index. The published custom Read 2 primer is
`AGAGGTGAGTGAGTGATGGTTGAGGATGTGTGGAG` and custom i7 primer is
`CTCCACACATCCTCAACCATCACTCACTCACCTCT`.

🟢 The Methods explicitly distinguish the two read layouts: the custom-primer
run reads the UMI at the start of Read 2; the standard-primer alternative reads
the eight-base cell barcode at Read 2 bases 1–8 and the UMI at bases 44–48.
The current schematic shows that standard-primer alternative.
