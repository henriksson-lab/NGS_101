# SNuBar-ATAC — research notes

Defining source: Wang *et al.*, [10.1016/j.molcel.2021.09.026](https://doi.org/10.1016/j.molcel.2021.09.026).

## Evidence used

- 🟢 Each sample receives a SNuBar oligo during tagmentation before samples are pooled for 10x single-cell ATAC.
- 🟢 The oligo contains a sequence complementary to the transposome universal tail, a sample barcode, UMI and PCR handle.
- 🟢 The sample tag is molecularly distinct from the later 10x droplet cell barcode.
- 🟢 The [open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC8611914/) identifies the SNuBar-ATAC barcode oligo (`SBO`) as a CST arm complementary to the Tn5 transposome universal tail, a sample barcode plus UMI, and a PCR handle. Its separate ATAC and RNA co-assay oligo (`SBOC`) adds a poly(A) capture tail; that tail is not evidence for the ATAC-only SBO.
- 🔴 The article's Methods S1 and supplementary figures need a full sequence-level check before any exact SBO, CST arm, barcode inventory or PCR primer can be marked source-verbatim. The author code repository is for demultiplexing and does not itself establish synthesis modifications.

## Protocol boundary

SNuBar-ATAC adds pre-pooling sample identity to the standard 10x ATAC library; the parent 10x chemistry remains separately documented.
