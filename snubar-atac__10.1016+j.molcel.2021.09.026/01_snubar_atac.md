# SNuBar-ATAC — research notes

Defining source: Wang *et al.*, [10.1016/j.molcel.2021.09.026](https://doi.org/10.1016/j.molcel.2021.09.026).

## Evidence used

- 🟢 Each sample receives a SNuBar oligo during tagmentation before samples are pooled for 10x single-cell ATAC.
- 🟢 The oligo contains a sequence complementary to the transposome universal tail, a sample barcode, UMI and PCR handle.
- 🟢 The sample tag is molecularly distinct from the later 10x droplet cell barcode.

## Protocol boundary

SNuBar-ATAC adds pre-pooling sample identity to the standard 10x ATAC library; the parent 10x chemistry remains separately documented.
