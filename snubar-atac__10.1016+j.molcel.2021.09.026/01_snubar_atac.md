# SNuBar-ATAC — research notes

Defining source: Wang *et al.*, [10.1016/j.molcel.2021.09.026](https://doi.org/10.1016/j.molcel.2021.09.026).

## Evidence used

- 🟢 Each sample is first tagmented with Illumina TDE1 and then receives an SBO oligo before nuclei are pooled for 10x single-cell ATAC (original STAR Methods, JATS lines 185–191).
- 🟢 The oligo contains a sequence complementary to the transposome universal tail, a sample barcode, UMI and PCR handle.
- 🟢 The sample tag is molecularly distinct from the later 10x droplet cell barcode.
- 🟢 The [open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC8611914/) identifies the SNuBar-ATAC barcode oligo (`SBO`) as a CST arm complementary to the Tn5 transposome universal tail, a sample barcode plus UMI, and a PCR handle. Its separate ATAC and RNA co-assay oligo (`SBOC`) adds a poly(A) capture tail; that tail is not evidence for the ATAC-only SBO.
- 🟢 The original article's STAR Methods (JATS lines 187–193) print the ATAC-only SBO family as `CCTTGGCACCCGAGAATTCCANNNNNNNN-(N)18-CTACTCTGCGTTGATACCAGACGCTGCCGACGA`: an eight-base UMI, eighteen-base sample-barcode window and source-known CST arm. They also print the TruSeq RPIX barcode-library primer `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`, P5 primer `AATGATACGGCGACCACCGAGATCTACAC`, custom Read 1 and i5 primers. The oligo family and two PCR primers are displayed on the public page.
- 🟢 After GEM cleanup the ATAC and SBO barcode templates are separated, PCR-amplified as distinct libraries, and mixed 8:2 for sequencing (JATS lines 191–193). The public final-library panel shows the ATAC branch; the reaction graph also shows the separate SBO library.
- 🔴 Table S1, which identifies individual eighteen-base sample barcodes, and detailed Methods S1 are linked but presently gated in the PMC download route. The article establishes the family architecture but does not select an exact sample barcode for the generic model. The author code repository is for demultiplexing and does not itself establish synthesis modifications.

## Protocol boundary

SNuBar-ATAC adds pre-pooling sample identity to the standard 10x ATAC library; the parent 10x chemistry remains separately documented.
