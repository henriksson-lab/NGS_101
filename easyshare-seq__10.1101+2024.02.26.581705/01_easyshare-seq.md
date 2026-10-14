# easySHARE-seq — source intake

Defining preprint: Ma et al. (2024), DOI
[10.1101/2024.02.26.581705](https://doi.org/10.1101/2024.02.26.581705).

🟢 RNA and accessibility outputs share two 7-nt ligation barcodes and a third PCR index.
The later eLife version is an additional source, not a second protocol entry.

🟢 Supplementary oligo workbook `media-6`, sheets “BC oligos R1” and “BC oligos R2”,
prints each barcode together with its scaffold and explicitly says that annealing leaves
single-stranded ligation overhangs. Representative Round1_001 is
`/5Phos/GGTTAGCGTCTCGT` with scaffold
`CCGAGCCCACGAGACGCTAACCTCT`; the terminal `TCT` is retained as bases 8–10 of Index 1.

🟢 The same workbook, sheet “PCR oligos”, gives the exact branch-finishing pair.
`i7-Tru-Seq-long` is `CAAGCAGAAGACGGCATACGAGAT`; the N5XX family supplies the
opposite primer and the eight-base third cell index. Representative N501 is
`AATGATACGGCGACCACCGAGATCTACACAGCGCTAGTCGTCGGCAGCGTC`.

🟢 Both the scATAC-seq and scRNA-seq methods use that pair. The source sequenced
Read 1 150 cycles, Index 1 17 cycles, Index 2 8 cycles and Read 2 150 cycles on a
NovaSeq 6000. Index 1 is the two seven-base ligation barcodes separated by their
three-base linker; Index 2 is the N5XX PCR-well index.

🟢 The published final library is Nextera-ended, not TruSeq-ended. Its Read 2
sequencing-primer site is assembled from the `GTCTCGT` tail of the round-1 barcode
and the `GGGCTCGGAGATGTGTATAAGAGACAG` handle introduced by Tn5-B2S or the RT
primer. Together these form the canonical Nextera Read 2 primer sequence.
