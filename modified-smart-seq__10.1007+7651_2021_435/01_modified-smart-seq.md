# Modified SMART-Seq — research notes

Defining protocol: Zheng *et al.* (2022), DOI [10.1007/7651_2021_435](https://doi.org/10.1007/7651_2021_435). Author-hosted full text was read.

## Evidence

- 🟢 Table 1 gives the barcoded RT primer as `CTACACGACGCTCTTCCGATCT[8-bp well barcode]NNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTVN` and the TSO as `AAGCAGTGGTATCAACGCAGAGTACATrGrGrG`.
- 🟢 Table 1 also prints both cDNA PCR primers (`CTACACGACGCTCTTCCGATCT` and `AAGCAGTGGTATCAACGCAGAG`) and the final library PCR pair: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC` and `CAAGCAGAAGACGGCATACGAGAT[8-bp plate barcode]GTCTCGTGGGCTCGG`. The named plate barcode is a selected eight-base index, not a random UMI.
- 🟢 The methods specify reverse transcription and template switching in each well, cDNA PCR, Exonuclease I cleanup, pooling, Tn5 tagmentation and Nextera-style index PCR.
- 🟢 Read 1 contains the 8-bp well barcode followed by the 10-bp UMI; Index 1 is the 8-bp plate index; Read 2 contains transcript sequence.

This is a named 3′-counting modification, not a replacement for the separate SMART-seq, SMART-seq2 or SMART-seq3 pages.
