# scRCAT-seq — source intake

Defining paper: Liu et al., *Nature Communications* (2020), DOI
[10.1038/s41467-020-18976-7](https://doi.org/10.1038/s41467-020-18976-7).

🟢 Reverse transcription tags both cDNA boundaries with a UMI and cell identifier;
tagmentation and targeted PCR enrich 5′- and 3′-end fragments.

## Exact library oligos

🟢 The open article Methods, “scRCAT-seq library construction for a single cell,”
give the RT primer as
`AAGCAGTGGTATCAACGCAGAGTN8[16 bps of cell barcode]T30VN`.
The eight random bases are the UMI; the cell barcode is a separate 16-base
field. The amplified cDNA is tagmented with Nextera XT.

🟢 Selective endpoint amplification uses the complete P5 primer
`AATGATACGGCGACCACCGAGATCTACAC[8 bps of index]TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGGTATCAACGCAGAGT`
and P7 primer
`CAAGCAGAAGACGGCATACGAGAT[8 bps of index]GTCTCGTGGGCTCGG`.
These tails fix more of the completed library than the current compact schematic
draws. The article refers to the Smart-seq2 template-switch and amplification
steps rather than printing their complete oligos here; their inherited sequence
should be checked against the cited Smart-seq2 protocol before replacing the
schematic's inferred states. The Methods specify paired-end 150 sequencing.
