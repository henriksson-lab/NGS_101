# scONE-seq

**Defining paper:** Yu et al., *Science Advances* (2023), doi:10.1126/sciadv.abp8901. The open article and supplement were read from PMCID PMC9812385.

## Molecular boundary

🟢 Custom Tn5 adapters label genomic fragments with a six-base DNA-branch barcode and eight-base UMI. Internally primed RT labels total-RNA cDNA with a distinct six-base RNA-branch barcode and eight-base UMI.

🟢 Poly(C) tailing installs the opposite cDNA amplification end. DNA and RNA products are then amplified together in one tube and separated after sequencing by their branch identifiers.

## Evidence boundary

🟢 Supplementary Table S2 prints the exact adapter and RT oligos. These include
the biotinylated One-Tn5 adapter
`GTCTCGTGGGCTCGGTCATGNNNNNNNNAGATGTGTATAAGAGACAG`, the P5-Tn5
adapter `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, phosphorylated Mosaic
oligo `CTGTCTCTTATACACATCT`, three biotinylated RNA-priming variants,
the 3′ adapter, and `GATGTGTGGAGGTCTCGTGGGCTCGG` amplification primer.

🟢 The same original table publishes custom Illumina sequencing primers:
ONE-Illu Read 1 `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, ONE-Illu Read 2
`CCGAGCCCACGAGACCTGTCTCTTATACACATCT`, and ONE-Illu i7
`AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG`. It also prints the full P5 and
P7 index-primer families with symbolic `[i5]` and `[i7]` windows.

🟢 The ONE-Illu i5 index-primer row is marked with an asterisk directing users
to the sequencing-equipment manufacturer's instructions and does not print a distinct
sequence. That is the remaining source boundary; it does not justify hiding the other
published primers.
