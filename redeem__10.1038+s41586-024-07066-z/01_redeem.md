# ReDeeM — research notes

Defining source: Weng *et al.*, [10.1038/s41586-024-07066-z](https://doi.org/10.1038/s41586-024-07066-z).

## Evidence used

- 🟢 ReDeeM modifies droplet-based 10x Multiome using whole cells and conditions that retain mitochondrial DNA while preserving RNA and ATAC libraries.
- 🟢 Tiled mtDNA probes hybrid-capture a third, deeply sequenced library carrying the same cell barcodes.
- 🟢 Cell barcode plus mtDNA fragment endpoints act as endogenous molecular identifiers for consensus error correction.
- 🟢 The article's [Supplementary Data 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10937407/) is an XLSX with a `Primers` sheet and eight named mtDNA-probe sheets (four human and four mouse). It publishes the individual capture-probe sequences and the `Nextera_N701`–`Nextera_N706`, P5, P7 and `Hash-Spikin primer` sequences. The full probe inventory belongs in that source workbook, not a copied table here.
- 🟢 The separate supplementary protocol specifies four staggered probe pools for each mtDNA capture, post-capture PCR with P5/P7, and paired-end sequencing of the mtDNA branch (R1 150, i7 8, i5 24, R2 150 cycles). It also prints optional HTO-additive and Hash-primer oligos and a generic `N701`–`N712` primer pattern.
- 🟡 The supplementary protocol predates the final workbook: it says probe sequences were forthcoming and describes four sets of 137 100-nt probes. The published workbook instead contains 120-nt probe sequences (with a few shorter rows) and different row counts across its human sheets. Use the workbook for sequence identities; its plain base strings do not by themselves establish terminal capture-probe modifications.
- 🔴 The underlying 10x Multiome bead oligos and proprietary `SI PCR primer B` are identified by reagent name, not printed in the ReDeeM supplement. Their full molecular structures require the corresponding 10x source, rather than inference from the capture-probe table.

## Protocol boundary

ReDeeM captures mtDNA from Multiome material; MAESTER instead enriches mitochondrial transcripts from scRNA-seq cDNA.
