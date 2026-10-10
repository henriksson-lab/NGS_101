# MUSIC — source notes

Defining paper: Zhao et al., *Nature* (2024), [doi:10.1038/s41586-024-07239-w](https://doi.org/10.1038/s41586-024-07239-w).

- 🟢 PMC full text, Methods and supplementary notes were read.
- 🟢 RNA linker: `5′-OH-CGAGGAGCGCTTNNNNNrArUrArGrCrArUrUrGrC-OH-3′`; the five random DNA bases are a UMI.
- 🟢 The Y-shaped DNA linker has printed top and bottom strands and a 5-nt UMI. A single-base sticky end ligates to dA-tailed genomic DNA.
- 🟢 Three sequential 14-bp barcode sets, each with 96 members, encode the cell. The complex identity combines a 16-bp 10x GEM barcode with an 8-bp i7; the adjacent 12-bp 10x UMI deduplicates molecules.
- 🟢 The run is 28-bp Read 1, 8-bp index and 150-bp Read 2. Read 2 contains cell barcodes, analyte linker and RNA or DNA insert.
