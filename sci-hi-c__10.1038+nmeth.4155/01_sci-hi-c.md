# sci-Hi-C

## 1. What it is

Ramani et al.'s combinatorial-indexing Hi-C method, in which the intersection of two plate barcodes identifies a nucleus without physical isolation into one well at the start.

## 2. Sources and evidence

- 🟢 Ramani et al., *Nature Methods* 2017, DOI [10.1038/nmeth.4155](https://doi.org/10.1038/nmeth.4155); full text [PMC5330809](https://pmc.ncbi.nlm.nih.gov/articles/PMC5330809/).
- 🟢 DpnII-digested nuclei receive biotinylated indexed bridge adaptors in the first 96-well plate.
- 🟢 Nuclei are pooled, proximity-ligated, redistributed, and receive a second index through custom Illumina Y adaptors before streptavidin enrichment and PCR.
- 🟡 The page preserves the published 8-base first-round barcode length but uses placeholders for barcode identity and unresolved bridge bases.

## 3. Distinguishing chemistry

The first barcode is internal to the contact molecule. The second barcode resides in the sequencing adapter. Both are required to assign read pairs to a nucleus, so the two barcode rounds are never collapsed into a generic dual-index library.
