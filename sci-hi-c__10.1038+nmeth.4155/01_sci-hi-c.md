# sci-Hi-C

## 1. What it is

Ramani et al.'s combinatorial-indexing Hi-C method, in which the intersection of two plate barcodes identifies a nucleus without physical isolation into one well at the start.

## 2. Sources and evidence

- 🟢 Ramani et al., *Nature Methods* 2017, DOI [10.1038/nmeth.4155](https://doi.org/10.1038/nmeth.4155); full text [PMC5330809](https://pmc.ncbi.nlm.nih.gov/articles/PMC5330809/).
- 🟢 DpnII-digested nuclei receive biotinylated indexed bridge adaptors in the first 96-well plate.
- 🟢 Nuclei are pooled, proximity-ligated, redistributed, and receive a second index through custom Illumina Y adaptors before streptavidin enrichment and PCR.
- 🟢 Supplementary Data print the bridge-adapter oligos and all 96 Y-adapter pairs. The
  model uses A1 representatives for both rounds: BC1 `GAGCTCGA` in the
  `GATC–BC1–GAATT` biotinylated bridge strand, and BC2 `TGACCTTG` in the SeqShort HI/LO
  Y adapter. The published internal biotin-dT, 5′ phosphate, terminal phosphorothioate,
  DpnII-compatible end, EcoRI overhang and dA/dT ligation end are retained.

## 3. Distinguishing chemistry

The first barcode is internal to the contact molecule. BC2 is the first eight bases of
both genomic reads (followed by three fixed bases), not an Illumina index read. Both
barcode rounds are required, so they are never collapsed into a generic dual-index library.
