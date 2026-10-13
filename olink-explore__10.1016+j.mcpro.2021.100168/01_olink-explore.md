# Olink Explore / NGS-PEA — source notes

Defining source: Wik *et al.*, *Molecular & Cellular Proteomics* (2021),
[doi:10.1016/j.mcpro.2021.100168](https://doi.org/10.1016/j.mcpro.2021.100168).

## Evidence read

- 🟢 Full article and detailed Experimental Procedures, PMC8633680.
- 🟢 Figure 1 architecture and the PCR1, PCR2, sequencing and data-processing sections.
- 🟢 Supplementary figures and tables for controls, abundance blocks and validation.

## Molecular workflow

🟢 Each protein assay uses a matched pair of antibody-conjugated oligos. One arm
contains P5 and the Read 1 sequencing-primer site; the other contains a common
reverse preamplification site. Both carry assay-specific barcodes and complementary
hybridization arms. Co-binding to the same protein raises their effective local
concentration, allowing the arms to hybridize.

🟢 PCR1 begins with 20 minutes at 50 °C for polymerase extension across the paired
arms, followed by 25 cycles with universal primers. Products from four concentration
blocks are pooled for each sample. PCR2 uses the same forward primer and one of 96
reverse primers to add an inline sample index and P7 during ten cycles.

🟢 The finished amplicon is 148 bp. Its published element order is P5, Rd1SP,
forward assay barcode, hybridization-derived junction, reverse assay barcode,
sample index and P7.

🔴 The paper does not publish the probe oligos, universal PCR primers, assay barcode
whitelist, hybridization arms or sample-index sequences. The schematic therefore
shows their identities and order without assigning bases or unsupported lengths.

## Sequencing

🟢 Each of four panel libraries is loaded in one NovaSeq 6000 XP lane. A single
66-cycle Read 1 covers both assay barcodes and the inline sample index. There are no
separate index reads: exact barcode combinations are counted to estimate protein
abundance after control normalization.
