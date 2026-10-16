# Olink Explore / NGS-PEA — source notes

Defining source: Wik *et al.*, *Molecular & Cellular Proteomics* (2021),
[doi:10.1016/j.mcpro.2021.100168](https://doi.org/10.1016/j.mcpro.2021.100168).

## Evidence read

- 🟢 Full article and detailed Experimental Procedures, PMC8633680.
- 🟢 Figure 1 architecture and the PCR1, PCR2, sequencing and data-processing sections.
- 🟢 Supplementary figures and tables for controls, abundance blocks and validation.
- 🟢 Olink patent family WO2021191442A1, whose worked protocol matches the
  defining paper's abundance-block, PCR1, PCR2 and panel-pooling workflow. Its
  sequence listing publishes the P5 and P7 adapter portions of the PCR primers.
- 🟢 Two readers checked the linked Olink software landing page and current vendor
  preprocessing guide. The software page itself says its contents are password
  protected, so its archives remain inaccessible. The guide is NPX Map Preprocessing
  v6.2.2 (2026-05-19); it covers `ngs2counts` and run handling without an oligo
  table. The linked CRUK Cambridge `combine-olink` example files still need source
  review. These later materials cannot prove the 2021 probe sequences absent.

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

🟢 WO2021191442A1 SEQ ID 1 is the 20-nt P5 adapter portion
`AATGATACGGCGACCACCGA`; the patent states that the PCR1 forward common primer
contains it and the same primer is reused in PCR2. SEQ ID 2 is the 24-nt P7
adapter portion `CAAGCAGAAGACGGCATACGAGAT`; each PCR2 reverse primer contains
it plus a sample-index barcode. These are published portions, not the complete
primer sequences.

🟡 The rest of the exact standard Illumina P5/Rd1SP end shown on the schematic
comes from Illumina's adapter documentation. The combined end is 58 nt. Together
with the patent's 24-nt P7 adapter, the published 148-bp total leaves exactly 66 nt
between Rd1SP and P7. This agrees with the published 66-cycle read, which starts at
Rd1SP and covers the inline assay and sample identifiers. The individual lengths
within those 66 nt have not yet been verified from an original oligo table.

🔴 The original article, its four supplements and the checked patent do not supply
the complete probe oligos, universal PCR primers, assay-barcode whitelist,
hybridization arms or sample-index sequences in the inspected sections. This is a
provisional source gap: the software archive is password protected, the checked
preprocessing guide is from 2026, and third-party example counts have not yet
received an independent second source check.
The schematic therefore shows these regions by identity and order without inventing
bases or individual lengths. No raw-read accession was found in the article's data
availability statement; this absence claim also remains under review.

## Additional public material checked

- [Olink software downloads](https://olink.com/software/download)
- [WO2021191442A1 — Method for detecting analytes of varying abundance](https://patents.google.com/patent/WO2021191442A1/en)
- [Olink Explore preprocessing technical information](https://7074596.fs1.hubspotusercontent-na1.net/hubfs/7074596/01-User%20Manuals%20for%20website/1344-Olink-explore-preprocessing.pdf)
- [CRUK Cambridge `combine-olink`](https://github.com/crukci-bioinformatics/combine-olink)

## Sequencing

🟢 Each of four panel libraries is loaded in one NovaSeq 6000 XP lane. A single
66-cycle Read 1 covers both assay barcodes and the inline sample index. There are no
separate index reads: exact barcode combinations are counted to estimate protein
abundance after control normalization.
