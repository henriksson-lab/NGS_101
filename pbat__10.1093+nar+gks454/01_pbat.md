# PBAT — research notes

## 1. What it is

Post-bisulfite adaptor tagging reverses conventional WGBS order: destructive conversion
comes first, and two random-primer extensions install the sequencing handles afterward.

## 2. Sources read

- 🟢 Miura et al. 2012, DOI [10.1093/nar/gks454](https://doi.org/10.1093/nar/gks454),
  complete article and supplement route (PMC3458524).

## 3. Molecular order

🟢 Bisulfite-treated single-stranded DNA is copied from a 5′-biotinylated adapter primer
ending in random N4. The product is captured on streptavidin beads, residual primer is
removed, and alkaline denaturation exposes the immobilized first strand. 🟢 A second
adapter primer, also ending in N4, makes the strand released for sequencing and size
selection. The original demonstration supports PCR-free library production.

## 4. Evidence boundary

🟢 N4 architecture, biotin capture, extension order and directionality come from the
paper. 🟢 The main paper prints BioPEA2N4 as
`5′-biotin-ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNN-3′`, PE-reverse-N4 as
`5′-CAAGCAGAAGACGGCATACGAGATNNNN-3′`, and Primer-3 as
`5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`.
The original readout is unindexed single-end sequencing of the generally G-poor strand
complementary to bisulfite-converted DNA.
