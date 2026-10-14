# SAMOSA-Tag — research notes

## Sources read

- 🟢 Dennis et al. 2024, DOI [10.1038/s41588-024-01748-0](https://doi.org/10.1038/s41588-024-01748-0), full text, Methods and Supplementary Note 1.

## Molecular order

🟢 Permeabilized nuclei are treated with the nonspecific adenine methyltransferase EcoGII
and SAM. 🟢 Hairpin-adapter-loaded Tn5 then fragments the marked chromatin in situ while
transferring PacBio-compatible hairpin/mosaic-end structures. 🟢 Purified DNA is gap
repaired and sequenced without PCR. The closed template preserves endogenous CpG
methylation and exogenous m6dA accessibility marks on the same molecule.

## Interpretation boundary

🟢 Supplementary Table 5 publishes 26 complete 5′-phosphorylated SMRT-A barcode
hairpins. The page uses SMRT-A_bc001 as a representative exact molecule: reverse-
complement mosaic end, reverse-complement 8-nt barcode, 45-nt PacBio hairpin loop,
barcode, and mosaic end. The shared hairpin constructor verifies that the 27-nt
barcode-plus-mosaic-end arms self-pair before rendering the nucleotide-resolved fold.

🔴 The separate PacBio sequencing-primer oligo remains proprietary. This does not make
the published hairpin sequence unknown. The page now distinguishes those two facts.
Tn5 transfer still leaves 9-nt gaps; polymerase/ligase gap repair is a separate reaction
before exonuclease selection of closed templates.
