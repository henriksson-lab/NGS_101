# CRISPR CIRCLE-seq — research notes

## Sources read

- 🟢 Tsai et al. 2017, DOI [10.1038/nmeth.4278](https://doi.org/10.1038/nmeth.4278), PMC5924695, main text, Methods and supplementary protocol.

## Molecular order

🟢 Genomic DNA is sheared to about 300 bp, repaired, dA-tailed and ligated to the
uracil-containing oSQT1288 stem-loop. Exonucleases remove open molecules; USER plus PNK
opens the loop and dilute T4 ligation circularizes individual fragments. Plasmid-Safe
DNase removes residual linear DNA. Cas9–guide complexes linearize target-bearing circles.
The fresh ends are dA-tailed, ligated to a USER-openable hairpin adapter, opened, PCR
amplified and sequenced 2×150.

## Disambiguation

This method creates circles from ordinary genomic fragments to select Cas9 cleavage.
eccDNA Circle-Seq starts with naturally circular molecules and uses exonuclease selection
plus phi29 amplification. The shared word “Circle” does not make them one protocol.
