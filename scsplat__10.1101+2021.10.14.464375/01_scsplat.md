# scSPLAT — research notes

Defining source: Raine *et al.*, [10.1101/2021.10.14.464375](https://doi.org/10.1101/2021.10.14.464375); published as [10.1038/s41598-022-09798-2](https://doi.org/10.1038/s41598-022-09798-2).

## Evidence used

- 🟢 Single cells undergo bisulfite conversion in 384-well format.
- 🟢 Random second-strand synthesis installs part of P5 and a six-base inline cell barcode before up to 32 wells are pooled.
- 🟢 A bulk splinted-ligation reaction installs the adapter at the other end; the cell barcode begins Read 1.
- 🟢 The published [Scientific Reports Methods](https://doi.org/10.1038/s41598-022-09798-2) print two modified SPLAT-adapter oligos: `5'AmMC6/GACGTGTGCTCTTCCGATCTNNNNNN/3'AmMo` (`ss1a`) and `5'Phos/AGATCGGAAGAGCACACGTC/3'AmMo` (printed as `ss1b`). The first has a six-base random splint; the second has a 5′ phosphate.
- 🟢 The same Methods specify the barcoded strand-synthesis family `SSO_H_X` in Supplementary Table S1 and NEBNext dual-index PCR primers. Supplementary Table S1 is the next source to inspect for every barcode-bearing primer and adapter variant.
- 🟡 The published Methods say the SPLAT adapter was annealed from `ss1a` and `ss2b`, yet print the second oligo as `ss1b`. Supplementary Table S1 likewise lists `ss1b` and no `ss2b`; the intended partner name remains uncertain. The schematic does not assert an exact adapter duplex.

## Protocol boundary

The page retains the cell-first barcode/pool/splint order that distinguishes scSPLAT from PBAT and snmC-seq2.
