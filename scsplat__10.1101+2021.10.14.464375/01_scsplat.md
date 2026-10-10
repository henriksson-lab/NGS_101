# scSPLAT — research notes

Defining source: Raine *et al.*, [10.1101/2021.10.14.464375](https://doi.org/10.1101/2021.10.14.464375); published as [10.1038/s41598-022-09798-2](https://doi.org/10.1038/s41598-022-09798-2).

## Evidence used

- 🟢 Single cells undergo bisulfite conversion in 384-well format.
- 🟢 Random second-strand synthesis installs part of P5 and a six-base inline cell barcode before up to 32 wells are pooled.
- 🟢 A bulk splinted-ligation reaction installs the adapter at the other end; the cell barcode begins Read 1.

## Protocol boundary

The page retains the cell-first barcode/pool/splint order that distinguishes scSPLAT from PBAT and snmC-seq2.
