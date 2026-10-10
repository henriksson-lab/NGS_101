# stLFR — research notes

Defining source: Wang *et al.*, [10.1101/gr.245126.118](https://doi.org/10.1101/gr.245126.118).

## Evidence used

- 🟢 stLFR means single-tube Long Fragment Read and co-barcodes subfragments originating from individual 20–300-kb DNA molecules.
- 🟢 Transposons place capture sequences along the long DNA; bead-bound adapters then transfer one of millions of barcodes to the captured subfragments.
- 🟢 The published implementation completes DNB libraries for short-read sequencing; the barcode reconnects reads computationally to their long source molecule.

## Protocol boundary

The page describes the molecular co-barcoding workflow, not the downstream haplotyping or assembly algorithms.
