# stLFR — research notes

Defining source: Wang *et al.*, [10.1101/gr.245126.118](https://doi.org/10.1101/gr.245126.118).

## Evidence used

- 🟢 stLFR means single-tube Long Fragment Read and co-barcodes subfragments originating from individual 20–300-kb DNA molecules.
- 🟢 Transposons place capture sequences along the long DNA; bead-bound adapters then transfer one of millions of barcodes to the captured subfragments.
- 🟢 The published implementation completes DNB libraries for short-read sequencing; the barcode reconnects reads computationally to their long source molecule.
- 🟢 The paper-cited detailed protocol (Cheng *et al.*, [10.1038/protex.2018.116](https://doi.org/10.1038/protex.2018.116)) links a supplementary XLSX with exact common-bead, RNA bridge, transposon, PCR, 3′-branch-ligation, BGI sequencing and Illumina-compatible oligos.
- 🟢 That workbook also publishes 1,536 exact 10-nt sequences for each of Barcode1, Barcode2 and Barcode3. The three sets give the stated approximately 3.6 billion combinations; the workbook is the authoritative whitelist and is not copied into this repository.
- 🟢 The detailed protocol states paired-end 100-base genomic sequencing and says its demultiplexer expects 42 bases of barcode plus common-adapter sequence associated with the first read.
- 🔴 The exact adapter or splint geometry used to close the amplified stLFR product into the single-stranded DNB template is not printed in the defining paper or its cited detailed protocol.

## Protocol boundary

The page describes the molecular co-barcoding workflow, not the downstream haplotyping or assembly algorithms.
