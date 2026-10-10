# SDR-seq

**Defining paper:** Lindenhofer et al., *Nature Methods* (2025), doi:10.1038/s41592-025-02805-0. The open article and detailed methods were read on the publisher site; the associated detailed protocol is doi:10.17504/protocols.io.6qpvr9q43vmk/v1.

## Molecular boundary

🟢 Fixed permeabilized cells undergo in-situ RT with custom poly(dT) primers that add a UMI, sample barcode and capture sequence to cDNA. gDNA remains in the same cell.

🟢 In Tapestri droplets, multiplex PCR copies a bead cell barcode onto targeted gDNA and cDNA amplicons through complementary capture-sequence overhangs.

🟢 Version 2 reverse primers distinguish the branches: gDNA primers carry Nextera R2N and RNA primers carry TruSeq R2, enabling separate sequencing-ready libraries.

## Evidence boundary

🟢 Target panels and sample-barcode oligos are published in supplementary tables. They vary with the chosen loci and samples; the page models the invariant UMI/barcode/overhang architecture.
