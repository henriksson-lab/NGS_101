# SDR-seq

**Defining paper:** Lindenhofer et al., *Nature Methods* (2025), doi:10.1038/s41592-025-02805-0. The open article and detailed methods were read on the publisher site; the associated detailed protocol is doi:10.17504/protocols.io.6qpvr9q43vmk/v1.

## Molecular boundary

🟢 Fixed permeabilized cells undergo in-situ RT with custom poly(dT) primers that add a UMI, sample barcode and capture sequence to cDNA. gDNA remains in the same cell.

🟢 In Tapestri droplets, multiplex PCR copies a bead cell barcode onto targeted gDNA and cDNA amplicons through complementary capture-sequence overhangs.

🟢 Version 2 reverse primers distinguish the branches: gDNA primers carry Nextera R2N and RNA primers carry TruSeq R2, enabling separate sequencing-ready libraries.

## Evidence boundary

🟢 Target panels and sample-barcode oligos are published in supplementary tables. They vary with the chosen loci and samples; the page models the invariant UMI/barcode/overhang architecture.

## Oligo source audit

🟢 Supplementary Table 1 (proof-of-concept v1) and Supplementary Table 2 (v2) are in `s41592-025-02805-0_supp_41592_2025_2805_MOESM3_ESM.xlsx`, sheets `Supplementary_Table_1` and `Supplementary_Table_2`. These tables give actual primer sequences, not only an assay diagram. The invariant v2 droplet-amplification forward primer is `GTACTCGCAGTAGTCGACACGTC`; the v2 RNA reverse primer is `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-cDNA_specfic` as printed, where the locus-specific suffix is a placeholder. The corresponding v2 gDNA primers are `GTACTCGCAGTAGTC-gDNA_specific` and `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-gDNA_specific`.

🟢 The same table prints complete v2 RT-primer examples, including `RT_primer_1_8` as `GTACTCGCAGTAGTCGACACGTCTCGCCTTANNNNNNNNTTTTTTTTTTTTTTTTVN`. Its sample-barcode choice changes across the table. Separate i5 primers and branch-specific i7 primers are listed there as full sequences; these are choices from the published index sets, not a single immutable library end.

🟡 The diagram's role-level representation remains appropriate for the variable target panel and barcode/index choices. The source does resolve the common primer handles above; any exact read-placement model should select a stated v1/v2 and sample-index instance rather than treating the variable rows as one molecule.
