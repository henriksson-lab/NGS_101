# R-ChIP — in-vivo R-loop capture with inactive RNASEH1

## Sources read

- 🟢 Chen *et al.*, “R-ChIP Using Inactive RNase H Reveals Dynamic Coupling of R-loops with Transcriptional Pausing at Gene Promoters,” *Molecular Cell* (2017), DOI `10.1016/j.molcel.2017.10.008` (defining paper; PMCID: PMC5957070).
- 🟢 Chen *et al.*, “R-ChIP for genome-wide mapping of R-loops by using catalytically inactive RNASEH1,” *Nature Protocols* (2019), DOI `10.1038/s41596-019-0154-6` (stepwise protocol; PMCID: PMC6604627).
- 🟢 The detailed protocol's Table 1 and steps 69–74 were checked in its open XML text (`PMC6604627_PMC6604627.xml.txt`, lines 1475–1485 and 657–740); the adapter strands and PCR primers are printed there, superseding the older role-only adapter note.

## Molecular interpretation

- 🟢 Cells stably express V5-tagged catalytic mutant RNASEH1 D210N; the binding-defective WKKD mutant supplies the negative control. Formaldehyde fixes the protein–hybrid complex before sonication and anti-V5 ChIP.
- 🟢 Captured DNA is denatured and random-primed with `5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′`.
- 🟢 Extension makes dsDNA and leaves a 5′ overhang at the tailed-primer end. After dA tailing, the Illumina adapter ligates only to the opposite end; PCR therefore preserves strand orientation.
- 🟢 The detailed protocol gives the annealed adapter strands: oligo A `5′-/Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3′` and oligo B `5′-AGACGTGTGCTCTTCCGATCT-3′`. Reverse complementing oligo B gives `AGATCGGAAGAGCACACGTCT`; its 12-base `GATCGGAAGAGC` span pairs with oligo A, leaving a one-base 3′-T overhang on oligo B. The 5′ phosphate on A is part of the ligating reagent.
- 🟢 The same table gives PCR primer `5′-CAAGCAGAAGACGGCATACGAG-3′` and barcoded primer `5′-AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`; `NNNNN` is the flexible sample barcode.
- 🟢 The schematic now draws the source-written adapter heteroduplex, its 3′-T ligating end and the dA-tailed substrate as checked molecular geometry. 🟡 Its sealed library still uses a role-level adapter-body segment; the full PCR-product reconstruction remains open. The sequencer's primer supplied with the Illumina platform is not specified as a separate synthesis oligo in this protocol.

## Final readout

The sequenced fragment is the DNA component of an in-vivo RNASEH1-bound R-loop. It is a strand-specific, one-ended library and is not interchangeable with S9.6 DRIP-seq or RNA-component DRIPc-seq.
