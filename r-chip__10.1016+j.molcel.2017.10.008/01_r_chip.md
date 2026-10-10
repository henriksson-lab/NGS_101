# R-ChIP — in-vivo R-loop capture with inactive RNASEH1

## Sources read

- 🟢 Chen *et al.*, “R-ChIP Using Inactive RNase H Reveals Dynamic Coupling of R-loops with Transcriptional Pausing at Gene Promoters,” *Molecular Cell* (2017), DOI `10.1016/j.molcel.2017.10.008` (defining paper; PMCID: PMC5957070).
- 🟢 Chen *et al.*, “R-ChIP for genome-wide mapping of R-loops by using catalytically inactive RNASEH1,” *Nature Protocols* (2019), DOI `10.1038/s41596-019-0154-6` (stepwise protocol; PMCID: PMC6604627).

## Molecular interpretation

- 🟢 Cells stably express V5-tagged catalytic mutant RNASEH1 D210N; the binding-defective WKKD mutant supplies the negative control. Formaldehyde fixes the protein–hybrid complex before sonication and anti-V5 ChIP.
- 🟢 Captured DNA is denatured and random-primed with `5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′`.
- 🟢 Extension makes dsDNA and leaves a 5′ overhang at the tailed-primer end. After dA tailing, the Illumina adapter ligates only to the opposite end; PCR therefore preserves strand orientation.
- 🔴 The paper calls the second oligo an Illumina adapter without printing its sequence, so that end remains a named role on the schematic.

## Final readout

The sequenced fragment is the DNA component of an in-vivo RNASEH1-bound R-loop. It is a strand-specific, one-ended library and is not interchangeable with S9.6 DRIP-seq or RNA-component DRIPc-seq.
