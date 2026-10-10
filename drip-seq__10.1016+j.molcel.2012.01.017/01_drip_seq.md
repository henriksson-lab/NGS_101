# DRIP-seq — DNA:RNA immunoprecipitation sequencing

## Sources read

- 🟢 Ginno *et al.*, “R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters,” *Molecular Cell* (2012), DOI `10.1016/j.molcel.2012.01.017` (defining paper).
- 🟢 Sanz and Chédin, “High-resolution, strand-specific R-loop mapping via S9.6-based DNA:RNA ImmunoPrecipitation and high-throughput sequencing,” *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1` (authors’ detailed DRIP/DRIPc protocol).

## Molecular interpretation

- 🟢 Genomic nucleic acid is extracted gently and digested with HindIII, SspI, EcoRI, BsrGI and XbaI so RNA:DNA hybrids remain intact.
- 🟢 S9.6 immunoprecipitates the hybrid. An RNase H-treated aliquot is the essential specificity control.
- 🟢 DRIP-seq sequences the DNA component using a non-strand-specific dsDNA library; unlike DRIPc-seq it does not preserve which genomic strand carried the RNA.
- 🟢 The captured material is sonicated, end-repaired, dA-tailed, adapter-ligated and PCR-amplified. The original implementation used Affymetrix Prep2Seq adapters and single-end Illumina sequencing.
- 🔴 Complete adapter and sequencing-primer oligo sequences are not printed in the defining paper, so the schematic names their roles instead of inventing bases.

## Final readout

One read locates an enriched DNA fragment that contained an S9.6-bound RNA:DNA hybrid. Resolution is limited by restriction-fragment and post-capture shearing geometry, and the read itself is not strand-specific.
