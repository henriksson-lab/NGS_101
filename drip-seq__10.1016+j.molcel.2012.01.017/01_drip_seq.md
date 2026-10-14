# DRIP-seq — DNA:RNA immunoprecipitation sequencing

## Sources read

- 🟢 Ginno *et al.*, “R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters,” *Molecular Cell* (2012), DOI `10.1016/j.molcel.2012.01.017` (defining paper).
- 🟢 Sanz and Chédin, “High-resolution, strand-specific R-loop mapping via S9.6-based DNA:RNA ImmunoPrecipitation and high-throughput sequencing,” *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1` (authors’ detailed DRIP/DRIPc protocol).
- 🟢 The open detailed protocol (`PMC6615061_PMC6615061.xml.txt`, reagent list lines 606–621) was checked for its library primers and adapter naming; its 2019 TruSeq implementation is distinct from the defining 2012 Affymetrix Prep2Seq library.

## Molecular interpretation

- 🟢 Genomic nucleic acid is extracted gently and digested with HindIII, SspI, EcoRI, BsrGI and XbaI so RNA:DNA hybrids remain intact.
- 🟢 S9.6 immunoprecipitates the hybrid. An RNase H-treated aliquot is the essential specificity control.
- 🟢 DRIP-seq sequences the DNA component using a non-strand-specific dsDNA library; unlike DRIPc-seq it does not preserve which genomic strand carried the RNA.
- 🟢 The captured material is sonicated, end-repaired, dA-tailed, adapter-ligated and PCR-amplified. The original implementation used Affymetrix Prep2Seq adapters and single-end Illumina sequencing.
- 🟢 The later protocol names Illumina TruSeq Single index adapters and prints its PCR primers P5 `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′` and P7 `5′-CAAGCAGAAGACGGCATACGAGAT-3′`. They document the later TruSeq branch and must not be substituted for the 2012 Prep2Seq ends.
- 🔴 The defining 2012 paper does not print the complete Affymetrix Prep2Seq adapter or sequencer primer bases. Its page therefore keeps those ends by role; the 2019 PCR pair does not resolve the historical library.

## Final readout

One read locates an enriched DNA fragment that contained an S9.6-bound RNA:DNA hybrid. Resolution is limited by restriction-fragment and post-capture shearing geometry, and the read itself is not strand-specific.
