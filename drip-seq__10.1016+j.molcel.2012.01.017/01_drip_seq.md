# DRIP-seq — DNA:RNA immunoprecipitation sequencing

## Sources read

- 🟢 Ginno *et al.*, “R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters,” *Molecular Cell* (2012), DOI `10.1016/j.molcel.2012.01.017` (defining paper; read as PMC author manuscript PMC3319272 JATS text: Results line 74–77, Experimental Procedures "Genome-wide R-loop profiling using DRIP-seq and DRIVE-seq" line 113, data line 117). Its two supplements (NIHMS362480 xlsx/doc) are behind the PMC download gate and remain unread.
- 🟢 Sanz and Chédin, “High-resolution, strand-specific R-loop mapping via S9.6-based DNA:RNA ImmunoPrecipitation and high-throughput sequencing,” *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1` (authors’ detailed DRIP/DRIPc protocol).
- 🟢 The open detailed protocol (`dripc-seq__…/PMC6615061_PMC6615061.xml.txt`: applications lines 130–140, library overview 420–446, reagent list 600–621, Procedure step 50 note line 1259) was checked for its library primers and adapter naming; its 2019 TruSeq implementation is distinct from the defining 2012 library, whose adaptors are only called "Illumina barcoded adaptors".

## Molecular interpretation

- 🟢 Total nucleic acid is extracted by SDS/proteinase K at 37 °C, phenol-chloroform and ethanol precipitation, then digested with HindIII, EcoRI, BsrGI, XbaI and SspI, pre-treated or not with RNase H overnight (2012 line 113).
- 🟢 S9.6 immunoprecipitation follows the MeDIP procedure (Weber et al. 2005, not read) with the denaturation step omitted. The RNase H-pre-treated sample and input DNA are separate control libraries; peaks require >10-fold enrichment over both.
- 🟢 DRIP-seq sequences the DNA component using a non-strand-specific dsDNA library; unlike DRIPc-seq it does not preserve which genomic strand carried the RNA.
- 🟢 2012: the pulled-down material (± RNase H) and input DNA were sonicated, size-selected and ligated to Illumina barcoded adaptors, then sequenced on Illumina GAIIx and HiSeq (4.8–6 million mapped reads per library). Read mode and length are not stated there.
- 🟡 End repair, dA-tailing and library PCR are not stated in the 2012 Methods; the page draws them from the authors' 2019 DRIP-seq branch (end repair, dA-tailing, TruSeq adapter ligation, step 50 skipped, PCR step 51) and marks them inferred.
- Model correction (round 2026-10-17-breadth-first-09): an earlier version of this note said the 2012 library used Affymetrix Prep2Seq adapters and single-end sequencing. Neither claim appears in the inspected 2012 text; both were removed. Prep2Seq may come from the unread supplements and is left as an open question.
- 🟢 The later protocol names Illumina TruSeq Single index adapters and prints its PCR primers P5 `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′` and P7 `5′-CAAGCAGAAGACGGCATACGAGAT-3′`. They document the later TruSeq branch and must not be substituted for the 2012 Prep2Seq ends.
- 🔴 The inspected 2012 text does not print adapter or sequencing-primer bases (supplements unread, so this is not confirmed absence). Its page therefore keeps those ends by role; the 2019 PCR pair does not resolve the historical library.

## Final readout

Each read locates an enriched DNA fragment that contained an S9.6-bound RNA:DNA hybrid. Resolution is limited by restriction-fragment and post-capture shearing geometry, and the read itself is not strand-specific.

## Open questions

- 2012 supplements (NIHMS362480-supplement-01.xlsx, -02.doc): adapter identity, size window, read mode.
- MeDIP predecessor (Weber et al. 2005) IP and wash conditions.
