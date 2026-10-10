# mNET-seq

**Defining paper:** Nojima et al., *Cell* (2015), doi:10.1016/j.cell.2015.03.027. The open author manuscript and its extended methods were read from PMCID PMC4410947.

## Molecular boundary

🟢 Native chromatin is solubilized by MNase, and elongation complexes are immunoprecipitated with antibodies recognizing Pol II or particular CTD phosphorylation states. RNA protected in those complexes is purified and size-selected.

🟢 The informative nucleotide is the RNA 3′ end. It may be the nucleotide in the Pol II active site or a cotranscriptional cleavage intermediate, so the page does not equate every peak with polymerase position.

🟢 RNA 5′ ends are phosphorylated before small-RNA adapter ligation. The original implementation used the Illumina TruSeq Small RNA workflow after end repair.

## Evidence boundary

🟡 The paper establishes adapter roles and read direction but does not make one complete adapter/index oligo set an invariant of mNET-seq. Those regions are therefore role-labelled rather than reconstructed from unrelated modern kit sequences.
