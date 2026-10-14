# mNET-seq

**Defining paper:** Nojima et al., *Cell* (2015), doi:10.1016/j.cell.2015.03.027. The open author manuscript and its extended methods were read from PMCID PMC4410947.

## Molecular boundary

🟢 Native chromatin is solubilized by MNase, and elongation complexes are immunoprecipitated with antibodies recognizing Pol II or particular CTD phosphorylation states. RNA protected in those complexes is purified and size-selected.

🟢 The informative nucleotide is the RNA 3′ end. It may be the nucleotide in the Pol II active site or a cotranscriptional cleavage intermediate, so the page does not equate every peak with polymerase position.

🟢 RNA 5′ ends are phosphorylated before small-RNA adapter ligation. The original implementation used the Illumina TruSeq Small RNA workflow after end repair.

## Evidence boundary

🟡 The paper establishes adapter roles and read direction but does not make one complete adapter/index oligo set an invariant of mNET-seq. Those regions are therefore role-labelled rather than reconstructed from unrelated modern kit sequences.

## Oligo source audit

🟢 Checked the defining author manuscript (`j.cell.2015.03.027_PMC4410947.xml.txt`, especially lines 83–85 and 669–671) and its extended methods (`j.cell.2015.03.027_supp_mmc2.pdf.txt`, lines 2500–2510). The mNET-seq library uses the TruSeq Small RNA preparation after 5′ phosphorylation and gel selection; the authors sequence from the reverse primer to expose the RNA 3′ end. The extended methods name the kit but do not print its 5′/3′ RNA adapters, RT primer, PCR/index primers or reverse sequencing primer sequence. The scraper's other oligos are siRNA and locus-specific validation primers, not mNET-seq library reagents. Exact kit-era adapter/index sequences remain unresolved without the corresponding historical Illumina guide or a run-specific reagent record.
