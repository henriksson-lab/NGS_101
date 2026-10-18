# DRIPc-seq — strand-specific sequencing of R-loop RNA

## Defining source and inspected scope

- 🟢 Sanz et al., *Molecular Cell* (2016), DOI `10.1016/j.molcel.2016.05.032`, author manuscript PMC4955522. The current migration re-read Experimental Procedures, “DRIP-seq and DRIPc-seq” (JATS 137).
- 🟢 Sanz and Chédin, *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1`, PMC6615061. Reagent/primer passage 610–623 and second-strand setup 1098–1124 were inspected as a **later procedure/comparator**, not as assignments for the 2016 experiment.
- The mandatory bundle scraper returned 20 hits in the 2019 document. Its two PCR primers were verified in that original. Their presence does not establish the 2016 primer identities. The 2016 Supplemental Experimental Procedures, NIHMS800865 supplement-1, is a named original absent from the local bundle; its manual locator and the existing blocked source-review status remain.

## Current 2016 molecular path

1. 🟢 The defining method says DRIP followed Ginno et al. 2012. Extraction, fragmentation and hybrid capture are represented locally with the detailed 2016 conditions unassigned; no family constructor supplies a reagent mixture. S9.6 capture retains the RNA:DNA hybrid and displaced DNA, rather than prematurely removing DNA.
2. 🟢 The eluted material receives **4 U DNase I, 35 min at 37 °C**. The RNA is then ethanol-precipitated. DNA hydrolysis products remain an unresolved molecular fraction; degradation is not physical discard, and their recovery during precipitation is unassigned.
3. 🟢 iScript reverse transcription supermix (Bio-Rad), containing poly(T) and random hexamers, makes first cDNA. RNA/template and cDNA remain distinct molecular populations in the model.
4. 🟢 Second-strand synthesis substitutes dUTP for dTTP. The marked strand is opposite the first cDNA, not an extra contiguous segment of that strand. Specific 2016 enzyme preparations and subsequent RNA fate are unassigned.
5. 🟢 Enrichment at R-loop-positive versus negative loci is checked by qPCR on newly synthesized cDNA. The graph separates an unassigned-volume QC aliquot from material retained for library preparation.
6. 🟢 cDNA is sonicated to approximately 200 bp before library construction. The 2016 method mentions size selection but does not supply an inspected order or complete end-construction implementation; these remain an explicit boundary rather than importing the later 2019 route.
7. 🟢 UDG treatment precedes PCR and provides strand specificity. 🟡 Glycosylase base excision does not by itself prove immediate strand cleavage or physical removal. The affected second strand and excised uracil remain unresolved; the retained first-cDNA precursor is a single-strand role-level state with unknown ends.
8. 🟢 PCR, BioAnalyzer QC and HiSeq 2000 sequencing complete the reported route. The actual incoming PCR precursor is preserved by unavailable geometry; a completed duplex is not substituted as its template. Exact 2016 adapter/primer assignments and acquired read roles, lengths and orientation remain unassigned.

The article states that hybrids must exceed 50 bp to be detected because of size selection. This is a detection limitation, not an assigned library-fragment cutoff or acquired read length.

## Separate 2019 procedure

The later protocol uses **6 U DNase I for 45 min**, names RNase H/E. coli ligase/DNA polymerase I second-strand synthesis and supplies end-repair, dA-tailing, TruSeq Single index, UNG, PCR and selection details. Those are a separate detailed implementation, not a published-source conflict that makes the known 2016 dose uncertain. The 2016 native graph uses its own 4-U/35-min value.

The page lists these exact 2019 primers as comparison inventory only:

- PCR primer 1.0 P5: `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′`.
- PCR primer 2.0 P7: `5′-CAAGCAGAAGACGGCATACGAGAT-3′`.

They are not molecular inputs to the 2016 graph. The later SE100 example is not a 2016 acquisition declaration. Historical 2026-10-17 bounded readings remain in the source-review JSON; their combined-version graph has been superseded by this source-owned 2016 path.

## Remaining questions

The named 2016 supplemental methods and Table S1, exact referenced DRIP preparation, original adapter/primer implementation, second-strand/UDG physical transformations and recovery, and GSE70189 run metadata remain deferred. Later vendor or laboratory revisions do not fill these values. This capped two-original migration is not a complete chemistry review or confirmation of missing published information.
