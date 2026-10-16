# PRO-cap — original 2013 cap-selection branch

[The defining Kwak 2013 paper](https://doi.org/10.1126/science.1229386) derives PRO-cap from precision nuclear run-on. All four biotin-NTPs are supplied together; the informative coordinate is the capped nascent-RNA 5′ end.

## Original source and molecular order

🟢 The [author-hosted PDF](https://www.cecerelab.com/_files/ugd/594009_926d2d49ebf94caa8e6cf707c2f4846a.pdf) is the complete original Science 2013 supplement, containing both PRO-seq and PRO-cap. Its printed title/date/DOI identify it independently of the gated PMC download route.

The page follows the original PRO-cap branch:

1. Four-biotin-NTP run-on, RNA extraction/desalting and first streptavidin enrichment, **without** PRO-seq fragmentation.
2. Standard Illumina 3′ RNA-adapter ligation and second streptavidin enrichment.
3. Antarctic phosphatase before TAP; PNK is expressly omitted.
4. Standard 5′ RNA-adapter ligation and third streptavidin enrichment.
5. Standard Illumina RT/PCR, Phusion amplification and PAGE library purification.

The adjacent reverse PRO-seq adapters are assigned to the other branch. They are not silently used as PRO-cap reagents.

## Uncapped material and index boundary

🟡 AP/TAP chemistry permits an uncapped, biotin-labelled 5′-OH outcome. The graph shows this as an inferred logical class with unresolved downstream fate, not as a physically separated fraction. Streptavidin recognizes the run-on label; capture alone does not prove that every uncapped molecule or incompletely adapter-ligated molecule was removed.

The original historical standard oligo identities and sample-index configuration remain unassigned in the molecular model. A shared helper's previous eight-base index was unsupported and has been removed. The local graph declares RT, PCR and purification without importing a generic cDNA-finishing recipe. These are repository modelling corrections, not publication errors.

## Later and dated vendor comparators

🟢 [Mahat 2016](https://doi.org/10.1038/nprot.2016.086), Table 1, prints a concrete later PRO-cap implementation:

- RA3: `{{= pro_cap.RA3_2016 }}-/Inverted dT/`; the table separately states its 5′ phosphate.
- RA5: `{{= pro_cap.RA5_2016 }}`.
- RTP: `{{= pro_cap.RTP_2016 }}`.
- RP1: `{{= pro_cap.RP1_2016 }}`.
- RPI-n: `{{= pro_cap.RPI_2016 }}`.

These declarations are version-scoped comparisons, not inputs to the original 2013 graph. The later protocol adds a 5′-phosphate-dependent exonuclease before phosphatase treatment. That step is not added to the original AP/TAP branch.

🟢 The dated Illumina *TruSeq Small RNA Sample Preparation Guide*, 15004197 Rev C, March 2011, names RA3/RA5/RTP/RP1/RPI and six-base indexes. It describes truncated RNA ligase 2 for its 3′ reaction and a separate 5′ ligation/RT/index-PCR route. A matching reagent name alone does not establish equal modifications, a chosen 2013 kit, or equivalence to the later author-modified RNA adapters.

## Readout and remaining audit scope

Original SOM analysis uses the 5′ ends of **non-reverse-complemented** reads for initiation mapping, unlike PRO-seq. The source names GenomeAnalyzer 2/HiSeq 2000; the exact original loaded run primers and acquired cycles are not declared here. PCR-site architecture is not a loaded-primer assignment.

The audit remains in progress. Remaining original SOM figure/table/script scope and named precursor/source links are recorded in the review. An unassigned historical reagent choice is not an invented unavailable manual, and no publication-wide missing-information claim is accepted.
