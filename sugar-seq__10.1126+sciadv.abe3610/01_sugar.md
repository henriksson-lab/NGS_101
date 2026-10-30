# SUGAR-seq

Defining original: **SUGAR-seq enables simultaneous detection of glycans, epitopes, and the transcriptome in single cells**, DOI[10.1126/sciadv.abe3610](https://doi.org/10.1126/sciadv.abe3610).

Capped source/model pass,2026-10-29. Private and incomplete;structural graph validity is separate from source completeness. No missing-information claim has passed the two-reader gate.

## Originals inspected

- `sciadv.abe3610_PMC7895430.xml`: Published article Results78–90/Fig1A and Methods160–162: biotinylated/nonbiotinylatedL-Pha, cold washes, anti-biotin/epitope/hashtag staining, 5prime Chromium GEX/TCR/combinedADT–HTO library branches and named indexkit. Separate validation assays not imported into SUGAR acquisition.

## Chemistry and scope

🟢 This preparation converts glycan affinity recognition into an antibody-tag readout. Surface glycans bind L-Pha; biotinylated lectin is then detected by an anti-biotin oligo-antibody. Glycan, lectin, biotin, antibody and DNA tag remain distinct physical components. The source does not establish a covalent glycan–DNA or glycan–RNA product, one DNA tag per glycan, complete binding occupancy or perfect specificity. The local qualitative reporter complex preserves those limits.

🟢 FACS-sorted TILs are incubated with biotinylated/nonbiotinylatedL-Pha at1:5,1µg/ml,20min on ice. Three coldPBS washes precede cellhashing/CITE-seq antibody staining. Source Results identifies anti-biotin,PD-1,TIM3 reporters;stained cells are washed,counted and loaded at~1000cells/µl. Unbound lectin/antibodies are removed only at actual washes. Secondary reporter binding is separate from original lectin–glycan affinity binding.

🟢 The main specifies 5′ SingleCellV(D)J withFeatureBarcoding for GEX,TCR and surface features. This is not a license to reuse a3′ CITE-seq size selection or a modern kit recipe. The source separately prepares GEX,TCR enrichment with the mouseTcellV(D)J kit,and one combined5′HTO/ADT library with the5′FeatureBarcodeLibraryKit, indexed using i7MultiplexKitNSetA. These remain source-named operational continuations with exact capture/RT/primer/fragment/cleanup geometry unassembled. Branch allocations are not ideal modality purification.

The separate genome-wideCRISPR screen,75bp bulkRNA-seq and75bp ATAC validation in other main Methods sections are not SUGAR acquisition settings. The selected SUGAR subsection does not itself assign those run cycles or sequencers to its libraries; this pass leaves acquisition geometry unassembled rather than declaring absent information. Computational doublet filtering does not remove DNA in the protocol graph.

Source limitation: the JATS supplementary-material links exactly `http://advances.sciencemag.org/cgi/content/full/7/8/eabe3610/DC1`. The standard PMC supplementary ZIP contains only F1–F5 image aliases;it is not the chemistry supplement. One declared standard retry of DC1 returnedHTTP403. No guessed supplement filename,challenge bypass or newer revision supplied missing details.

Component ancestry: Methods162 cites CITE-seq staining(ref1) and the named historical5′ manufacturer workflow. The local lectin/secondary-antibody complex is the distinct preparation interface;neither tumor application nor the shared platform establishes complete recipe identity. Historical guides and linked chemistry SI remain pending.

## Deferred originals and questions

- Exact published supplement index http://advances.sciencemag.org/cgi/content/full/7/8/eabe3610/DC1 — actual standard retrieval HTTP403
- Cited CITE-seq staining source and historical 5primeVDJ/FeatureBarcode/MouseTCR/i7MultiplexKitNSetA guides with exact revisions
- Whole antibody/hashtag/capture polymers, identifier widths, primer/end/read/run geometry
- Remaining linked chemistry originals and independent review

## Validation

Local source transcription and physical branch/state checks are recorded in the source-review JSON. Full chemical products,library completion and complete source coverage are not established.
