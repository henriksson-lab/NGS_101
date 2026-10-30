# Circulome-Seq

## Sources read and boundaries

🟢 `g3.117.300141_PMC5633380.xml.txt`: PublisheddefiningoriginalMethods90–106/Fig1 and initialResults; selectednativecircleenrichment/tagmentation/PCR, noRCA.

🟢 `g3.117.300141_supp_3295FigureS1.pptx`: OriginalPPTXslideXML caption: exoV/gradient/combinedenrichmentcomparison and technical/biologicalreplicate coverage. Plots not numerically reanalyzed.

🟢 `g3.117.300141_supp_3295FigureS2.pptx`: OriginalPPTXslideXML caption: synthetic3kb/412bpcirclecontrol, explicit421bp filledinsert=412+9duplication. No synthetic-template recipe imported.

Cappedpass: at most10minutes/3originals; definingmain/chemistrysupplements before linkedoriginals. Ordinary deferredhistoricalgeometry is not an accessbarrier.

## Native circle enrichment

🟢 Circulome-Seq works on native extracted DNA and does not use RCA. The original offers ExoV, CsCl-ethidiumbromide topology separation, or both as branches. ExoV digestion enriches circles relative to linear DNA; surviving linear species and digested residues are not assumed absent. The selectedC.elegans example uses400U/µg over3days with70°C30min heat inactivation after rounds. The human5day/500–1000Uµg cases are source conditions, not separate methods. Heat inactivation does not physically remove enzymes or nucleotides. Exact transfer/cleanup before tagmentation is unreviewed and remains an opaque boundary.

🟢 The gradient route uses CsCl1.55g/ml,500237g2.5h and a separate plasmid-reference tube. Collection targets the approximate covalently closed circle band. The other fraction contains linear DNA **and nicked/relaxed circles**; it is not reclassified as pure linear waste or silently thrown away by the graph. CsCl/TE-saturated1-butanol removes ethidium, then TE dialysis2days4°C follows. Absolute topology specificity/purity and complete molecular partition are not guaranteed.

## Tagmentation, copying and acquisition

🟢 Source low-input NexteraXT uses0.5µltagmentase37°C30min; the singly-tagmented-circle enrichment branch attenuates to3min without shaking. The genomic-control1ng/1.5µl30min recipe is separate. Source PCR-cycle choice is based on gel monitoring, with10–12cycles reported sufficient. Exact kit primer/adaptor bodies, fills, cleanup and raw run schedule remain unreviewed. The model does not import canonical current Nextera arms merely because the brand matches.

🟢 MainFig1B and originalFigureS2 establish a conditional singly inserted circle product: Tn5's9nt offset is filled during copying, creating a9bp duplicate across the linear insert boundaries. The local insert-only relation rotates an explicitly supplied circle and appends the first9bases, yielding circlelength+9. This is a selected post-fill-in class, not proof that all native circles receive exactly one event or predict exact insertion positions, full adapters or ligation ends. Opposite read directions require orientation handling; identical top-strand boundary strings do not imply identical5′read strings. Native circles, gapped tagged intermediates and linear copied libraries remain distinct.

The actual source uses paired-end read mapping and the duplication signature. Computational circle classification/deduplication/enrichment does not physically remove other molecules. Separate untreated genomic controls and the synthetic3kb/412bp validation substrates are not silently merged into the experimental eccDNA construct. Exact source platform/read cycles/primers remain unassembled.

## Inspected supplements and remaining scope

The cappedpass inspected three originals: main and FigureS1/S2 originalPPTX captions, with no recursive source acquisition. FigureS1 compares methods/replicates; FigureS2 explicitly says421bp=412bp+9. OtherPPTX figures remain unreviewed. FileS1/S2ZIP member inventories identify positional/repeat summary tables; their data values were not chemically reviewed and the inventory is not an absence audit. Mandatory scraperzerohits reflects availabletext conversions, not proof that all original tables lack sequences.

The defining article cites genericTn5/Nextera and ExoV/topology predecessors; no component-level recipe inheritance is established solely by recorder similarity. Historical NexteraXT guide and exact ExoV-transfer/cleanup geometry remain direct deferred targets. Circle-Seq/CIDER source chemistries are not imported. No absence claim or public promotion.

## Actual retrieval history

- https://www.biorxiv.org/content/10.1101/128686v1.supplementary-material — supplementary page; HTTP 0
- https://www.biorxiv.org/content/10.1101/128686v1.full.pdf — biorxiv full text, v1; host cool-down after HTTP 429 (3.2 h left; not requested)

Availablepublishedsupplements were read where stated; failedPMC aliases and supersededpreprints are retrievalhistory, not requests to reacquire available/superseded originals. No unsupported source-absence conclusion.

## Local validation

Privategraph-onlyimplementation; no page builder or completionmarker. Focusedtests cover selected chemical identity/copying/fates. Fullhistoricalmolecular/end/readgeometry remains unassembled. Exact sequenceguard results and independent critical review will be recorded in the source-review JSON.
