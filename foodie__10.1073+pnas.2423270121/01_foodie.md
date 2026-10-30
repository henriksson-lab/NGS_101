# FOODIE

[He et al., Genome-wide single-cell and single-molecule footprinting of transcription factors with deaminase](https://doi.org/10.1073/pnas.2423270121), PNAS2024, PMC11670102.

🟢 Published main184–188/282–300 and publishedSI TableS2 actualPDFp6 (printed38) inspected. Mandatory whole-bundle scraper pointers verified against original table. Two originals; automatic preprintv3 was not opened or used. Protein production is upstream of the selected graph. The graph models fresh-cell bulk library preparation; the single-cell adapter component is represented separately, without assigning its identities to the bulk commercial stock.

## Source preparation and chemical meaning

🟢 Start30,000freshK562/GM12878cells, no selected heat-shock/fixation alternative. Source450g5min, coldPBSwash/500g5min4°C. Cold permeabilization10mMTris7.4/10mMNaCl/3mMMg/.1%IGEPAL/.1%Tween/.1%digitonin;550g5min4°C, aspirate supernatant. Tagment with25µl2×TD,1µlVazymeS111,24µlPBS/.1%digitonin,37°C30min800rpm. The2×TDstock is20mMTAPS8.5/10mMMg/20%DMF. Add100µlice-coldRSB,550g5min4°C and retain nuclear pellet. Actual source transposome stock contributes DNA adapters and catalyst; whole loaded oligos are unassembled rather than imported from another Nextera recipe.

🟢 PreparedDddB is the selected footprinting catalyst. Source uses its prepared enzyme/DRB reaction at37°C10min for fresh bulk material, followed by ZymoDNAclean/concentration. Chemical cytosine-to-uracil conversion and protected/incompletely converted template classes remain distinct from thymine subsequently reported by DNA copying/sequencing. Conversion does not establish perfect protection specificity or quantitative TF occupancy. Fixed-cell and DddA comparison branches are not mixed into the selected route.

🟢 Amplify with source1×Q5Ustock plusBst3.0/.08Uµl and1.25µMNextera forward/reverse primers:65°C5min/80°C12min/98°C2min, then98°C15s/60°C30s/72°C1min. At fivecycles, take an actual aliquot for20cycleqPCR, preserving the remaining reaction. The source generally uses10–12TOTALproductioncycles, not10–12additionalcycles. Co-retained primers, catalysts and nucleotides in the remaining reaction are not fresh additions. The qPCR fraction is QC; exact assay mixture/volume remains unassembled. Select above200bp with ZymoSelect-a-Size. Ordinarybulk libraries use2×150PE; the separate cooperativity library uses2×250PE. Whole bulk primers, index choices and read landing remain unassembled.

## Separate single-cell adapter component

🟢 PublishedTableS2 actualp6/printed38 gives P5universal tail `{{= chem.P5_TAIL }}` plus the canonical{{= len(nextera.ME) }}nt mosaic end; P7Tn5-1 tail `{{= chem.P7_TAIL }}`, barcode `{{= chem.TN5_BC1 }}` ({{= len(chem.TN5_BC1) }}nt), then the same mosaic end. The opposing oligo is source5Phos plus `{{= revcomp(nextera.ME) }}`. Whole source orderings are retained and local antiparallel pairing is computed withScene. Tail and barcode are cytosine-depleted; the canonical mosaic end retains cytosines. No protection mechanism for every cytosine, full integration product or resistance to all deamination is inferred from the label “cytosine-depleted”. Barcode meaning is structured as a first-round cell identifier, not a UMI.

🟢 The high-throughput branch uses96tagwells then a physical pool before deamination and a second sorting/index-PCR round; ten nuclei per second well does not equal one nucleus per well. TableS2 also lists many P5/P7PCR primers and custom sequencing primers. Main300 prints10mM P5/P7primer values and a25µl PCR; these literal descriptions are not silently normalized or transferred into the selected bulk graph. Full high-throughput combinatorial assignments, index layouts and source concentration interpretation remain outside this partial implementation. Bulk and single-cell variants are one method entry, not extra protocol counts.

Computed local duplex for the separate single-cell TableS2 component (not the bulk kit):

```
{{scene: chem.single_cell_adapter('p7_1')}}
```

## Limits

🟡 Whole occupied chromatin substrates, adapter integration/end processing, copy/readout geometry, deamination distributions and all unreacted classes remain unassembled. Naked-DNA controls and computational mapping/deduplication do not physically purify the production library. Historical kit documents and predecessor methods are ordinary unread sources; no source-absence audit was accepted. Source comparison/control branches were inspected as boundaries, not inherited recipes.

The directly linked published CommentIn [10.1073/pnas.2425203122](https://doi.org/10.1073/pnas.2425203122) remains an ordinary unread original. Its content, relevance and effect on interpretation have not been assessed; identity alone does not make it a correction, an irrelevant source or a separate protocol. It is a deferred source-review target, not an access-failure request.
