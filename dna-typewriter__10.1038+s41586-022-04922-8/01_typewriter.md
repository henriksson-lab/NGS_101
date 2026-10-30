# DNA Typewriter

Defining published original: **A time-resolved, multi-symbol molecular recorder via sequential genome editing**, DOI[10.1038/s41586-022-04922-8](https://doi.org/10.1038/s41586-022-04922-8).

Private,incomplete capped source pass (2026-10-29);no missing-information claim accepted. Structural graph validity is separate from source completeness.

## Sources inspected

- `s41586-022-04922-8_PMC9352581.xml`: PublishedVORMethods200–240:cloning/sourceconfigurations,lysisreference,bulk2/3PCRroutes,PacBio route,singlecellCG000205RevD2.3A/B recovery. Results66–100 tape/keyconcept andenhancedPEvariantboundaries.
- `s41586-022-04922-8_supp_41586_2022_4922_MOESM4_ESM.xlsx`: PublishedSupplementaryTables4–5 selectedTAPE1constructs,fullamplicon notation andnamedTAPE_amp_fwd/rev;sequenceinventory not universalallele/product.

## Chemistry and scope

🟢 DNA Tape contains tandem partial CRISPR targets;prime editing installs a variable symbol plus a constant key that enables a subsequent target. The local recorder context does not predict perfect directional editing,efficiency,copy number or the sample allele. PE2/pegRNA and enhancedPEmax/epegRNA/hMLH1dn experiments remain separate source configurations. Recorder creation is upstream of recovery.

🟢 Bulk Illumina: cell lysate1.5µl is allocated to25µl KAPARobust targetPCR25–28cycles;6%TBEassessment,1xAMPure recovery,dual-index/flowcellPCR5–10cycles,AMPure recovery,TapeStation anddenaturation. The source says firstPCR primers include sequencingadaptors at their3′ends;this wording is retained as source wording rather than repaired by inference. Complete product geometry remains unassembled. Fig1single-end vsFig2/3paired-end are explicit source run contexts.

🟢 Separate bulkUMI route starts with a singleforwardprimer10ntUMI extension,1.5xAMPure,then upstreamUMIhandle/locusreversePCR,1xAMPure,indexPCR5–10cycles,AMPure/denaturation. Source calls this twoPCRcycles but details repeated anneal/extendblocks;the model uses single-primer extension,not exponentialPCR. Full chosen UMIprimer remains unassembled;the displayed standardbulkprimer pair is not substituted into this branch.

🟢 Long-read route uses sample-indexed longamplicon primers,KAPARobust25–28cycles with3minextension,0.6xAMPure,TapeStation andSequelCCS. NoIlluminaindexPCR or inventedSMRTbellgeometry is inserted.

🟢 Singlecell source explicitly citesCG000205RevD (NextGEM3′v3.1CRISPRscreening) throughStep2.3. Source16ntCBC/12ntUMI label transcripts. Both2.3A and2.3B eluates containTape (the latter unexpectedly);bothTapecomponents are combined for paired-end targetreadout alongside geneexpressionlibrary. They are not discarded by assuming size equals modality. Whole historic kitrecipe/targetallocation andlibrarygeometry remain unassembled. Component inventories do not define wholebeadprimer order.

🟢 PublishedTable5 `TAPE_amp_fwd` is `GCGTCAGATGTGTATAAGAGACAGGGCTTTATATATCTTGTGGAAAGGACGAAACACC`;`TAPE_amp_rev` is `TTCAGACGTGTGCTCTTCCGATCTttgtctcaagatctagttacgccaagc`. These are standardbulk inventory,not PacBio/UMIbranch substitutions. Table4 referenceconstruct/amplicon strings are not editedsample alleles. InheritedTape,TargetBC integrationidentity,currentcellCBC,UMI andsampleindexes retain distinct meanings.

Ancestry: publishedMethods200 explicitly states the lentiviral TargetBC–Tape–pegRNA construct is based on CROP-seq vector Addgene86708,then describes specific source modifications (GFP–TargetBC–Tape–CaptureSequence1 andInsertBC–pegRNA). This is vectorarchitecture component inheritance,not generic guidecapture recipe inheritance. PublishedMethods240 explicitly cites historical ChromiummanualCG000205RevD;2024DNATypewriterprotocol10.1038/s41596-024-01003-0 is supportinglaterdetail,not a separatemethod or automatic2022recipeequivalent.

## Deferred

- Historical Chromium CG000205 RevD and cellpreparation CG00053 RevC originals
- Linked celllysis/protease predecessor ref2,fullUMI/PacBio/indexPCRoligos andwholeproduct geometry
- 2024 detailed Molecular recording using DNA Typewriter DOI10.1038/s41596-024-01003-0;applicability/revision not assumed
- Remaining originalSI andindependent source review

## Validation

Protocol-local source transcription and branch/order checks;sourceguard and scoped molecular audit recorded in source review. No full base-resolved products or chemistry completion claimed.
