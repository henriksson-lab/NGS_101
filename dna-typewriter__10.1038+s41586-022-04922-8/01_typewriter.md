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


## Bounded publication adjudication — 2026-10-30, batch 13

🟢 Selected bulk route:1.5ul cell lysate to25ul KAPA Robust targetPCR25–28cycles,6percentTBE assessment/1xAMPure,5–10cycle dual-index PCR, cleanup/TapeStation/denaturation and NextSeq/MiSeq. Source distinguishes bulk single-primer10ntUMI extension, long-amplicon PacBio and single-cell recovered tape from BOTH2.3A/B eluates. Cell-level CBC, molecule UMI, integrated TargetBC/Tape and sample indexes have different meanings. Recorder editing is upstream of sequencing recovery, not repeated by a generic PCR.

🟡 Original required unassembled remains uncovered: Unreviewed linked originals and full molecular geometry. Source LIMIT explicitly includes actual edited alleles, full library duplexes, historical CG000205RevD details, UMI primer identity and PacBio kit geometry. The bulk7step model does not fulfill that original breadth; bulk final-index typing cannot replace or narrow it. Genuinely cited cell-lysis/protease predecessor and historical single-cell guide versions remain uninspected in this pass. Four generic construct-changing steps plus non-library terminal also fail normal publication gates. Runtime graph-enumerating source coverage is not an independent source inventory.

Inspected sources: Published PMC9352581 main JATS/PDF: Version/title header; Methods200–218 recorder production and cited cell lysis; Methods218–240 full bulk2PCR, separate UMI/PacBio routes and actual single-cell two10kcell reactions throughCG000205RevD2.3A/B; source experiment-specific single/paired read contexts; Published MOESM4 XLSX: Actual original sharedStrings XML TAPE_amp_fwd/rev whole orders plus SupplementaryTables4–5 reference construct/standardPCR/lineage/long-amplicon role rows.

Remaining actual source/model obligations: Document original full linked-source/geometry hold. A narrower bulk scope may be useful later but cannot silently rename or exclude the existing broad required item. Source-known KAPA Robust catalyst and its actual nucleotide/buffer supplies need explicit inventory and typed PCR with honest unknown geometry; currently targetPCR lacks any catalyst, while indexPCR hides primer/catalyst inventory in one label. No invented master-stock packaging. Keep real lysate allocation/remainder; source TBE/TapeStationQC aliquots and actual library contribution/loading remainder must be separate physical fates, not just selection-reagent labels. Source says firstPCR primers include sequencing adaptors at3primeends; preserve this source wording and actual Table5 whole orders without silently rearranging the primer or calling the paper erroneous solely from a generic expected architecture. Respect distinct bulkUMI/PacBio/single-cell routes. Single-cell source actually specifies2reactions at10000targeted cells each, both2.3A/B tape fractions recovered, and0.04percentBSA/PBS preparation. Broad source items cannot be satisfied by a bulk terminal. A later2024 detailed protocol is genuinely referenced in existing review history but its applicability/revision is not automatically2022 equivalence. No new retrieval or inference from platform names.

No module changes are made in this adjudication. Fixing the terminal type or changing generic operations to opaque transformations would not establish the original required molecular scope. A future model pass must use an independently reviewed literal source inventory, rather than treating existing runtime enumeration of graph nodes as source evidence. The broad required item remains uncovered; unknown components are not renamed, excluded or covered by a library label. Genuine linked originals and full geometry remain unreviewed/unassembled, not confirmed absent. No sources were acquired and no architecture expansion was undertaken. Existing private builders and normal publication holds remain intact.
