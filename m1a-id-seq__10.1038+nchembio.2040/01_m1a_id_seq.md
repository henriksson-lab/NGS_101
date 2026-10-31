# m1A-ID-seq

Li et al.(2016), *Transcriptome-wide mapping reveals reversible and dynamic N1-methyladenosine methylome*. DOI https://doi.org/10.1038/nchembio.2040 .

## Sources and boundary

🟢 Publisher main title/abstract and Fig1caption establish antibody enrichment, native m1A RT stalling and demethylation comparison. The full main Methods were not obtained: the standard source tool records the defining DOI as manual/no open fulltext. This is not an inspected main protocol.

🟢 Published SI MOESM157, title/authors and Table1 p2 inspected; Fig2–5 pp4–7 source109–169 inspected. Table1 contains qRT-PCR validation primers, not established sequencing-library adapters. Mandatory scraper run; several hits concatenate neighboring table rows, so they are not used as a library sequence. The workbook opening is a peak-coordinate/results table; the workbook is otherwise unreviewed.

The rendered original SI p1 (not merely its extraction) prints footer DOI `10.1038/nbt.2040`, while matching title/authors and the publisher's supplement URL associate the file with `10.1038/nchembio.2040`. This bibliographic footer anomaly does not establish a chemical disagreement or justify switching sources. No public finished schematic is asserted.

## Positive molecular scope

🟢 SI Fig3 explicitly optimizes demethylation of immunoprecipitated m1A-containing RNA: pH6.5 and protein/RNA2:1; BSA is omitted because commercial BSA contains RNA contamination. The selected SI passage does not by itself identify the exact demethylase stock/cofactor/temperature schedule. The abstract identifies ALKBH3 biology; that does not establish that the library demethylation reagent is ALKBH3. The model leaves its identity unassembled.

🟢 SI Fig5 compares AMV and SuperScriptIII, selects AMV for signal/background and reports0.2mMdNTP/10mMMg optimized conditions. The local graph starts at the source IP-RNA boundary, then represents the source demethylation comparison and AMV copying. Intervening recovery and downstream library assembly/run geometry remain unresolved because main Methods are blocked. It does not import the contemporaneous m1A-seq NEBNext recipe.

🟡 `SelectedRNA.demethylated(positions)` is a selected-site chemical relation: supplied m1A residues become A while RNA coordinates remain. It neither predicts enzyme selectivity/conversion fraction nor removes other modifications. `selected_stop(index)` copies the supplied downstream suffix into complementary firstDNA orientation; it does not predict universal stalling, a full molecular RT/PCR product, or physical RNA cleavage.

## Distinctness, lineage and open work

This source-supported demethylation/AMV comparison differs from Dominissini2016's published SI alkaline m1A→m6A comparison. These two independent originals therefore have a chemical boundary beyond their names. Native m1A can cause multiple RT outcomes; no universal mismatch is assigned. Later m1A-MAP is not imported as a library recipe. Complete main Methods, exact protein/cofactor recipe, adapters/end processing/primer identities and acquisition remain blocked/unreviewed, not confirmed absent. No two-reader source/linked-material absence audit.

## Bounded publication hold, 30 October 2026

SIestablishes IP-RNA demethylationcomparison, pH6.5/protein:RNA2:1 andnoaddedBSA, AMV0.2mMdNTP/10mMMg optimizedbackground. Exact demethylase/cofactor/fullreactionorder andactualNGSlibraryproducer are notestablished byselectedSI. Table1is qRTPCRvalidation, not NGSadapterinventory. Keep blocked_mainrequired/uncovered; do notinventRT/PCR/libraryfinish fromcontemporaneousm1Aseq/MAP. Existing printedSIfooter nbt.2040 versus title/authors/publisherattachmentnchembio.2040 qualification retained; not a chemistry conflict.

Required item **blocked_main** remains unchanged and uncovered: Real required definingMethods unavailable/unreviewed; selectedSIcomparisons do not supplyactualNGSlibraryworkflow. No source-obligationexclusion or absenceclaim. The source investigation stopped at the capped original scope; no source/body expansion, acquisition, gate exclusion or public promotion occurred. Current model/renderer code was left unchanged.
