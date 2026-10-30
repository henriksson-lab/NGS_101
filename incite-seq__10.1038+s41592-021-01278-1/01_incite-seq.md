# inCITE-seq: intranuclear protein/RNA preparation

Chung et al., *Joint single-cell measurements of nuclear proteins and RNA in vivo*, DOI10.1038/s41592-021-01278-1(2021). Current selected recipe is associatedNIHMS1745895/PMC8532076.1; full final-body equivalence is unverified. PublishedMOESM1 is a reporting summary, not an oligo table, and is recorded as such after inspection.

🟢 Main173–209 establishes intranuclear fixed/permeabilized antibody-DNA staining with optimized dextran-sulfate blocking and simultaneous nucleus hashing before historical10X3primev3 RNA/antibody-tag capture. This is a preparation adaptation beyond changing a protein panel: surface-only CITE handling is not substituted. Antibodies are sourceTotalSeq-Aformat; whole molecular conjugate/primer/index/barcode orderings remain unassembled, not source-absent.

🟢 EZnucleus extraction uses2mllysis/douncing24A+24B, additional3mllysis,5minice/pellet500g5min4°C, repeat5mllysis/ice/pellet, then1mlPBS/Mg/RNaseinhibitor and35umfilter. Source3ml1.33%FA-NT is added to1mlnuclei,10min4°C;3ul1Mglycine and20umfilter precede850gpelleting.500uloptimizedblockingbuffer15min4°C is followed by200ulstainingbuffer1h4°C. Source1%BSA/FcX carrier/blocking components belong to those actual buffer additions; they are not silently omitted or added twice. Two PBST/RNaseinhibitor washes precede sequencing; secondary fluorophore/DAPI validation is a separate branch.

🟢 Simultaneous reverse-crosslinking/RNA RT and DNA-tag copying occur in sourceGEM workflow53°C45min→85°C5min. AMP mix includes HTOadditive0.1uM and ADTadditive0.2uM.0.6xSPRI split physically retains WTA on beads and antibody-tag products in supernatant. WTA is eluted40ul;10ul selected for historicalGEX construction, remaining30ul retained. Supernatant1.4xcleanup yields24ul; separate6ul aliquots enterHTO andADT indexPCRs, with12ul remainder. Each uses sourceNEBNext2xMastermix,21cycles and2.0xcleanup/elution20ul. No whole contemporary10X or CITE constructor is imported.

🟢 GEX/ADT/HTOlibraries receive Qubit/TapeStationQC, then selected-library pools are sequenced onNextSeq500 with sourceR1:28/I1:8/R2:55. The source75cyclekit descriptor is not forced into a new cycle mask. Full oligos, sourcebarcodevalues and primer landings remain unassembled. No source absence/full chemistry validation or public schematic readiness is claimed.

Sources inspected: associated `s41592-021-01278-1_PMC8532076.xml` selectedMethods173–209 and declaredTotalSeq-Aformat197; published `s41592-021-01278-1_41592_2021_1278_MOESM1_ESM.pdf` reporting-summary identity/antibody/study scope only. Historical10X/vendor/nucleus-hashing predecessor references remain ordinary unreviewed.

## Selected source scope

- s41592-021-01278-1_PMC8532076.xml: AssociatedNIHMS1745895identity; selectedMethods173–209, TotalSeq-Aformat197, directnuclearstaining/capture/physicallibraryfractions.
- s41592-021-01278-1_41592_2021_1278_MOESM1_ESM.pdf: Published reporting-summary identity/study/antibody/flowcytometryscope. Verified reportingsummary, NOTan oligotable or central chemistryrecipe.
- s41592-021-01278-1.html: Published finalNatureidentity/abstract/EDcaptions read; actualfullMethods notexposed: subscriptionpreview.

Capped positive/component review only; ordinary unread references and unavailable sources do not establish source absence.
