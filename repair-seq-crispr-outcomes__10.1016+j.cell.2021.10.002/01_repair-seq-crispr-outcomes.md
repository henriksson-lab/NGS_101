# Repair-seq: pooled-screen outcome-library branch

Hussmann et al., *Mapping the Genetic Landscape of DNA Double-strand Break Repair*, DOI10.1016/j.cell.2021.10.002(2021). Current recipe authority is associated NIHMS1750972/PMC9074467.1, not a verified final-body equivalent. The superseded preprint is not merged.

🟢 Main264–270/292 describes a distinct outcome-readout preparation beyond a routine CRISPR screen: NotI digestion and wide genomic target-band recovery retain the linked repair outcome/sgRNA region; HindIII removes one NotI end; annealed oBA676/oBA677 adapters attach a12ntUMI at the remaining4ntNotI-compatible end before PCR. SourceTableS2 is actual linked supplement-9.xlsx, currently gated/TDM. Its whole adapter/primer orderings and chemical modifications are unassembled, not absent. Generic table notes about indexNs are not imposed on uninspected individual oligos.

🟢 gDNA extraction uses50–100millioncells perNucleoSpinBloodXLMaxicolumn. NotI→0.8%agarosewide-band purification→HindIII→0.8xSPRIselect precedes selected1ug DNA in110ulligation(30ulbuffer/10ulligase/200:1adapter:insert;4°Covernight/lid30°C). Sourcecalculated0.65xcleanup retains the PEG-dependent effective-ratio caveat, followed by0.8xcleanup.

🟢 Multiple30–50ng aliquots in50ulPCRs use0.6uM eachP7/indexedScreenindexing primer,3%DMSO and1xKAPAHiFiReadyMix,16or20cycles. The source says nearlyentire stock, not all stock; the model preserves an unresolved remainder, parallel source-derived aliquots and preparedPCRstocks with structural nucleotide components. PostPCR0.8xcleanup, analyticalQC and actual indexed-library pooling remain separate physical boundaries. Pool volumes and target recovery/yield are unassigned.

🟢 Selected pooled-screen NovaSeq6000 run: I1=12ntUMI,I2=8ntsampleindex,R1=45ntCRISPRisgRNA,R2=258ntoutcome. ArrayedMiSeq and companionPerturb-seq preparations are separate branches, not inherited. Computational outcome/UMI filtering is not physical library cleanup.

🟡 The partial adapter-end view retains the known source4ntNotI-compatible interface with unassembled base/strand/UMI geometry; its5prime placement is qualified from NotI cleavage class and inferred by default. It does not invent a full adapter or sealed product. Whole integrated constructs/editing treatments lie upstream of this selected library-prep boundary. Actual TableS2, final-body equivalence and directly linked historical gel guide remain open; no complete source/chemistry/absence audit or public promotion.

Source inspected: external `j.cell.2021.10.002_PMC9074467.xml` selectedMethods264–270/292 and supplemental roles469–477. Mandatory source scraper run; sequence discovery is not source absence evidence.

## Selected source scope

- j.cell.2021.10.002_PMC9074467.xml: AssociatedNIHMS1750972identity; selectedMethods264–270/292, UMIs106 and supplementalTableS2role469–477.

Capped positive/component review only; ordinary unread references and unavailable sources do not establish source absence.

🟡 Selected implementation review: source main270 does not specify whether parallelPCR products are combined before or after0.8×SPRI cleanup. The graph keeps both actual source-derived parallelPCR preparations in an opaque workup boundary; it does not invent a physical pre-cleanup pool. Subsequent mixing of indexed samples is explicit source preparation.
