# CUT&RUN

## Source and selected implementation

[Skene and Henikoff, eLife (2017)](https://doi.org/10.7554/eLife.21856), published v3 Methods and Appendices 2–3, define several implementations. The schematic selects the soluble yeast branch with an optional spike-in and pellet QC. Mammalian, whole-reaction total-DNA and high-salt variants are separate source branches.

🟢 The yeast branch binds antibodies and pA–MNase to bead-bound nuclei, activates cleavage, stops digestion with chelators, releases fragments and separates the supernatant. Extracted supernatant DNA enters library preparation. Appendix 2 explicitly offers pellet extraction for gel analysis; the pellet is not obligatorily discarded. Failure of high-molecular-weight DNA to amplify is not physical removal.

🟢 Methods line 376 names KAPA library preparation, no size-selection, eight or more PCR cycles and 25-cycle paired-end HiSeq 2500 acquisition. The total-DNA variant uses negative half-volume AMPure selection: retain supernatant, remove bead-bound large DNA. Its 14-cycle PCR preferentially amplifies shorter fragments.

## Library evidence and model boundary

🟢 Kasinathan 2014 Supplementary Fig. 1c shows end repair, A-addition, “Illumina Tru-seq adapters”, PCR and paired-end sequencing. Henikoff 2011 Methods names Illumina PE-930-1001 and modified cleanup. These are cited predecessor implementations, not a selected 2017 adapter/index/ordered-primer assignment.

🟡 The declaration therefore leaves ordered library ends, adapter-ligation geometry, PCR binding geometry and index acquisition unassembled. It asserts the source's 25/25 acquired reads without a fabricated index width or member. Earlier repository claims of an established generic TruSeq shell and eight-base single index were modelling errors, not publication errors. The linked KAPA URL currently returns an unrelated Solix FAST qPCR shop page; that rejected response is quarantined, not treated as a historical kit manual.

## Deposited protocol differences

🟢 All seven original GEO GSE84474 RTF bodies were independently read twice. The file named `Cut-and-Run_total_DNA_protocol` instructs collecting the supernatant before extraction and negative AMPure selection, whereas article Appendix 2's total-DNA branch extracts the total reaction. This unresolved label/body difference is displayed on the schematic; no typographical error is assumed. The displayed workflow uses the article's soluble yeast branch.

The extract RTF adds its Drosophila spike after separation; Appendix 2 offers a spike in STOP. These are distinct implementations, not one combined spike-in schedule. PEG and spike-in RTFs include negative AMPure selection. The human size-selection RTF explicitly extracts the whole reaction before retaining unbound DNA below approximately 700 bp. Native ChIP is a separate comparator method.

## Source checklist and remaining work

The located two-reader checklist is maintained in `catalogue/source_reviews/cut-run__10.7554+eLife.21856.json`. Current originals include eLife v3, its appendices, Kasinathan's original supplement, Henikoff's original main article, GEO metadata and all seven deposited RTFs. Henikoff's identified original SI remains route-blocked; the current KAPA target fails document identity. GEO's cross-link extraction filename, actual analysis-repository children and comparative dataset scope remain explicit linked obligations. No missing-information claim is accepted from these incomplete source checks.
