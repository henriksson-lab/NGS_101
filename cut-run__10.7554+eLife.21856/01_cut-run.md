# CUT&RUN

## Source and selected implementation

[Skene and Henikoff, eLife (2017)](https://doi.org/10.7554/eLife.21856), published v3 Methods and Appendices 2–3, define several implementations. The schematic selects the soluble yeast branch with an optional spike-in and pellet QC. Mammalian, whole-reaction total-DNA and high-salt variants are separate source branches.

🟢 The yeast branch binds antibodies and pA–MNase to bead-bound nuclei, activates cleavage, stops digestion with chelators, releases fragments and separates the supernatant. Extracted supernatant DNA enters library preparation. Appendix 2 explicitly offers pellet extraction for gel analysis; the pellet is not obligatorily discarded. Failure of high-molecular-weight DNA to amplify is not physical removal.

🟢 Methods line 376 names KAPA library preparation, no size-selection, eight or more PCR cycles and 25-cycle paired-end HiSeq 2500 acquisition. The total-DNA variant uses negative half-volume AMPure selection: retain supernatant, remove bead-bound large DNA. Its 14-cycle PCR preferentially amplifies shorter fragments.

🟢 Re-read 2026-10-17 in the eLife v3 PDF twin, Methods "CUT&RUN for yeast" (lines 781–794) and Appendix 2 (lines 1427–1470): digestion is started with 3 µL 100 mM CaCl2 (2 mM final) at 0 °C in ice-water and stopped with an equal volume of 2XSTOP — 200 mM NaCl, 20 mM EDTA, 4 mM EGTA, 50 µg/ml RNase A, 40 µg/ml glycogen — optionally containing 5–50 pg/ml mononucleosome-sized spike-in DNA from crosslinked MNase-treated *Drosophila* chromatin. Excised fragments are released at 37 °C for 20 min, the beads are centrifuged 5 min at 13,000 rpm at 4 °C and the supernatant is taken off on a magnet stand. Supernatant DNA is deproteinized with 0.1% SDS and Proteinase K at 70 °C for 10 min, extracted with phenol–chloroform–isoamyl alcohol in a phase-lock tube, re-extracted with one volume of chloroform, precipitated with glycogen and 2–2.5 volumes of ethanol and dissolved in 25 µL 0.1× TE8. The optional pellet route takes the beads up in 1XSTOP without RNase or glycogen and ends with a 10 min spin to pellet insoluble bead material. These conditions are now carried by the declaration's reagent labels rather than by prose alone.

🟢 The library paragraph (lines 839–846) was re-verified: libraries follow Kasinathan 2014 / Henikoff 2011 but without size-selection, use the KAPA DNA polymerase library preparation kit protocol with eight or more cycles, and are sequenced 25 cycles paired-end on HiSeq 2500. For total-DNA samples, half a volume of AMPure XP depletes large fragments (supernatant retained, beads discarded) and PCR is raised to 14 cycles **with adapter concentrations increased accordingly** — that adapter-concentration clause is now shown in the amplification-preference notice on the page.

## Library evidence and model boundary

🟢 Kasinathan 2014 Supplementary Fig. 1c shows end repair, A-addition, “Illumina Tru-seq adapters”, PCR and paired-end sequencing. Henikoff 2011 Methods names Illumina PE-930-1001 and modified cleanup. These are cited predecessor implementations, not a selected 2017 adapter/index/ordered-primer assignment.

🟡 The declaration therefore leaves ordered library ends, adapter-ligation geometry, PCR binding geometry and index acquisition unassembled. It asserts the source's 25/25 acquired reads without a fabricated index width or member. Earlier repository claims of an established generic TruSeq shell and eight-base single index were modelling errors, not publication errors. The linked KAPA URL currently returns an unrelated Solix FAST qPCR shop page; that rejected response is quarantined, not treated as a historical kit manual.

## Deposited protocol differences

🟢 All seven original GEO GSE84474 RTF bodies were independently read twice. The file named `Cut-and-Run_total_DNA_protocol` instructs collecting the supernatant before extraction and negative AMPure selection, whereas article Appendix 2's total-DNA branch extracts the total reaction. This unresolved label/body difference is displayed on the schematic; no typographical error is assumed. The displayed workflow uses the article's soluble yeast branch.

The extract RTF adds its Drosophila spike after separation; Appendix 2 offers a spike in STOP. These are distinct implementations, not one combined spike-in schedule. PEG and spike-in RTFs include negative AMPure selection. The human size-selection RTF explicitly extracts the whole reaction before retaining unbound DNA below approximately 700 bp. Native ChIP is a separate comparator method.

## Source checklist and remaining work

The located two-reader checklist is maintained in `catalogue/source_reviews/cut-run__10.7554+eLife.21856.json`. Current originals include eLife v3, its appendices, Kasinathan's original supplement, Henikoff's original main article, GEO metadata and all seven deposited RTFs. Henikoff's identified original SI remains route-blocked; the current KAPA target fails document identity. GEO's cross-link extraction filename, actual analysis-repository children and comparative dataset scope remain explicit linked obligations. No missing-information claim is accepted from these incomplete source checks.
