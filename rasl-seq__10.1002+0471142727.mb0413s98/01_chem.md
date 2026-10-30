# RASL-seq: selected purified-RNA preparation

Li, Qu and Fu (2012), **RASL-seq for Massively Parallel and Quantitative Analysis of Gene Expression**, DOI10.1002/0471142727.mb0413s98. The actual original read is the associated author manuscript PMC3325489.1/NIHMS368182; no complete equivalence to a final publisher body is asserted.

## Scope and sources

🟢 Associated manuscript planning45–47, materials63–85 and BasicProtocol95–197 inspected. One original, investigation closed. Mandatory oligo scraper ran; whole PCR orderings were checked in the main original. No whole target-probe panel, complete RNA/DNA intermediate or linked-source audit is claimed. The private model selects isolated total RNA, not the separate MELT cell-lysate alternative.

## Physical preparation boundary

🟢 Two DNA probes hybridize adjacently to a target transcript, often across an exon junction, while biotinylated oligo(dT)25 and streptavidin beads capture polyadenylated RNA. After washing, T4DNA ligase joins the DNA probes on the RNA template. Heat releases probes; selected supernatant feeds indexed PCR, then pooled products are column-purified and sequenced. There is no reverse-transcription step in this selected route, and the sequenced probe DNA is not asserted to be a direct RNA copy.

🟢 Planning45 assigns the phosphate to the probe binding the transcript's3′-side target, whose3′ part corresponds to the index-primer arm; planning47 puts the common P5 arm at the5′ end of the other probe. Materials63 calls the5′-phosphorylated species “upstream”; the model does not infer an opposite orientation from that name. Actual gene-specific20+20 target arms and full sealed products remain unassigned. No invented transcript sequence or example probe is labelled a published instance.

## Verbatim historical primer orderings

🟢 Materials85: barcoding primer `5′-CAAGCAGAAGACGGCATACGAG(XXXXXXX)TAGCTGATGCTACGACCACAGG-3′`.
The literal P7-labelled prefix is `CAAGCAGAAGACGGCATACGAG` (not the longer modern canonical P7), the middle window is a seven-base sample code, and the3′ index-primer sequence is `TAGCTGATGCTACGACCACAGG`. Common P5: `AATGATACGGCGACCACCGAGAT`. Canonical constants are sliced only where the complete source ordering verifies them. The actual sample code is unassigned locally and is not a UMI.

🟢 Oligo(dT)25 is biotinylated, but the source's inspected text does not assign an attachment position in this model; a5′ modification is not invented. Both DNA probe and oligo(dT)/bead components are represented as constituents of the single binding-cocktail addition. Known5′ phosphate on the ligation donor is retained, while exact RNA-splinted ends are explicitly unassembled rather than filled from a DNA-templated DASL recipe.

## Step order and physical fates

🟢 Approximately1µg RNA in20µl → add20µl2×binding cocktail →65°C10min/45°C60min → magnetically aspirate, leaving approximately10µl →five50µl washes →one50µl ligation-buffer wash →add20µl1×ligation buffer with5WeissU T4DNA ligase →37°C60min →three washes →resuspend in30µl water →65°C5min →magnetic stand →transfer5µl supernatant to fresh PCR plate. Original plates, including untransferred material, are stored−20°C as backup. Complete eluate purity is not assumed, and the backup is not mislabelled as discarded template.

🟢 Add10µl PCR cocktail with source AmpliTaqGold, both primers, buffer/Mg/dNTP →94°C10min then20–25cycles94°C30s/58°C30s/72°C30s →gel-check a subset of wells →pool PCR products from1–4plates →NucleoSpinExtractII →analytical Qubit allocation →sequence. Stock versus final concentrations in the described additions are not silently recalculated as source facts. Bulk mix components are produced once and contribute their enzyme/primer/buffer roles without duplicate stock additions.

🟢 The physical pool includes other actual PCR wells as an unassembled contributing stock. Random or unrelated ligation products can remain through preparation; later counting of intended ligation products is computational filtering, not perfect DNA removal.

## Endpoint and remaining limits

🟢 Source193–197 describes40nt sequencing of the ligated-probe region, followed by stripping and reprime with the index primer to read the sample code on historical Solexa/HiSeq2000. A modern paired-end/index architecture is not substituted. The source's114nt expected gel band is reported as an experimental expectation, not a checked complete molecular product here.

Complete target-specific probe DNA, RNA-splinted geometry, PCR landing/amplicon strands, actual sample-index choices, historical flowcell/primer implementation and linked precursor/vendor versions remain unassembled or unreviewed. No missing-information claim is accepted. The manuscript's physical preparation is positive evidence for this scoped method; source version and reaction/read limits remain explicit.
