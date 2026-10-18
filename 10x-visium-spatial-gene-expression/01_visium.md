# 10x Visium Spatial Gene Expression

## 1. What it is

Visium captures polyadenylated RNA from a permeabilized tissue section on a glass array whose oligo-dT primers carry known spatial coordinates and UMIs.

## 2. Sources and evidence

- 🟢 10x Genomics, *Visium Spatial Gene Expression Reagent Kits User Guide*, CG000239 Rev F.
- 🟢 The guide defines imaging followed by permeabilization, on-slide reverse transcription, on-slide second-strand synthesis (Step 2.1), denaturation with 0.08 M KOH for 10 min and transfer of the eluate into 1 M Tris pH 7.0 (Step 2.2 a–g; the slide may then be discarded), qPCR cycle determination, cDNA amplification and SPRIselect cleanup (Steps 3.1–3.3), combined fragmentation/end repair/A-tailing, double-sided size selection, adaptor ligation, post-ligation cleanup, sample-index PCR and a final double-sided size selection (Steps 4.1–4.6).
- 🟡 The guide does not state slide attachment chemistry. The model infers that the KOH eluate carries the second strand printed in the Appendix p. 69, because the first strand is primed by the surface oligo.
- 🟢 The sequencing layout is Read 1 28 cycles, i7 10 cycles, i5 10 cycles and Read 2 90 cycles (Introduction p. 17 library figure); Read 1 contains a 16-nt spatial barcode and 12-nt UMI.
- 🟢 The **Oligonucleotide Sequences appendix** (CG000239 Rev F, pp. 70–71) prints the slide capture primer as `5′-CTACACGACGCTCTTCCGATCT-N16-N12-(T)30VN-3′`, with spatial barcode `N16` and UMI `N12`. The model's partial Read 1 is `{{= visium.SLIDE_READ1 }}` ({{= len(visium.SLIDE_READ1) }} nt); round 2026-10-17-breadth-first-06 corrected a repository slice that had added one extra 5′ C to the slide and cDNA forward primers. It prints the TSO `5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′`, preserving three terminal RNA guanosines.
- 🟢 The same appendix prints second-strand primer `AAGCAGTGGTATCAACGCAGAG`, cDNA PCR forward `CTACACGACGCTCTTCCGATCT`, and cDNA PCR reverse `AAGCAGTGGTATCAACGCAGAG` (all 5′→3′). The reverse primer is the same sequence as the second-strand primer.
- 🟢 The ligation adapter is shown as top `5′-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3′` and bottom `3′-TCTAGCCTTCTCG-5′`, leaving a 3′ T overhang on the bottom strand for the dA-tailed insert. The appendix also gives the fixed flanks of both 10-nt sample-index PCR primers; the variable index values depend on the chosen Dual Index TT well.
- 🟢 The [original 10x Dual Index TT Set A table](https://cdn.10xgenomics.com/raw/upload/v1655151898/support/in-line%20documents/Dual_Index_Kit_TT_Set_A.json) publishes the i7 and both i5 workflow orientations for all 96 wells. The guide does not specify which well an individual run used.
- 🟢 [10x Space Ranger spatial-output documentation](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs) states that the spot barcode inclusion list and slide-coordinate mapping are bundled with Space Ranger and written to `tissue_positions.csv`; it lists 4,992 barcodes for the 6.5-mm Visium array. The barcode set is therefore accessible as analysis reference data, although individual spot bases are not recopied into this schematic.
- 🔴 The guide does not identify the slide attachment chemistry or a single fixed spatial-barcode/UMI species; the appendix gives the complete *sequence pattern* and the random/spot-specific windows. [Space Ranger documentation](https://www.10xgenomics.com/support/software/space-ranger/3.1/analysis/spatial-outputs) locates the barcode inclusion list in its distribution and reports the spot-to-barcode map in `tissue_positions.csv`; a physical slide serial and selected TT well remain run-specific.

## 3. Molecular path

The spatial coordinate is installed during reverse transcription, before cDNA leaves the slide. Second-strand synthesis also takes place on the slide; KOH denaturation then moves the cDNA into a tube. Subsequent library construction is bulk, but the inline barcode preserves the originating capture spot.

## 4. Open questions

- Slide attachment chemistry and the strand left on the slide are not printed in CG000239 Rev F.
- qPCR is represented as a consumed 1 µL analytical aliquot with its cycle-selection role; its complete fluorescence/threshold reaction is outside the native graph.
- Native first-cycle PCR geometry remains unassigned for the inferred released single strand and the forked adapter-ligated precursor. Printed oligos and source intermediate patterns are explicit; they do not establish validator coverage of these substrates.

## 5. Native graph migration, 2026-10-18

Read the printed Rev F guide: Appendix oligonucleotide/intermediate panels (text lines 4960–5118), second-strand synthesis/denaturation, qPCR and cDNA amplification (2880–3510), index PCR and final size selection (4253–4415). Ran the bundle oligo scraper and followed its relevant pointers into this original. This pass inspected one original; TT tables and Space Ranger material were not re-read.

The local graph now follows first-strand RT/template switching → second-strand duplex → conditional KOH-released single strand → amplification. The source qPCR aliquot is distinct from the main approximately 35 µL amplification pool. Fragmentation/end preparation produces a population of fragment ends, rather than preserving a full cDNA as an exact cleavage product. The ligation product carries the printed partial Read 1 and forked Read 2 adapter, without P5/P7 or sample indexes; the printed index primers create the indexed final product. Both PCR declarations retain the actual incoming molecular state and explicitly leave native first-cycle geometry unresolved. The final Step 4.6 0.6X/0.8X double-sided selection now records first-cut beads and second-cut supernatant as physical waste, with the final eluate sequenced at 28/10/10/90 cycles. No shared protocol finisher selects these products.

This is a bounded migration check, not an absence audit. Slide attachment and native PCR validation remain open model/source questions; no missing-information conclusion is accepted.
