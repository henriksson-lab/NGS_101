# Combined quantification of intracellular (phospho-)proteins and transcriptomics from fixed single cells

Defining DOI: https://doi.org/10.1038/s41598-018-37977-7

## Published source scope

🟢 Gerlach etal., Scientific Reports9(2019), DOI10.1038/s41598-018-37977-7, published PMC6365588 final article. Read RNA-barcode production/conjugation, fixation/immunostaining and library Methods62–90; computational identification96; workflowResults122. Published SI MOESM1: FigS2/caption and Detailed RAID Library Preparation Protocol for Keratinocytes, explicitly dated2018-05-15, pp8–10 through first-strand cleanup; downstream second-strand/IVT headings were inspected only for scope. Original renderedp8 checked for primer fields/colour highlighting. Mandatory scraper run; compensation primer/ARC marker verified in the article.

This is an RNA-reporter antibody method, not a DNA-tagged antibody recipe. It uses capped polyadenylated synthetic RNA with an explicit15nt reporter UMI and10nt antibody-identity region. The reporter UMI is inherited from a diversified template population; separate IVT copies need not have different reporter-UMI sequences. No unique-per-transcript quantitative efficiency is predicted.

## Reporter production and covalent conjugation

🟢 Source mMessage-mMachineT7 production creates capped reporter RNA, followed by RNA-column purification. FigS2 assigns a28nt polyA tail. Full DNA-template/transcription boundaries, cap atom structure and exact linker geometry are unassembled. The graph begins with this pre-produced reporter, retaining explicit RNA/cap semantics. The12nt identification marker is printed as `ATCAGTCAACAG` in DNA-template/readout spelling; the RNA molecular scene derives `AUCAGUCAACAG` and does not assert ribothymidine residues. Local ligation geometry is explicitly unreviewed/unassembled, not claimed unreported by the original.

🟢 T4RNA-ligase labeling incorporates source pCp-Azide at the RNA3′end overnight16°C. RNA purification removes excess reagent. DBCO-PEG12-TCO is added20fold molar excess overnight at roomtemperature, followed by another RNA purification. The resulting RNA-TCO couples to antibodies previously functionalized with NHS-SS-Tetrazine,1h roomtemperature in borate-buffered salinepH8.4. Upstream antibody functionalization is a named unassembled source input; it is not manufactured by a shared default recipe. pCp/linker chemistry is retained as typed end-label state rather than guessed atom/base-level geometry; the drawing is a partial reporter-body view, not an assertion that the added end label is a nucleotide-free block.

## Fixed-cell handling and first strands

🟢 Keratinocytes are reversibly fixed with2.5mMDSP plus2.5mMSPDP45min, quenched with Tris/salt, blocked/permeabilized with PFBB/tRNA/RNase-inhibitor/0.1%TritonX100, and stained overnight with ARCs. Six gentle10ml washes physically remove unretained material; neither quenching nor staining is treated as purification.

🟢 FACS collects cells into separate primer/mineral-oil-preloaded384well plates. DTT reverse-crosslinking25°C45min then65°C5min releases cellular RNA and the RNA-side ARC label. The selected reducibleSS-containing connection is broken; reporter RNA and5′cap are not deleted. Reduced protein/fixative/linker species remain in the well.

🟢 MaximaHminus produces first strands from endogenous and reporterRNA without covalently fusing those RNA populations. Barcode Compensation Primer `GGGAGACCGACGAAACTGTTAACG` is added with heat/steppedcooling before pooling. Full binding/protection geometry is unassembled. Source first-strand wells are pooled, mineral-oil upper phase removed, lower aqueous phase transferred, then ExonucleaseI37°C20min. DNA-digestion products persist until the source AMPure/PEGNaCl/ethanol physical cleanup. The local graph stops at that first-strand preparation, with possible RNA/reagent carryover unquantified.

## Source primer/run compatibility remains open

The published SIp8 prints an N8 random window followed by highlighted8nt cell barcode `CGTCTAAT`; its full cell-code list is referred to as “supplementalfileF5”. The published main reports63cyclesRead1,14cyclesRead2 and6cyclesIlluminaindex, with transcripts/ARC inRead1 and cellbarcode/UMI inRead2. Without full sequencing-primer coordinates, N8+cell8 cannot simply be reconciled to14cycles. The model does not silently substitute a familiar CEL-seq2 N6 primer, assign an exact RT-UMI width or declare a completed run plan. ReporterUMI15 is a separate explicit identifier and is unaffected. Pipeline read-name swapping is computational, not a reversal or transformation of physical DNA.

The SI IVT recipe also prints “TTP” in its T7-nucleotide table. That label is left uninterpreted; no thymidine-containing RNA or corrected nucleotide recipe is constructed. Full fileF5 identity, second-strand/IVT chemistry, compensated secondRT, two PCR stages, index/read geometry and source-linked SORT-seq/CEL-seq2 originals remain deferred. No absence finding or whole-recipe inheritance is authorized. The source explicitly attributes the frontend to modified SORT-seq/CEL-seq2 and related antibody-conjugation strategies; these are scoped lineage candidates, not exact version equivalence.

## Bounded publication follow-up (2026-10-30)

The published detailed protocol supplies downstream steps that the current first-strand graph does not yet implement: second-strand synthesis (SI383–435), IVT (438–460), compensated/random-octamer RT (501–535),12-cycle preamplification (538–568),6-cycle indexing PCR and final cleanup/Qubit/Bioanalyzer (600–666). The actual RNA→DNA-producing finish is available for modeling; the broad required `unassembled_source` item remains uncovered, not excluded. Full occupied polymers, pCp/linker structure and historical cell-primer/read compatibility remain unassembled. The source's IVT “TTP” label is unresolved, and N8 plus cell-code8 versus14-cycleRead2 is not silently replaced by a familiar CEL-seq2 architecture.

This bounded pass makes only one molecular-inventory correction: source ATP used during pCp-Azide RNA labeling is a molecular structural contribution distinct from the ligase catalyst. Packaging/separate-stock count and dose remain unassigned; this does not invent a new fresh addition or complete the donor/junction chemistry. Two existing published originals only; no acquisition or broader source expansion. Notes status remains a model-scope hold.
