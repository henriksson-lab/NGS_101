# FragSeq

**Defining source:** Underwood *et al.*, *Nature Methods* (2010),
doi:10.1038/nmeth.1529. The article and method overview were read.

## Evidence used

- 🟢 Purified nuclear RNA is refolded and partially digested with nuclease P1.
- 🟢 Nuclease P1 preferentially cleaves single-stranded nucleic acid, so mapped fragment
  termini report accessible RNA positions.
- 🟢 Digested and undigested-control fragments are converted with a small-RNA sequencing
  workflow and compared computationally.
- 🟢 The supplementary library methods specify ABI SOLiD Small RNA Expression Kit adapters ligated to fragments with 5′ phosphate and 3′ hydroxyl, then one sample-specific 3′ PCR primer (Supplementary Methods, “Library generation and high-throughput sequencing”). Supplementary Data File 4 is the authors' actual read-filter FASTA. Its named `P1` pattern is `NNNNNCCACTACGCCTCCGTTTCCTCTCTATGGGCAGTCGGTGATNNNNN`; `adaptor_bc1_P2` is `NNNNNGCGGAACCGGCATGTCGTCCCCGAATCTCTTACTCCTTGGGCCCCGTCNNNNN`. It also lists 15 other six-base barcode variants. The flanking `N` characters are filter-context wildcards, not claimed synthesized bases.
- 🟡 Those filter patterns establish adapter-associated sequence used by the read-processing pipeline, but do not establish the complete ordered, modified ligation oligos or the SOLiD sequencing primer. The schematic's unresolved oligo roles remain warranted.

## Construct logic

FragSeq is the closest direct RNA analogue of a nuclease-accessibility assay: structure
determines enzyme access, and sequencing locates the resulting molecular ends. Unlike
ATAC-seq, the enzyme cleaves RNA rather than installing adapters in one reaction.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC3247016 JATS67–75, FragSeq methodology and three parallel treatments.

Added the published PNK/ATP endogenous-break control as a third branch beside P1-treated and undigested samples; recorded20–100nt RNA selection and preserved SOLiD platform.

Deferred: Original supplementary detailed methods, explicit gel-selection steps, P1/PNK end-chemistry geometry, vendor oligos and run acquisition. This pass does not establish a complete source review or publication-wide absence.

## Native ownership migration, 2026-10-18

Sources read: Published PMC3247016 JATS Results FragSeq methodology69–73 and Fig1caption317; published supplementaryPDF Fig9 and Note1 mapping/trim-read definition215–231,377–485; published supplementary software filterFile.fasta records1–38 (analysis mask patterns, not oligo orders).

Removed rna_end_readout step/PCR recipes and unused generic workflows. Locally owns P1,undigested andPNK branches,20–100nt gel fractions, P1-generated5′P/3′OH and PNK end-processing semantics, end-compatible RNAadapter ligation, actual single-strandRTproducts and unavailablePCRgeometry on those realprecursors. Physicaloffsizegel fractions distinguished fromunligatedmaterial, whose recovery is unresolved. HistoricalSOLiD endpoint/individualPCRbarcoding remains source-limited with no fabricatedreagent sequences/cycles. FilterFASTA Ns are computationalpadding and were excluded from actual reagent states.

Model boundaries: HistoricalSOLiD adapter/RT/PCR species, precise multi-end ligationorder and nativegeometry, selected barcode/run cycles and linked detailed methods not followed beyond capped originals. Mandatory bundle oligo scraper was run and relevant original passages checked. This capped pass accepts no publication-absence claim.

## Independent review corrections, 2026-10-18

Independent reader inspected Published PMC3247016 Results69–73; published nmeth.1529_supp_NIHMS245940-supplement-supplement_1.pdf.txt Methods35–87. Published SupplementaryMethods places PNK before FlashPAGE20–100nt; specifies post-PCR6%nativeTBEgel selecting25–100nt inserts, equimolarpool, emulsionPCR and ABI SOLiD3 50bp reads. Actual first-strand PCR precursor/unligated fates are sound, but RT template RNA fate needs explicit unresolved output.

Each RT now outputs an explicit molecular RNA-template class with unresolved recovery, while PCR follows the cDNA target component. Moved PNK before FlashPAGE; added branch-specific post-PCR nativegel selection and physical excluded gel fractions. Separate barcoded species are pooled equimolarly, followed by a source-qualified emulsionPCR machinery boundary. Read1 has source50bp biological length; color-space cycle recipe/index acquisition remain unassigned. This correction used that bounded review and the existing source scope; no additional original was read by the implementation reviewer. No publication-absence claim is accepted.
