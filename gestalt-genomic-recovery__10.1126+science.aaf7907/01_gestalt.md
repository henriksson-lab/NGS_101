# GESTALT — original genomic reporter recovery

McKenna et al. (2016), *Whole-organism lineage tracing by combinatorial and cumulative genome editing*, [10.1126/science.aaf7907](https://doi.org/10.1126/science.aaf7907).

## Actual source scope

One original was read: **accepted author manuscript PMC4967023**, Results 65–87, especially the v1 cell-culture recovery/UMI paragraph 79 and the v6/v7 genomic-recovery paragraph 87. Its manifest identifies an author manuscript/TDM source. The standard fetcher and full-bundle oligo scraper were run; no exact primer sequences were obtained from this original. A zero-hit scraper does not establish absent disclosure.

The actual named manuscript **supplemental figures and methods DOCX** and **supplemental tables XLSX** both returned PMC download gates. One bounded standard attempt at the identified institutional published-main PDF returned HTML rather than a valid PDF. Final published main/Methods and supplementary-version reconciliation remain pending. The fetched superseded preprint and its supplements are excluded, not silently used for primers or cleanup. This pass is not a complete source or missing-information audit.

The private graph is explicitly a **conditional accepted-manuscript v1 cell-culture recovery** model. It does not claim a complete final published recipe, infer cleanup or enzyme conditions from later CARLIN, or import scGESTALT's inDrops RNA construction.

## Recorder creation is distinct from recovery

🟢 In this source, the v1 wild-type array contains ten CRISPR/Cas9 targets with 3-base linkers and is 257 bp long, downstream of EGFP. Cell transduction/editing and clonal expansion create the genomic alleles before the represented library-preparation input. The same source separately describes zebrafish v6/v7 arrays with ten independently targeted sites in a DsRed 3′ UTR. Their target identities, exact allele lengths and genomic sequences are not equated with v1.

An edited array is a structured **lineage barcode**: inherited relatives may share it, different lineages may collide, and some earlier mutations may be obscured by subsequent edits. It is not a unique present-day cell barcode or per-copy UMI. The local recorder descriptor records source design context without predicting edits, copying an unreviewed full array, or assigning edited alleles the wild-type length.

## Known genomic recovery

🟢 The source describes adding a **10-base UMI by a single round of polymerase extension**, then amplifying and sequencing the genomic reporter. The recovery UMI is distinct from the inherited recording allele. Full primer ends, adapter arrangement, extension enzyme/program, PCR conditions, intermediate cleanup and final read/index geometry are unassembled because the relevant detailed original is blocked.

These source-known operations appear in the graph, whose initial material is already prepared genomic DNA. Physical reaction species remain unresolved mixtures; no source-unsupported purification, disposal or repeated new UMI assignment is inserted. Reported barcode sequencing is an incomplete readout endpoint, not a claim that crude PCR material is a finished sequencing-ready construct. The source checklist keeps its intervening workup/geometry open.

## Lineage and remaining sources

The existing scGESTALT page covers a later expressed-array/inDrops transcriptome-and-lineage preparation, not this genomic recovery. No constructor or recipe is inherited from it by name. CARLIN's accepted-manuscript Methods later explicitly cites this original for its no-repeat-denaturation UMI-extension strategy; it also states that six source guide designs were retained while four were redesigned. Those are narrow component ancestry statements, not evidence of identical primers or products.

Deferred: actual defining supplemental Methods DOCX and table workbook, final published article, document-version applicability, exact recorder/primer/end/full-library geometry, all directly relevant linked sources and an independent source review. No missing sequence or reagent is confirmed absent. No finished schematic is published.


## Public-page hold adjudication, 2026-10-30

**Required detailed-method/producer hold.** The actual normal graph still
reports `authoritative source item 'detailed_methods' is not accounted for`
and rejects sequencing output `amplified` because it is molecular material,
not a declared library. Merely changing its kind to `K.LIBRARY` would not
resolve the unassembled finishing/source stages. No module or required
coverage changed in this documentation-only adjudication.

The actual inspected original is the associated NIHMS789616manuscript,
[PMC4967023](https://pmc.ncbi.nlm.nih.gov/articles/PMC4967023/), external
`science.aaf7907_PMC4967023.xml(.txt)`, Results75–89/especially79. The
[v1 cell-culture recovery](https://doi.org/10.1126/science.aaf7907) reports a
single polymerase-extension round introducing a10ntUMI before targeted
amplification and NGS. This is qualified accepted-manuscript evidence, not
complete final-version procedural equivalence. The separate v6/v7zebrafish
recovery is not substituted for this v1genomic branch. The downloaded
superseded preprint was not read/adopted.

Actual observed original routes remain unresolved: the
[institutional published-main PDF](https://shendure-web.gs.washington.edu/documents/mckenna_science_2016.pdf)
returned HTML rather than a validated PDF; the named
[Methods/Figures DOCX](https://pmc.ncbi.nlm.nih.gov/articles/instance/4967023/bin/NIHMS789616-supplement-supplemental_figures_and_methods.docx)
and [supplemental Tables XLSX](https://pmc.ncbi.nlm.nih.gov/articles/instance/4967023/bin/NIHMS789616-supplement-supplemental_tables.xlsx)
have actual download-gate manifest rows. Their detailed recipes, enzyme/
nucleotide/buffer supplies, oligos, workups and acquired read/primer geometry
remain unreviewed. No attachment URL is guessed or related-recorder recipe
used in place of these defining originals.

The inherited edited recorder is a **lineage barcode**, separate from the
random10nt recoveryUMI. Wild-type recorder length is not edited-allele length.
The source's molecular/cell-count interpretation of a UMI–allele combination
is not an intrinsic guarantee that every UMI is a unique present-day cellID.
The known UMI component does not determine its complete primer placement.
Recorder creation/editing is upstream input context: genomic recovery is not
blocked on constructing that recorder, but it is blocked on its actual
recovery/library source obligation. Neither later scGESTALT nor CARLIN supplies
missing historical chemistry.

If source completion is revisited, the existing automatic
`tuple(g.materials)+tuple(g.steps)` coverage must become an independently
curated literal inventory. Actual library finishing, controls, physical
stock/aliquot/QC/loading fates need source support before promotion; an
existing PCR amplicon or generic sequencing label is insufficient.

Mandatory scraper `/tmp/unlinked-gestalt11-scraper.log` and actual source
checks preceded adjudication. Six preflight checks passed; chemistry3unresolved/
0contradictory. Code/gates are unchanged, so no tests/builds repeated and no
source acquisition occurred. No absent-disclosure or full independent
source/chemistry claim is accepted.
