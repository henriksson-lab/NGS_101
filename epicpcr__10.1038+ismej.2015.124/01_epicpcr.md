# Massively parallel sequencing of single cells by epicPCR links functional genes with phylogenetic markers

Defining publication: https://doi.org/10.1038/ismej.2015.124.

## Molecular distinction and scope

🟢 epicPCR links an environmental cell's functional gene to its16S phylogenetic marker by compartmentalized overlap-extension PCR. This local model follows dsrB–16S, not the distinct synthetic barcode–16S control branch. Retention in polyacrylamide hydrogel is source-declared qualitative behavior; no diffusion, occupancy, specificity or recovery calculation is implied.

## Sources inspected

- Published article PMC4737934, Methods: bead formation/lysis, emulsion fusionPCR, nested blockingPCR and sequencing; Results/ Figure1 linkage overview.
- Published supplement `ismej2015124x1.pdf`: SupplementaryMethods, TablesS2/S4/S5. Rendered original printed page18 (PDFpage19) checked for primer colours and C3 modifications.
- Published `ismej2015124x2.pdf`: supplement index/README only. It identifies locations, not molecular evidence.

## Step-by-step and oligos

🟢 Model begins with source-prepared lysed cells in polyacrylamide beads. Emulsion fusionPCR contains F1 and R2 at1µM each and a limiting10nM R1–F2′ bridge. The functional product's added16S homology supports subsequent extension onto the16S template. Actual genomes/landing coordinates are unassigned; a source bridge sequence does not establish complete pairing against an unknown substrate.

🟢 TableS2 dsrB-F1: `GTGTAGCAGTTACCGCA`; dsrB-R1_519R: `GWATTACCGCGGCKGCTGTGCCTSAAYATGTGYGGYG`. Degenerate IUPAC segments are retained as source placeholders, not a selected concrete allele.

🟢 After fusion, emulsions are broken, aqueous material transferred and DNA purified with AMPure. NestedPCR uses forward/reverse0.3µM and3.2µM of each blocking primer. TableS4 oligos are:

- i_dsrB-F3: `ACACGACGCTCTTCCGATCTYRYRVAGVATSGCGATRTCGGA`
- U519R-block10: `TTTTTTTTTTGWATTACCGCGGCKGCTG/3SpC3/`
- U519F-block10: `TTTTTTTTTTCAGCMGCCGCGGTAATWC/3SpC3/`

🟢 Blockers carry3′ three-carbon spacers, not3′ phosphate or extension-capable DNA ends. Their suppression of unfused material is a source-declared design. The model neither computes polymerase competition nor removes DNA merely because it is blocked. Reaction products remain until physical purification. qPCR chooses nestedPCR cycle count; cleaned products then undergo TableS5 Illumina adapterPCR and purification.

🟢 YRYR is sequence-diversity sequence, not a UMI. TableS5 numbered reverse primers carry sample indices. The environmental branch has no modelled cell barcode: physical gene–16S fusion supplies the linkage. MiSeq2×250bp is source-declared.

## Deferrals and component ancestry

Full bead-formation/lysis route, all exact primer pools, qPCR aliquot details, selected sample index, complete fusion products and read-primer landing remain unassembled. No missing-information claim is accepted. Main Methods explicitly attributes the blocking design to Wetmur2005 and Turchaninova2013, with C3 modification and Cradic2004 rationale. This supports a component-level methodological ancestry candidate, not wholesale reagent/recipe inheritance. No import from SiC-seq's different agarose chemistry.

## Completion and verification

Private incomplete source-local graph. Required source inventory explicitly retains unassembled full-strand chemistry; a valid local reaction contract does not establish source-review completion. Mandatory primer scraper was run against fetched sources and relevant hits checked in originals. No two-reader missing-information pass has been completed.
