# Targeted Perturb-seq enables genome-scale genetic screens in single cells

Defining publication: https://doi.org/10.1038/s41592-020-0837-5.

## Molecular distinction and scope

🟢 TAP-seq enriches selected genes and CROP-seq guide-encoding transcripts from already cell-barcoded cDNA. It does not modify capture beads like DART-seq. This model follows the10x3′ reagentkitv2/CROP-seq route; the separately described Drop-seq and mouse bone-marrow/no-CROP variants remain distinct and deferred.

## Sources inspected

- Published study's accepted-manuscript PMC7610614, Methods: Cell barcoding and reverse transcription; PCR1/PCR2; PCR3; Illumina sequencing (text189–207).
- Published SupplementaryInformation MOESM1, TAP-seq step-by-step PCR1–PCR3 section. Its compatibility discussion mentionsv2/v3/NextGEM; that is not permission to substitute a newer kit into the main study'sv2 route.
- Published TableS6 (publisherMOESM8), universal primers and chemical modification comments. Manifest confirms exact MD5 identity with the PMC `Table_S6.xlsx` cache. TableS3 is enhancer-screen results, not the oligo table; it is not used for chemistry.

## Workflow and identities

🟢 Upstream10x3′v2 omits the template-switch oligo from its master mix. The graph begins with35µl cleaned GEM-RT cDNA rather than importing a contemporary Chromium recipe. A separate large-screen modified loading/SuperscriptIV route is outside this standard branch.

🟢 PCR1 uses gene-specific outer primers, CROPouter and PartialRead1 with KAPAHiFi,11cycles (12for a separate large-screen condition). Products receive1.5×AMPureXP cleanup, two80%ethanol washes and30µl elution.10ng enters PCR2 with nested inner gene/CROP primers and the same PartialRead1 reverse primer,8cycles; repeat cleanup.10ng then enters PCR3 with Targeted10X and an N7XX reverse index primer,8cycles, followed by the same cleanup.

🟢 CROPouter/CROPinner are added eightfold over each individual gene primer. They amplify the vector-derived PolII transcript in the same targeted library, not a fictitious separate guide-library branch. Existing cell barcode/UMI identity is preserved; no new UMI is introduced during targetedPCR. Untransferred aliquot remainders are explicitly unresolved, not discarded without separation.

## Source oligos and modifications

- PartialRead1: `CTACACGACGCTCTTCCGATCT`
- CROPouter: `AATGGACTATCATATGCTTACCGT`
- CROPinner: `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTGTGGAAAGGACGAAACACC`
- Targeted10X: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTC*C*G`
- N7XX: `CAAGCAGAAGACGGCATACGAGAT[X]8GTCTCGTGGGCTCGG`

🟢 Asterisks specify phosphorothioate linkages between the final three Targeted10X bases. They are retained in the reagent label/source test and notes; full bond-level chemistry is not simulated. N7's8base window is structurally tagged as sample index; no arbitrary index choice is supplied. Canonical P5/P7 and Nextera components are imported when assembled.

🟢 Main Methods assigns Read1 cell barcode+UMI26cycles, Read2 cDNA approximately58cycles, and i7 sample-index8cycles using standard Illumina read primers.58 is a representative approximate schedule, not a requirement for every implementation.

## Deferrals and lineage

Full historical10xv2 input architecture, target panels inTableS1, exact outer/inner landing and final read-primer coordinates remain unassembled. Relevant documents are available; no missing-information claim is made. Documented component ancestry candidates are10x3′v2 capture/barcoding (with explicit omission), Drop-seq in a separate implementation, and CROP-seq guide-transcript capture. None licenses wholesale shared-recipe inheritance.

## Completion and verification

Private incomplete source-local graph. Required source inventory explicitly retains unassembled full-strand chemistry; a valid local reaction contract does not establish source-review completion. Mandatory primer scraper was run against fetched sources and relevant hits checked in originals. No two-reader missing-information pass has been completed.
