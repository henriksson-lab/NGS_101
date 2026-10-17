# HiChIP — research notes

## 1. What it is

HiChIP combines in situ proximity ligation, protein-directed immunoprecipitation and
on-bead Tn5 library construction to enrich chromatin contacts involving a chosen target.

## 2. Sources read

- 🟢 Mumbach et al. 2016, DOI [10.1038/nmeth.3999](https://doi.org/10.1038/nmeth.3999),
  complete author manuscript and Online Methods (PMC5501173).
- 🟢 The defining paper's supplementary HiChIP protocol gives Nextera Ad1.1 (universal) and Ad2.x (barcoded) as the first PCR primers; the Online Methods calls them Ad1_noMX and Ad2.X. Supplementary Table 4 publishes the full Ad1_noMX `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG` and all 24 Ad2 indexed members; for example Ad2.1_TAAGGCGA is `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT`. The Online Methods also prints the optional post-pooling short primers verbatim: Nextera_i7Short `CAAGCAGAAGACGGCATACGAGAT` and Nextera_i5Short `AATGATACGGCGACCACCGAGATCTACAC`. These short primers perform additional amplification after indexed primer incorporation; they do not introduce the indexes.

## 3. Molecular order

🟢 Crosslinked chromatin is restriction-digested in nuclei, ends are filled with biotin,
and spatially adjacent ends are ligated. The chromatin is then sheared and
immunoprecipitated. 🟢 ChIP DNA is eluted, reverse-crosslinked, column-purified and captured on
streptavidin beads before Tn5 tagmentation. These beads carry biotin-bearing
purified DNA rather than the original immunoprecipitated chromatin.

## 4. Evidence boundary

🟢 The defining implementation uses MboI; the page derives its filled GATCGATC contact
scar and top-strand biotin coordinate from the shared restriction model. The final adapter
architecture is the canonical Nextera molecule implied by Tn5 library construction.
🟢 Supplementary Table 4 lists the complete Ad2.1–Ad2.24 primer family. The selected
Ad2.X member varies by sample; the published short-primer pair is fixed across samples.

## 5. Bounded source pass — 2026-10-17

Read published article archived manuscript PMC5501173 JATS Methods77–111:
fixation, contact generation, IP, elution/reverse crosslinking, biotin capture,
tagmentation, PCR and size selection. Read published Supplementary Table4
`MOESM258_ESM.xlsx.txt` heading and representative Ad1_noMX/Ad2.1 rows1–3.
Manifest identities distinguish these from the excluded superseded 073619
preprint. Ran the bundle oligo scraper; used hits were verified in published
JATS107–109/workbook2–3 rather than preprint copies.

🟢 The source order is IP → elution/proteinase K/reverse crosslinking/column
purification → streptavidin capture → Tn5. The former graph reversed capture
and tagmentation, incorrectly putting Tn5 on IP chromatin. This is a repository
correction. Release/recovery is a warned MAGIC boundary; its individual
molecular products remain unrepresented.

🟢 Table4 publishes unindexed Ad1_noMX and indexed Ad2.X. The local endpoint
now has no i5, with Ad2.1_TAAGGCGA as a representative member (written primer
index TCGCCTTA). Exact graph primer inputs derive from canonical arms and the
inspected sequences. The former shared dual-index endpoint did not match
this design. Article111 states paired-end 75-base acquisition; compatible
i7 bases alone do not establish acquired index cycles.

🟢 Fill-in uses Klenow large fragment; proximity ligation uses T4 DNA ligase;
PCR uses Phusion HF. The graph names these and selects the published 300–700 bp
PAGE route. The alternative two-sided AMPure route and optional low-input
pooling/short-primer amplification are not combined with it. Candidate
library-ineligible Tn5 products have an unresolved fate: non-amplification
does not establish physical removal or experimental observation.

Deferred: full published supplementary protocol/figures and remaining primer
family rows; in situ Hi-C/ATAC predecessors and vendor originals; applicability
of a separate correction to the chosen revision; individual elution, reverse
crosslink, column, EDTA and bead-wash states; loaded-end pool/off-product
geometry; sample-specific Ad2/index acquisition. This capped pass does not
validate complete experimental chemistry or accept source absence.
