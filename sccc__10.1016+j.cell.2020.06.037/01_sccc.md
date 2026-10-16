# scCC

**Defining paper:** Moudgil et al., *Cell* (2020), doi:10.1016/j.cell.2020.06.037. The open author manuscript and methods were read from PMCID PMC7510185.

## Molecular boundary

🟢 A self-reporting piggyBac transposon lacks its downstream polyadenylation signal, so reporter transcription continues through the terminal repeat into the genomic insertion flank.

🟢 The 10x first-strand product is split. The scCC branch is specifically amplified and circularized, moving the transcript's cell barcode and UMI next to the transposon–genome junction before linear-DNA removal, shearing and junction library construction.

## Evidence boundary

🟢 The paper's Table S4 defines custom adapters and primers. The page preserves their roles and links the source rather than treating this historical Nextera Mate Pair implementation as a universal kit geometry.

## Oligo source audit

🟢 The author preprint's supplementary PDF (`538553_media-1.pdf.txt`, Supplemental Table 4) is retrievable even though the published Cell supplement `NIHMS1608211-supplement-1.pdf` is gated at PMC. It prints the scCC library oligos verbatim, including modifications and a variable sample index. The two annealed adapter components are `scCC_P5_adapter`, `AATGATACGGCGACCACCGAGATCTTCACTCATTCCACACGACTCCTTGCCAGTCTC*T`, and `scCC_P7_adapter`, `/5Phos/GAGACTGGCAAGTACACGTCGCACTCACCATGA[index]ATCTCGTATGCCGTCTTCTGCTTG`. Their paired geometry should be computed before an exact junction is shown.

🟢 The same table gives the final PCR primers `scCC_P5_primer`, `AATGATACGGCGACCACCGAGATC`, and `scCC_P7_primer`, `CAAGCAGAAGACGGCATACGAGAT`; the custom sequencing primers are `scCC_PB_CustomRead2`, `CGTGTAGGGAAAGAGTGTGCGTCAATTTTACGCAGACTATCTTTCTAG`, and `scCC_CustomIndex1`, `GAGACTGGCAAGTACACGTCGCACTCACCATGA`. Capture primers include `Bio_Illumina_Seq1_scCC_10X_3xPT`, `/5Phos/ACACTCTTTCCC/iBiodT/ACACGACGCTCTTCCGA*T*C*T`, and `Bio_Long_PB_LTR_3xPT`, `/5Phos/GCGTCAATTTTACGCAGAC/iBiodT/ATCTTTC*T*A*G`.

🟢 The final publisher supplement has now been recovered directly as `1-s2.0-S009286742030814X-mmc1.pdf`. Its Table S4 (extracted lines 443–539) independently prints all eight scCC oligos shown on the page with the same bases and modifications. The variable `[index]` remains a source-declared symbolic window.

## Final-source workflow audit (2026-10-16)

🟢 The final companion `mmc2.pdf`, “Molecular biology of single cell calling cards”
(lines 183–199), and the article Methods (lines 220–234) establish the ordered split:
initial droplet RT omits TSO; the expression aliquot receives TSO and further RT;
the scCC aliquot receives a junction-enrichment PCR using the two phosphorylated,
internally biotinylated primers before circularization. After linear-DNA removal,
shearing and streptavidin capture, end repair, A-tailing and custom-adapter ligation
precede final PCR. The graph now includes the formerly omitted enrichment PCR,
expression-specific TSO restoration, end repair, A-tailing and adapter ligation.
These omissions were repository modelling gaps.

🟢 Article sequencing/analysis Methods (line 300) specifies standard Illumina BP10
for Read 1 and BP14 for Index 2, plus custom Read 2 and Index 1 primers. Read 1
contains the cell barcode and UMI; Read 2 begins with `GGTTAA` then enters genomic
sequence. Exact primer placement and run-cycle counts remain model/review tasks,
rather than claims that those oligos are unpublished.

The inherited Nextera Mate Pair original reference guide (15035209 v02, January 2016)
explicitly uses exonuclease to remove linear DNA and performs on-bead end repair,
A-tailing and adapter ligation. Its generic tagmentation and generic adapters are
not imported into scCC: the defining scCC Methods supplies the custom changes.
