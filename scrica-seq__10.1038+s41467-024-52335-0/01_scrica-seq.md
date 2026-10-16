# scRICA-seq — research notes

## Sources read

- 🟢 Jiang et al. 2024, DOI `10.1038/s41467-024-52335-0`, PMC article, Methods and supplements.
- 🟢 Supplementary Data 2 (`MOESM5`) oligo workbook, checked against the
  publisher-hosted workbook text.
- 🟢 The [Research Square v1 preprint](https://assets-eu.researchsquare.com/files/rs-3640234/v1/ab0538bc-4ec8-4a30-b79a-01c4a841cd14.pdf?c=1726297603), Methods pp. 8–11, and the [author scRCAT-seq2 pipeline](https://github.com/huyoujinlab/scRCAT-seq2), README and `run_star.sh`. The preprint's linked primer workbook and Extended Data figures have not been obtained; see `catalogue/source_reviews/` for the source checklist.

## Construct path

🟢 The RNA side is scRCAT-seq2: individual RNA/cDNA molecules receive UMI-bearing end labels, full-length cDNA is circularized, circles are amplified and repeated copies are randomly tagmented. Short reads with one UMI can then cover the TSS, TES and internal exons of one isoform. 🟢 scRICA-seq pairs this branch with scATAC-seq from the same cell or nucleus. Circularization belongs only to the RNA-derived branch.

🟢 The primary supplement publishes the complete low-throughput oligo set. The
oligo(dT) primer contains the SMART handle, a 6-base cell barcode, a 12-base UMI,
tag1 and `T30VN`; the TSO ends in three riboguanosines. The linker is
`revcomp(SMART handle) + TSO`, which supplies the homologous ends used by NEBuilder
for circle formation. The Phi29 primer is the TSO sequence itself. Separate printed
anti-tag1 and anti-TSO primer routes enrich TSS and TES fragments, followed by complete
N8-bearing N5/N7 primers. These are no longer role-only placeholders in the model.

🟢 The low-throughput Methods paragraph and Supplementary Fig. 3 show TSS and
TES enrichment as separate PCR1/PCR2 branches. Each indexed branch is then
size-separated into long and short retained libraries. The graph has four
sequencing terminals for those RNA libraries, alongside the separate ATAC
library. The intermediate PCR products have no exact reconstructed molecular
ends in the model, so their primer geometry remains explicitly unresolved.
Circularization uses NEBuilder HiFi DNA Assembly Master Mix, not an asserted
standalone ligase.

🟡 The public page uses the low-throughput route as the sequence-explicit
representative. The high-throughput route begins from 10x Multiome cDNA and has its own
printed downstream primers, but the proprietary 10x bead oligo remains covered by the
separate 10x protocol rather than being invented here.

## Published barcode-length discrepancy

🟢 Supplementary Data 2, `MOESM5`, prints `(N)6(N)12` in the low-throughput
oligo(dT) primer row. Its footnote calls `(N)6` a **16bp** cell barcode, but the
original Supplementary Fig. 3 (low-throughput molecular drawing, PDF p. 4)
labels N6 a **6bp** cell barcode. Supplementary Fig. 2 (high-throughput drawing,
PDF p. 3) instead labels N16 as 16bp. 🟡 The model uses six bases for the
low-throughput RNA oligo and library because the printed sequence and matching
figure agree. The physically paired ATAC fraction has no molecular cell-barcode
segment in this low-throughput view. The table footnote appears to have borrowed
the high-throughput length; no author correction has been verified.

## High-throughput read-layout limit

🟢 The author pipeline README assigns R1 to cDNA and R2 to cell barcode plus UMI.
Its TSS example extracts 28 bases, consistent with the 10x 16-base barcode and
12-base UMI. Its TES example instead extracts 18 bases, while `run_star.sh`
configures a 16-base barcode and 12-base UMI for TES. 🟡 The published
high-throughput TES processing layout is therefore unresolved from these code
examples. The public schematic represents the separate low-throughput branch;
the preprint's linked primer workbook and Extended Data figures still need to be read.
