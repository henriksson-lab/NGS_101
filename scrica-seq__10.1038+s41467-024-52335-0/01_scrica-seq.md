# scRICA-seq — research notes

## Sources read

- 🟢 Jiang et al. 2024, DOI `10.1038/s41467-024-52335-0`, PMC article, Methods and supplements.
- 🟢 Supplementary Data 2 (`MOESM5`) oligo workbook, checked against the
  publisher-hosted workbook text.

## Construct path

🟢 The RNA side is scRCAT-seq2: individual RNA/cDNA molecules receive UMI-bearing end labels, full-length cDNA is circularized, circles are amplified and repeated copies are randomly tagmented. Short reads with one UMI can then cover the TSS, TES and internal exons of one isoform. 🟢 scRICA-seq pairs this branch with scATAC-seq from the same cell or nucleus. Circularization belongs only to the RNA-derived branch.

🟢 The primary supplement publishes the complete low-throughput oligo set. The
oligo(dT) primer contains the SMART handle, a 16-base cell barcode, a 12-base UMI,
tag1 and `T30VN`; the TSO ends in three riboguanosines. The linker is
`revcomp(SMART handle) + TSO`, which supplies the homologous ends used by NEBuilder
for circle formation. The Phi29 primer is the TSO sequence itself. Separate printed
anti-tag1 and anti-TSO primer routes enrich TSS and TES fragments, followed by complete
N8-bearing N5/N7 primers. These are no longer role-only placeholders in the model.

🟡 The public page uses the low-throughput route as the sequence-explicit
representative. The high-throughput route begins from 10x Multiome cDNA and has its own
printed downstream primers, but the proprietary 10x bead oligo remains covered by the
separate 10x protocol rather than being invented here.
