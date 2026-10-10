# scRICA-seq — research notes

## Sources read

- 🟢 Jiang et al. 2024, DOI `10.1038/s41467-024-52335-0`, PMC article, Methods and supplements.

## Construct path

🟢 The RNA side is scRCAT-seq2: individual RNA/cDNA molecules receive UMI-bearing end labels, full-length cDNA is circularized, circles are amplified and repeated copies are randomly tagmented. Short reads with one UMI can then cover the TSS, TES and internal exons of one isoform. 🟢 scRICA-seq pairs this branch with scATAC-seq from the same cell or nucleus. Circularization belongs only to the RNA-derived branch.
