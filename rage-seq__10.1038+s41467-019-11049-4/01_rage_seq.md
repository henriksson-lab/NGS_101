# RAGE-Seq — research notes

Defining source: Singh *et al.*, [10.1038/s41467-019-11049-4](https://doi.org/10.1038/s41467-019-11049-4).

## Evidence used

- 🟢 RAGE-Seq combines short-read expression profiling with targeted long-read sequencing of full-length TCR and BCR transcripts.
- 🟢 Biotinylated capture probes enrich receptor cDNAs from barcoded droplet libraries.
- 🟢 Long reads retain the cell-identifying sequence and resolve full receptor nucleotide sequence and splice isoforms.
- 🟢 Before capture, and again after each of two Roche–NimbleGen capture rounds, the paper amplifies the full-length cDNA with TSO primer `AAGCAGTGGTATCAACGCAGAGT` and R1 primer `CTACACGACGCTCTTCCGATCT`. These sequences were present in the original Methods and are now shown on the schematic page.
- 🟢 The capture design targeted 644 exons spanning functional human immunoglobulin/TCR V, J and constant genes. Supplementary Table 7 gives the target-gene list.
- 🔴 Roche–NimbleGen designed and synthesized the SeqCap RNA Choice probe pool, but the paper and supplement do not provide its individual probe sequences.
- 🔴 Long-read libraries used Oxford Nanopore SQK-LSK108, or SQK-LSK308 for one sample. Those commercial motor-adapter sequences are not disclosed by the paper or kit documentation.

## Protocol boundary

The page covers targeted receptor capture and the long-read branch, not generic 10x gene-expression construction.

The matched short-read branch inherits the corresponding droplet scRNA-seq chemistry. The source-exact RAGE-Seq additions are the two capture-flanking amplification primers and the target definition; the capture-probe and Nanopore-adapter gaps are kept explicit.
