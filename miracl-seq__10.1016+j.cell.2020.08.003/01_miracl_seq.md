# MIRACL-seq — droplet overloading for rare-cell RNA-seq

## Sources read

- 🟢 Drokhlyansky *et al.*, “The human and mouse enteric nervous system at single cell resolution,” *Cell* (2020), DOI `10.1016/j.cell.2020.08.003`, PMCID `PMC8358727` (defining paper and STAR Methods).
- 🟢 10x Genomics GemCode Single Cell 3′ v2 chemistry, the commercial chemistry named by the paper.

## Molecular interpretation

- 🟢 MIRACL-seq means **Mining Rare Cells**. It applies droplet-based RAISIN RNA-seq to unsorted material while deliberately loading more nuclei than the nominal single-cell regime.
- 🟢 The defining experiments tested overloading coefficients from 2× to 8×. Rare cells are recovered more efficiently; because they are rare, rare-cell/rare-cell homotypic multiplets remain unlikely. Heterotypic contaminated droplets are identified computationally and removed.
- 🟢 Molecular capture itself uses the GemCode Single Cell Platform with v2 Gel Beads, Chip and Library kits: GEM lysis, bead-barcoded reverse transcription, amplification, shearing, 5′ adapter attachment and sample indexing.
- 🟡 MIRACL-seq does not introduce a new oligo architecture. The schematic therefore reuses the standard GemCode v2 endpoint and represents the named innovation at the partition-loading step.

## Final readout

Read 1 contains the 16-nt GemCode cell barcode and 10-nt UMI; Read 2 reports transcript-derived sequence and i7 identifies the library. A barcode can describe several nuclei in an overloaded GEM, so molecular barcoding alone does not resolve multiplets.
