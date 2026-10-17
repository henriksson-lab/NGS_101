# MIRACL-seq — droplet overloading for rare-cell RNA-seq

## Sources read

- 🟢 Drokhlyansky *et al.*, “The human and mouse enteric nervous system at single cell resolution,” *Cell* (2020), DOI `10.1016/j.cell.2020.08.003`, PMCID `PMC8358727` (defining paper and STAR Methods).
- 🟢 10x Genomics GemCode Single Cell 3′ v2 chemistry, the commercial chemistry named by the paper.
- Bounded pass 2026-10-17 (round breadth-first-13): published PMC8358727 JATS STAR Methods lines 1074–1103 (droplet RAISIN RNA-seq; MIRACL-seq overloading), 1127–1145 (Cell Ranger v2.0 processing), 1595–1625 (heterotypic-droplet removal), Results 442–470, Key Resources 2963–2972 (10X Chromium Single Cell 3′ Kit Cat#120237; NextSeq 500 high output v2 75-cycle kit); 10x CG000108 Rev A technical note (in the `10x-chromium-single-cell-3-ge-v2` source bundle) lines 1–36 and 60–139 for the v2 gel-bead oligo, TSO, cDNA primers, Read 2 adapter, sample-index PCR primers and final library. Figures S1–S7 and Tables S1–S6 remain PMC download-gated (unread).

## Molecular interpretation

- 🟢 MIRACL-seq means **Mining Rare Cells**. It applies droplet-based RAISIN RNA-seq to unsorted material while deliberately loading more nuclei than the nominal single-cell regime.
- 🟢 The defining experiments tested overloading coefficients from 2× to 8×. Rare cells are recovered more efficiently; because they are rare, rare-cell/rare-cell homotypic multiplets remain unlikely. Heterotypic contaminated droplets are identified computationally and removed.
- 🟢 Molecular capture itself uses the GemCode Single Cell Platform with v2 Gel Beads, Chip and Library kits: GEM lysis, bead-barcoded reverse transcription, amplification, shearing, 5′ adapter attachment and sample indexing.
- 🟢 Droplet RAISIN RNA-seq loads 7,000 RAISINs per channel (resuspended in ST buffer); MIRACL-seq overloads this 2×–8× (Methods lines 1079–1080, 1098). The per-experiment coefficients are in Table S2, which was not readable.
- 🟡 The paper does not list the kit steps beyond “amplification, shearing and 5′ adaptor and sample index attachment … following the manufacturer’s protocol”. CG000108 describes template switching, cDNA PCR, enzymatic fragmentation, size selection, Read 2 adapter ligation and sample-index PCR; which CG00052 user-guide revision was followed is not stated, so the graph keeps those kit intermediates behind a warned boundary and only the sample-index PCR endpoint is drawn.
- 🟡 MIRACL-seq does not introduce a new oligo architecture. The schematic therefore reuses the standard GemCode v2 endpoint and represents the named innovation at the partition-loading step.

## Final readout

Read 1 contains the 16-nt GemCode cell barcode and 10-nt UMI; Read 2 reports transcript-derived sequence and i7 identifies the library. 🔴 The paper prints neither the acquired read configuration for the droplet libraries (only “Illumina NextSeq 500” and Cell Ranger v2.0) nor the sample-index members used; an earlier model value of 26 Read 1 cycles was a vendor-style assumption and has been removed. A barcode can describe several nuclei in an overloaded GEM, so molecular barcoding alone does not resolve multiplets.
