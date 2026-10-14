# scRibo-seq

## Purpose

scRibo-seq couples ribosome-footprint generation to one-pot small-RNA library construction in individual wells.

## Source read

🟢 Van Insberghe et al., *Single-cell Ribo-seq reveals cell cycle-dependent translational pausing*, Nature (2021), [doi:10.1038/s41586-021-03887-4](https://doi.org/10.1038/s41586-021-03887-4), main text and Methods pp. 5–6. Supplementary Table 1 in `41586_2021_3887_MOESM3_ESM.xlsx` was read directly for the oligo-ordering strings.

## Oligos from Supplementary Table 1

🟢 `OMV630_miRNA4_3App`: `/5rApp/TGGAATTCTCGGGTGCCAAGG/3ddC/`; the 5′ adenylation and 3′ dideoxy-C are part of the ordered adapter.

🟢 `OMV572_miRNA4_RT`: `GCCTTGGCACCCGAGAATTCCA`.

🟢 `OMV632_miRNA5_5A_10U`: `rCrUrArCrArCrGrArCrGrCrUrCrUrUrCrCrGrArUrCrUr(N1:25252525)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)`. This retains the source's modified RNA and random-base notation; the ten random positions follow the fixed 5′ RNA sequence.

🟢 The well-indexed forward primer family is printed as individual `miRv6-PCR_F-cbc` rows (for example A1: `AATGATACGGCGACCACCGAGATCTACACGATTCACGACACACTCTTTCCCTACACGACGCTCTTCCG`). The plate-indexed reverse family is printed as `RPI` rows (for example RPI24: `CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`). The XLSX is the authority for each chosen well/plate member; these are examples, not universal index choices.

## Construct-changing path

🟢 Single cells are lysed under mineral oil and digested with MNase. Proteinase K plus EGTA/EDTA stops digestion. T4 PNK end repair precedes overnight ligation of the pre-adenylated `OMV630_miRNA4_3App` 3′ adapter with T4 RNA ligase 2 truncated KQ.

🟢 `OMV572_miRNA4_RT` is annealed before `OMV632_miRNA5_5A_10U` is joined with T4 RNA ligase 1. Maxima H-minus RT generates cDNA; well-specific `miRv6-PCR_F-cbc` and an `RPI-` primer index each single-cell library in ten PCR cycles.

🟢 Pooled libraries are PAGE-selected at 175–185 bp, corresponding to an approximately 30–40-nt insert. Reads retained computationally are 30–45 nt after adapter/UMI removal.

## Read layout

🟢 NextSeq 500 sequencing uses Read 1 75 cycles, i7 six cycles (plate index), and i5 ten cycles (cell index). The first ten bases of Read 1 are the UMI; the footprint follows. The paper trims `TGGAATTCTCGGGT` from the Read 1 3′ end.

🟡 The schematic's final duplex and sequencing-primer landing sites still require reconciliation with the complete ordering strings; the note above records the source oligos exactly, while the page retains its existing role-level final construct.
