# CROP-seq

## 1. What it is

CROP-seq makes a pooled CRISPR guide detectable in the same 3′ single-cell RNA-seq library
that measures the cell's transcriptional response. The defining implementation used Drop-seq.

Defining paper: Datlinger *et al.*, “Pooled CRISPR screening with single-cell transcriptome
readout,” *Nature Methods* (2017), doi:[10.1038/nmeth.4177](https://doi.org/10.1038/nmeth.4177).

## 2. Sources read

- 🟢 PMC full text (PMCID PMC5334791), especially Figure 1e, “Direct detection of gRNAs
  from single-cell transcriptome data,” and Online Methods.
- 🟢 Supplementary Protocol: pooled guide cloning, Drop-seq cDNA preparation, Nextera XT
  enrichment, the oligo-order table, and sequencing configuration.
- 🟢 Supplementary Table 1: vector construction, duplication validation and library-QC oligos.
- 🟢 Supplementary Data: the 10,214-bp CROPseq-Guide-Puro / Addgene #86708 GenBank map.

## 3. The vector operation

🟢 The paper places the hU6–guide cassette in the 3′ long terminal repeat of a
self-inactivating LentiGuide-Puro derivative. A packaged transfer RNA has `R–U5` at its
5′ end and `U3–R` at its 3′ end. During reverse transcription, the 3′ U3 rebuilds both
integrated `U3–R–U5` LTRs. The cassette inserted with that U3 is therefore duplicated.

This is represented by `lib/lentivirus.py`: callers supply U3 once, and the integrated
provirus constructor derives both copies. Equality is a construction rule, not a repeated test.

🟢 The downstream copy becomes part of the polyadenylated puromycin-resistance RNA
polymerase II transcript, while hU6 independently produces the functional,
non-polyadenylated RNA polymerase III guide. One DNA cassette creates two overlapping RNAs.

## 4. Guide cloning junction

🟢 The Supplementary Protocol gives each pooled 74-nt ssDNA oligo as:

```text
5'-TGGAAAGGACGAAACACCG-[20-nt guide]-GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC-3'
```

The 5′ prefix provides 18 bases of hU6 homology and adds G as the constant first base for
U6 transcription. The 35-base postfix is homologous to the guide backbone. The oligo pool
is assembled into the 8,333-bp BsmBI-digested backbone by Gibson/NEBuilder assembly.
For individual guides, the paper instead gives phosphorylated `CACCG(N)20` /
`AAAC(N)20C` oligos and BsmBI ligation.

🟢 The deposited empty 10,214-bp vector map contains a 1,885-bp filler between hU6 and the
guide scaffold; the BsmBI sites remove it. That filler belongs to the empty cloning backbone,
not a guide-bearing viral construct.

## 5. Single-cell library

🟢 Cells and Drop-seq barcoded beads are co-encapsulated. Polyadenylated RNAs are captured
by bead oligo-dT, reverse-transcribed with template switching, SMART-PCR amplified,
tagmented with Nextera XT, and selectively enriched with New-P5-SMART and N70x primers.
The authors use more cDNA at tagmentation to favor larger fragments and increase the chance
that a fragment retains the guide.

🟢 Ordered oligos printed in the protocol include:

```text
TSO:                 AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
SMART PCR primer:    AAGCAGTGGTATCAACGCAGAGT
New-P5-SMART hybrid: AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C
Custom Read 1:       GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
```

🟢 The run is 20 bases Custom Read 1, 8 bases Index 1, and 64 bases Read 2. Read 1 is
the 12-base cell barcode plus 8-base UMI. Read 2 supplies transcript sequence and identifies
guide-bearing molecules when it crosses the guide.

## 6. Separate pooled-library QC readout

🟢 The protocol also gives a one-PCR Illumina library for guide-pool QC or classical
pooled-screen genomic-DNA readout. Its staggered i5-side primer ends after the `CACCG`
cloning junction, so Read 1 starts at guide base 1. This is not the cell-barcoded Drop-seq
library and is intentionally separate from the main schematic.

## 7. Boundary of this page

This page documents original CROP-seq and its Drop-seq readout. Later dedicated guide-capture,
10x Feature Barcode, and multi-perturbation chemistries are separate protocols.
