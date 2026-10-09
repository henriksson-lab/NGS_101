# PairedTag

## 1. What it is

PairedTag jointly profiles a chosen histone modification and the transcriptome from the same single nuclei using antibody-targeted tagmentation and shared combinatorial cell barcodes.

## 2. Authoritative source

🟢 Zhu et al., *Joint profiling of histone modifications and transcriptome in single cells from mouse brain*, Nature Methods (2021), DOI `10.1038/s41592-021-01060-3`.

The PMC full text and Supplementary Table 1 (primer sequences) were read.

## 3. Construct-changing operations

1. 🟢 Antibody-bound pA–Tn5 transfers a barcoded mosaic-end adapter to the selected
   histone-mark DNA. In the same first-round well, barcoded oligo(dT) and random primers
   reverse-transcribe nuclear RNA.
2. 🟢 Nuclei are pooled, split across 96 round-2 ligation barcodes, blocked, pooled again,
   and split across 96 round-3 barcode/UMI oligos.
3. 🟢 TdT adds a dC tail; Anchor-FokI-GH supports linear amplification, followed by common
   PCR of the mixed DNA/cDNA material.
4. 🟢 The amplified material is split. DNA is cut with SbfI/FokI and receives a P5 adapter
   by ligation; RNA-derived material is cut with NotI and receives an N5 end by Tn5.
5. 🟢 Both final libraries use P7 indexing PCR and were sequenced as Read 1 + 7-base
   Index 1 + Read 2.

## 4. Key published oligos

- `pMENTs`: `/5Phos/CTGTCTCTTATACACATCTddC`
- `AdaptorA`: `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`
- `Linker-R02`: `CGAATGCTCTGGCCTCTCAAGCACGTGGAT`
- `Linker-R03`: `GGTCTGAGTTCGCACCGAAACATCGGCCAC`
- `P5-FokI`: `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`
- `PA-F`: `CAGACGTGTGCTCTTCCGATCT`
- `PA-R`: `AAGCAGTGGTATCAACGCAGAGT`

The plate contains all 12 × 96 × 96 barcode choices; the schematic uses role-labelled
six-base representatives rather than presenting one cell identity as universal.

## 5. Final read layouts

Read 1 maps the histone-mark DNA or transcript cDNA. Read 2 traverses the 10-base UMI and
the three nuclear barcodes; the paper places the first bases of BC3, BC2 and BC1 around
positions 10–13, 47–50 and 84–87, respectively. Index 1 contains the PCR sample index.
