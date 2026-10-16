# scGET-seq — source notes

Defining paper: Tedesco et al., *Nature Biotechnology* (2022), [doi:10.1038/s41587-021-01031-1](https://doi.org/10.1038/s41587-021-01031-1).

- 🟢 Supplementary methods and oligo table were read.
- 🟢 Ordinary Tn5 and HP1α-fused Tn5 (TnH) carry distinct 8-nt adaptor barcodes, separating accessible and H3K9me3-associated insertions.
- 🟢 Tagged nuclei are processed with the 10x Chromium single-cell ATAC workflow, adding cell identity in droplets.
- 🟢 Supplementary Table 10 prints Tn5ME-A.1 as `TCGTCGGCAGCGTC-TCCGATCT-CGTACTAG-AGATGTGTATAAGAGACAG` and TnHME-A.4 as `TCGTCGGCAGCGTC-TCCGATCT-TAAGGCGA-AGATGTGTATAAGAGACAG`. The page separates S5, connector, modality barcode and mosaic end explicitly.
- 🟢 Supplementary Code MOESM11, Step 1 and Step 4, uses eight inline barcodes for demultiplexing: TnH `TAAGGCGA`, `GCTACGCT`, `AGGCTCCG`, `CTGCGCAT`; ordinary Tn5 `CGTACTAG`, `TCCTGAGC`, `TCATGAGC`, `CCTGAGAT`. These agree with Supplementary Table 10 and verify the computational modality assignment. The code's `S:AGATATATATAAGGAGACAG` filter is a read pattern, not evidence for a supplied primer.
- 🟢 The adapters are shown loaded with the common 5′-phosphorylated mosaic-end reverse oligo, leaving their S5–connector–barcode arms single stranded, then transferred onto their respective genomic targets.
- 🟡 The final panel uses TnHME-A.4 as one representative of the printed barcode families. Tn5ME-A.1 is shown in parallel through transfer; its downstream Chromium processing has the same geometry.
- 🟡 The page exposes verified genomic Read 1, Index 1 and Read 2 sites.
- 🔴 The proprietary full 10x gel-bead oligo, dedicated cell-barcode sequencing primer and SI-PCR primers remain unresolved from the checked sources. The defining article and a versioned 10x guide are still unavailable, so absence in the publication is not confirmed; their exact sequences and binding geometry are not reconstructed.
