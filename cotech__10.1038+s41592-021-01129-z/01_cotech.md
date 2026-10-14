# CoTECH

**Defining paper:** Xiong et al., *Nature Methods* (2021), doi:10.1038/s41592-021-01129-z.

## Molecular boundary

🟢 Antibody-bound cells are distributed to wells containing uniquely barcoded protein-A–Tn5 adapters and a well-barcoded poly(dT) RT primer. Target-proximal DNA and mRNA therefore acquire a shared first-round well identity.

🟢 Cells are pooled and redistributed for preamplification, supplying the second combinatorial coordinate. Purified products are divided into separate chromatin and RNA libraries that remain linked by cell barcode.

## Evidence boundary

🟢 The paper establishes the barcode topology and distinct TruSeq-like DNA and Nextera-like RNA branches. Supplementary Table 1 exposes the custom adapter, RT and amplification oligos; the page currently represents their roles.

## Oligo source audit

🟢 The publisher's supplementary workbook `s41592-021-01129-z_41592_2021_1129_MOESM3_ESM.xlsx` is **Supplementary Table 1: custom primers**. Its `Custom_primers` sheet gives the PA-Tn5 mosaic end `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, phosphorylated common annealing oligo `5Phos/CTGTCTCTTATACACATCT` (modification notation as printed), template-switch oligo `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G`, IS/P2 primers, connector A/B and first-round TruSeq PCR primers.

🟢 The `PA-Tn5_barcodes` sheet lists full T5 and T7 adapter oligos for the well-specific transposome, while `barcoded_oligo_dT` lists full RT primers for each well. For example, T5-1 is `TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG`, and `RT-barcode_01` is `TCAGACGTGTGCTCTTCCGATCTAACGTGATDDDDDDDDttttttttttttttttttttttttt` as printed. The table identifies the D window as the UMI and lowercase t run as the poly(dT) capture portion.

🟡 The defining supplement does expose adapter and RT-primer bases. The public role-level page can be advanced using a selected published T5/T7/well combination; one row must not be presented as the universal barcode set. Keep the source's `D` and `+G` notation rather than replacing them with a guessed chemistry.
