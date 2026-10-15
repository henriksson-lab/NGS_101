# Multi-CUT&Tag — research notes

## Sources read

- 🟢 Gopalan et al. 2021, DOI `10.1016/j.molcel.2021.09.019`, PMC full text and STAR Methods.
- 🟢 Supplementary Table 1 in the bioRxiv supplement (`2021.04.27.441642_media-1.pdf`): complete target-coded A/B oligos, Tn5MErev and custom sequencing primers.

## Exact adapter design

🟢 Four A-family and four B-family top oligos are printed in Table S1. Each contains a custom sequencing-primer site, an 8-nt antibody barcode, a conserved connector and the 19-nt Tn5 mosaic end. The common bottom strand is `5′-/5Phos/CTGTCTCTTATACACATCT-3′`; it pairs only to the mosaic end, leaving the barcode-containing part single-stranded in the loaded transposome.

🟢 A-family structure: `TCGTCGGCAGCGTCTCCACGC` + 8-nt antibody barcode + `GCGATCGAGGACGGC` + `AGATGTGTATAAGAGACAG`.

🟢 B-family structure: `GTCTCGTGGGCTCGGCTGTCCCTGTCC` + 8-nt antibody barcode + `CACCGTCTCCGCCTC` + `AGATGTGTATAAGAGACAG`.

🟢 Each antibody is pre-coupled to pA–Tn5 loaded with its matched barcoded adapters. Magnesium activation transfers those adapters to nearby chromatin. A fragment may carry the same antibody code at both ends (homogeneous) or different codes at its ends (mixed, interpreted as co-localization).

## Completion and sequencing

🟢 Bulk libraries receive a 72 °C, 5-minute gap fill, then 17 PCR cycles with indexed i5 and i7 primers. The source run is Read 1 72 cycles, Read 2 72 cycles, Index 1 8 cycles and Index 2 8 cycles. Both genomic reads begin with 8-nt antibody barcode + 15-nt connector + 19-nt mosaic end; the paper trims these first 42 bases before alignment.

🟢 The custom sequencing primers are Read1 `TCGTCGGCAGCGTCTCCACGC`, Read2 `GTCTCGTGGGCTCGGCTGTCCCTGTCC`, Index1 `GGACAGGGACAGCCGAGCCCACGAGAC`, and Index2 `GCGTGGAGACGCTGCCGACGA`.

🟢 For single-cell Multi-CUT&Tag, the target-coded fragments enter the 10x Chromium Single Cell ATAC workflow. The paper uses one GEM amplification cycle and 16 Index 2 cycles. The named [10x CG000209 v1.1 guide](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevG.pdf) shows the gel-bead oligo's P5–16-base barcode–Read 1N architecture (Appendix A), and the vendor [Cell Ranger ATAC software](https://github.com/10XGenomics/cellranger-atac) distributes its barcode whitelist (`lib/python/barcodes/737K-cratac-v1.txt`). 🔴 The paper does not assign a particular whitelist barcode to a particular cell. The public schematic leaves that variable molecule abstract.
