# 10x Chromium Single Cell Multiome ATAC + Gene Expression — research notes

## 1. What it is

The commercial Multiome assay jointly measures open chromatin and polyadenylated RNA from the same nucleus. Bulk nuclear tagmentation precedes GEM partitioning; the GEM barcode is then incorporated into both the ATAC and gene-expression products.

## 2. Sources read

- 🟢 10x Genomics, *Chromium Next GEM Single Cell Multiome ATAC + Gene Expression User Guide*, CG000338 Rev G (2024): [current support page](https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide).
- 🟢 The guide's workflow, library diagrams and sequencing section establish bulk tagmentation, GEM barcoding, the two library branches, and their read roles.

## 3. Construct-changing steps

1. Tn5 tags accessible chromatin in intact nuclei.
2. Nuclei, gel beads and reaction mix form GEMs.
3. One gel-bead barcode is transferred to both ATAC fragments and oligo-dT-primed cDNA; RNA products also receive a UMI.
4. Cleanup separates the material into modality-specific amplification and library construction branches.

## 4. Evidence boundary

🟢 Barcode and UMI roles and the two final read layouts are printed in CG000338. 🔴 Complete gel-bead and bridge-oligo bases are proprietary; the schematic uses inferred placeholders only for those undisclosed regions.
