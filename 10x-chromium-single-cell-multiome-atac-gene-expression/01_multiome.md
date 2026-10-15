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

🟢 CG000338 Rev G Appendix pp. 80–81 prints the complete constant structure of both gel
bead oligos and both final libraries. The page now transcribes the GEX Read-1 arm,
16-base barcode, 12-base UMI, T30VN, dual 10-base sample indices, and the ATAC P5/barcode/
`CGCGTCTG`/s5/ME and s7/8-base-index/P7 shell exactly. Variable barcode, UMI and index
windows are structured placeholders, not inferred sequence. 🔴 The internal ATAC bridge
reagent and some reaction-mixture oligos remain proprietary.

### Guide-revision discrepancy

🟢 CG000338 Rev B Appendix p. 80 prints the gene-expression final library's 5′ P5 arm
as `AATGATACGGCGACCACCGAGATCT` before the `N10` i5 index. The ATAC gel-bead and
library lines on that same page include the following `ACAC`; Rev F and Rev G Appendix
pp. 80–81 also print the full `AATGATACGGCGACCACCGAGATCTACAC` in the
gene-expression library. 🟡 The model uses the full Rev G P5 arm, supported by the
other lines and canonical Illumina P5. 🔴 10x does not state whether Rev B's shorter
gene-expression line was a diagram error or represented a reagent variant.
