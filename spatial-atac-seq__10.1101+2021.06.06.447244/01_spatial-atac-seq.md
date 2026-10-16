# Spatial-ATAC-seq — source notes

## What it is

Spatial-ATAC-seq performs Tn5 tagmentation in a fixed tissue section and then ligates two
orthogonal microfluidic barcode sets. The ordered A×B combination records the tissue pixel
for each accessible-chromatin fragment.

## Sources read

- 🟢 Deng et al., *Nature* (2022),
  [doi:10.1038/s41586-022-05094-1](https://doi.org/10.1038/s41586-022-05094-1), is the
  defining published paper; the catalogue directory retains its earlier defining-preprint
  DOI [10.1101/2021.06.06.447244](https://doi.org/10.1101/2021.06.06.447244).
- 🟢 Farzad et al., *Nature Protocols* (2024),
  [doi:10.1038/s41596-024-01013-y](https://doi.org/10.1038/s41596-024-01013-y), gives the
  detailed shared Spatial-ATAC-seq/Spatial-CUT&Tag workflow. Its Supplementary Tables 2–3
  print the linker, barcode, PCR and custom sequencing-primer oligos.
- 🟢 The authors' `dyxmvp/Spatial_ATAC-seq` processing repository uses the same raw-read
  parser as Spatial-CUT&Tag: barcode B is raw Read 1 bases 1–8, barcode A is bases 39–46,
  and raw Read 2 is genomic sequence.

## Molecular reconstruction

🟢 The protocol figure gives the final order as P5, custom Read 1 site, barcode B,
ligation linker 2, barcode A, ligation linker 1, genomic DNA, S7, i7 and P7. The page uses
the shared exact spatial-epigenome constructor because the barcoding and final library are
the same for the two assays; only the targeting reaction differs.

🟡 The displayed read strand is derived by reverse-complementing the linker oligos as
ordered. The first and second barcode oligos determine the two covalent ligation boundaries
inside those linkers, which the model validates and marks.

🟢 The defining Spatial-ATAC Methods print the Tn5 loading oligo `Tn5ME-A` as
`/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG`. Farzad et al. 2024
Supplementary Table 1 prints a different `Tn5ME-A` for Spatial-CUT&Tag; the shared
spatial linker and barcode architecture does not make those loading oligos interchangeable.

## Read layout

🟢 Paired-end 150 sequencing was used. The custom Read 1 primer enters barcode B first,
then linker 2, barcode A and linker 1 before accessible genomic DNA. The i7 index is read
from the S7 side, and Read 2 enters genomic DNA from the opposite end.
