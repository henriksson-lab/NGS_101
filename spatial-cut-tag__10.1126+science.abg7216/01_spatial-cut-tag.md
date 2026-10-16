# Spatial-CUT&Tag / hsrChST-seq — source notes

## What it is

Spatial-CUT&Tag is the published name of the method introduced as
**high-spatial-resolution chromatin modification state profiling by sequencing
(hsrChST-seq)** in the authors' preprint. It is one protocol, not a family: antibody-tethered
pA–Tn5 tags histone-mark-associated DNA in a fixed tissue section, then two perpendicular
microfluidic ligations add an A×B spatial address.

## Sources read

- 🟢 Deng et al., *Science* (2022),
  [doi:10.1126/science.abg7216](https://doi.org/10.1126/science.abg7216), is the defining
  published paper.
- 🟢 Deng et al., bioRxiv (2021),
  [doi:10.1101/2021.03.11.434985](https://doi.org/10.1101/2021.03.11.434985), supplies the
  original hsrChST-seq name and describes antibody binding, pA–Tn5 activation, A and B
  ligations, crosslink reversal and library PCR.
- 🟢 Farzad et al., *Nature Protocols* (2024),
  [doi:10.1038/s41596-024-01013-y](https://doi.org/10.1038/s41596-024-01013-y), is the
  authors' detailed protocol for Spatial-CUT&Tag and Spatial-ATAC-seq. Supplementary
  Table 1 prints the CUT&Tag loading oligos and primers; Table 2 is the distinct
  Spatial-ATAC set, and Table 3 prints the shared spatial barcodes.
- 🟢 The authors' `dyxmvp/spatial-CUT-Tag` processing repository was read to verify the raw
  read coordinates: barcode B is raw Read 1 bases 1–8, barcode A is bases 39–46, and raw
  Read 2 is genomic sequence.

The 2021 bioRxiv PDF appends Supplementary Fig. S1 and other figure captions, but no
ordered oligo table or full wet-lab recipe. Its Data and materials availability section
links `dyxmvp/hsrChST-seq`, which did not resolve when checked; the accessible
`dyxmvp/spatial-CUT-Tag` repository may be a successor, but equivalence is unverified.
The 2024 supplementary PDF MOESM1 contains only a gene-score mapping figure; the
ordered oligos below are in the separate MOESM2 workbook. The defining 2022 Science
supplementary DOCX remains unavailable, so the 2024 sequences describe a later
published implementation and have not been independently established as the exact
sequences used in every 2022 experiment.

## Molecular reconstruction

🟢 The protocol figure gives the final order as P5, custom Read 1 site, barcode B,
ligation linker 2, barcode A, ligation linker 1, genomic DNA, S7, i7 and P7.

🟡 The page draws the sequenced strand by reverse-complementing the two linker oligos as
ordered. This is fixed by the barcode oligos themselves: `BC2_R1_1` contributes the first
15 bases of linker 2, `Tn5ME-A_BC1_1` contributes its second 15 bases, and their join is
the barcode-B ligation. The same first-barcode oligo contributes the 19-base barcode-side
part of linker 1; its join to the 15-base genomic-side part is the barcode-A ligation.
Both boundaries are consequently computed and marked rather than placed by eye.

🟢 The 50 first-barcode sequences and 50 second-barcode sequences are public. The page
shows one representative oligo from each plate and records the two 8-base regions as
separate parts of one combinatorial spatial coordinate.

🟢 Farzad et al. Supplementary Table 1 prints the CUT&Tag-specific Tn5ME-A loading
oligo, `/5Phos/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`. It differs from the
Spatial-ATAC loading oligo in Deng et al. (2022) Methods; the page displays the
CUT&Tag value and uses the shared linker and barcode construct for the sequenced strand.

## Read layout

🟢 The authors used HiSeq 4000 paired-end 150 sequencing. The custom Read 1 primer enters
barcode B first, then linker 2, barcode A and linker 1 before genomic sequence. The
Nextera-compatible index primer reads the i7 sample index, and the opposite Nextera read
primer enters genomic DNA from the S7 end.
