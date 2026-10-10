# sci-CAR-seq — research notes

## Sources read

- 🟢 Cao et al. 2018, DOI [10.1126/science.aau0730](https://doi.org/10.1126/science.aau0730), PMC6571013, main text and Methods.

## Molecular order

🟢 Intact nuclei first receive a well-indexed poly(dT) RT primer carrying a UMI, then a
well-indexed Tn5 transposome. Nuclei are pooled and redistributed by FACS. After cDNA
second-strand synthesis, each well is lysed and split: the RNA aliquot receives an
unindexed Tn5 site and a well-indexed PCR; the ATAC aliquot is amplified from its original
barcoded Tn5 ends with well-indexed PCR primers. RNA and ATAC libraries are pooled and
sequenced separately.

## Identifier semantics

The first RT and Tn5 indexes have different oligo sequences but encode the same first-well
coordinate. The RNA and ATAC PCR indexes likewise encode the same second/FACS well. The
two-coordinate lookup links both library types to one cell; no individual barcode segment
is physically shared between modalities.
