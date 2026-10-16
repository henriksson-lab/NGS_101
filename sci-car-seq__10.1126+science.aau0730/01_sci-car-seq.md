# sci-CAR-seq — research notes

## Sources read

- 🟢 Cao et al. 2018, DOI [10.1126/science.aau0730](https://doi.org/10.1126/science.aau0730), PMC6571013, main article workflow and figure legend. Its detailed Materials and Methods are in the currently inaccessible supplementary PDF and have not been read.
- 🟢 Amini et al. 2014, DOI [10.1038/ng.3119](https://doi.org/10.1038/ng.3119), Online Methods lines 116–122 and original publisher Supplementary Table 4 (`ng.3119_41588_2014_BFng3119_MOESM2_ESM.xlsx.txt`), correction notice lines 1–6 and oligo rows 9–13. This is the cited indexed-Tn5 predecessor; it does not identify sci-CAR's own oligo members.

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

## Oligo source review

🟢 The accessible article and Methods identify four structural oligo families:
well-indexed poly(dT)/UMI RT primers, first-well barcoded ATAC transposomes,
unindexed cDNA transposomes, and well-indexed PCR primers for each modality.
The RNA and ATAC oligos in a well encode the same coordinate but have distinct
sequences. 🔴 The article does not print their ordered bases or barcode lists;
the PMC supplementary PDF and workbook are blocked by the download gate in the
source manifest. The schematic's 10-base barcode and 8-base UMI tokens are
representative placeholders, not verified oligo lengths or orderable sequences.
The exact first-round and PCR index arms, transposome mosaic ends, and
sequencing-primer sites remain open until those supplements are obtained.

### Published predecessor correction

🟢 Amini 2014 Supplementary Table 4's **Correction Notice** (extracted lines 1–6)
states that its original online version omitted the first five rows; the corrected
table was posted 20 May 2015. The corrected rows 9–13 are `P5_i5_1` through
`P5_i5_5`, with respective 8-base index values `TATAGCCT`, `ATAGAGGC`,
`CCTATCCT`, `GGCTCTGA`, `AGGCGAAG`. The original values were absent, not
conflicting bases. 🟡 The sci-CAR model uses generic placeholders for its first
well Tn5 index because the Cao 2018 supplementary oligo table remains gated;
it does not assign these five Amini members to sci-CAR.
