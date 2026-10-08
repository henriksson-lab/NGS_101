# CITE-seq — research notes

## 1. What it is

CITE-seq converts surface-protein abundance into sequenceable antibody-derived tags while measuring the transcriptome from the same single cell.

## 2. Sources read

- 🟢 Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, Nature Methods 14, 865–868 (2017), DOI [10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380), open manuscript [PMC5669064](https://pmc.ncbi.nlm.nih.gov/articles/PMC5669064/).
- 🟢 Main Methods and supplementary protocol establish the poly(A)-tailed antibody tag, DNA-dependent extension during RT, size separation and independent ADT library PCR.

## 3. Construct-changing steps

1. Conjugate antibody identity oligos containing PCR handle, antibody barcode and poly(A) tail.
2. Label cells and wash away free tags.
3. Co-encapsulate cells and barcoded beads; oligo-dT captures both mRNA and ADT poly(A).
4. Reverse transcription/extension gives both molecule classes the cell barcode and UMI.
5. Size-separate large cDNA from short ADTs and make separate libraries.

## 4. Evidence boundary

🟢 The defining paper prints representative antibody oligos and establishes the structure. 🟡 The final page uses role placeholders for panel-specific antibody barcodes and normalizes platform arms for computed primer placement.
