# SUM-seq — research notes

## Sources read

- 🟢 Llorens-Rico et al. 2025, DOI `10.1038/s41592-025-02700-8`, PMC article, Methods and Supplementary Table 1.

## Construct path

🟢 Fixed nuclei receive sample identity twice by modality-specific routes: barcoded Tn5 for ATAC and barcoded oligo(dT) reverse transcription for RNA. RNA/cDNA hybrids are tagmented to add the primer site required for microfluidic indexing. 🟢 Sample-indexed nuclei are pooled and overloaded into 10x droplets, where a second droplet barcode is added. Cell identity is the sample-index × droplet-barcode combination; RNA and ATAC products are split only after droplet breakage and pre-amplification.
