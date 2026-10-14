# ScISOr-seq2 — research notes

## Sources read

- 🟢 Joglekar et al. 2024, DOI `10.1038/s41593-024-01616-4`, PMC JATS Methods and supplements.

## Construct path

🟢 ScISOr-seq2 starts from 10x Chromium 3′ cDNA, so the cell barcode and UMI already belong to the original RT product. Twelve Partial-Read1-only cycles preferentially amplify the barcode-bearing strand; six exponential cycles add the Partial TSO end. Exome probes enrich spliced cDNA (LAP-CAP). 🟢 The enriched pool is compatible with both PacBio SMRTbell and ONT ligation preparation. The two exact partial primers shown on the page are stated in Methods, and their complete terminal sites determine the LAP PCR geometry without requiring a transcript-specific coordinate.

🟢 The same printed Partial Read1 `CTACACGACGCTCTTCCGATCT` and Partial TSO `AAGCAGTGGTATCAACGCAGAGTACAT` also serve as blocking oligos during exome-probe hybridization, then as primers again for post-capture amplification. The source gives both roles; these are not additional undisclosed oligos.

🔴 The source names PacBio SMRTbell Express Template Prep Kit 2.0 and ONT SQK-LSK110 for the terminal platform libraries but does not publish their adapter or sequencing-primer bases. Those terminal commercial regions remain explicit placeholders.
