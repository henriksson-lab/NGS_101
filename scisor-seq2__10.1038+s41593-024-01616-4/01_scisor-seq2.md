# ScISOr-seq2 — research notes

## Sources read

- 🟢 Joglekar et al. 2024, DOI `10.1038/s41593-024-01616-4`, PMC JATS Methods and supplements.

## Construct path

🟢 ScISOr-seq2 starts from 10x Chromium 3′ cDNA, so the cell barcode and UMI already belong to the original RT product. Twelve Partial-Read1-only cycles preferentially amplify the barcode-bearing strand; six exponential cycles add the Partial TSO end. Exome probes enrich spliced cDNA (LAP-CAP). 🟢 The enriched pool is compatible with both PacBio SMRTbell and ONT ligation preparation. The two exact partial primers shown on the page are stated in Methods.
