# CoBATCH — research notes

Defining source: Wang et al., 2019, [doi:10.1016/j.molcel.2019.07.015](https://doi.org/10.1016/j.molcel.2019.07.015).

## Sources read

- 🟢 Author preprint 10.1101/590661 Methods: antibody targeting, PAT binding, barcoded tagmentation, pooling, redistribution and indexed PCR.
- 🟢 Preprint Supplementary Table 2 workbook (`590661_media-3.xlsx`): Mosaic-A/B, phosphorylated common ME oligo, PAT-T5 and PAT-T7 barcode series, Nextera i5/i7 PCR primers.
- 🔴 The defining *Molecular Cell* article and its final supplement were not retrieved by the source fetcher. Exact oligo transcription below is verified against the author preprint; final-version differences remain unassessed.

## Molecular path

🟢 Antibody-bound Protein-A–Tn5 deposits barcoded adapters at target chromatin in a first plate. The cells are pooled and redistributed; i5/i7 PCR indexes encode the second plate. A cell is identified by the combined T5, T7, i5 and i7 values.

🟢 The first-plate oligos are complete transposon top strands, not bare barcode windows. T5 has S5 + `TCCACGC` + 8-nt T5 barcode + `GCGATCGAGGACGGC` + ME. T7 has S7 + `CTGTCCCTGTCC` + 8-nt T7 barcode + `CACCGTCTCCGCCTC` + ME. Both anneal only over their terminal ME to the 5′-phosphorylated common bottom `CTGTCTCTTATACACATCT` before pA–Tn5 loading.

🟢 The source workbook contains 24 T5 and 25 T7 barcode oligos. The schematic draws the exact T5-1 and T7-2 representatives and attaches whitelist feature metadata; it does not render repeated placeholder bases as if they were a published barcode.

🟢 After targeted transfer and gap fill, the first-plate codes sit between their outer and inner connectors at both fragment ends. Indexed PCR adds P5–i5–S5 and S7–i7–P7, making the four-part cellular identity explicit in the completed duplex.
