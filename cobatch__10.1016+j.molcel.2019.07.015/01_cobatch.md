# CoBATCH — research notes

Defining source: Wang et al., 2019, [doi:10.1016/j.molcel.2019.07.015](https://doi.org/10.1016/j.molcel.2019.07.015).

## Sources read

- 🟢 Author preprint 10.1101/590661 Methods: antibody targeting, PAT binding, barcoded tagmentation, pooling, redistribution and indexed PCR.
- 🟢 Preprint Supplementary Table 2 workbook (`590661_media-3.xlsx`): Mosaic-A/B, phosphorylated common ME oligo, PAT-T5 and PAT-T7 barcode series, Nextera i5/i7 PCR primers.
- 🟢 Recovered original *Molecular Cell* Table S1 (`mmc2.xlsx`), sheets Mosaic-A&B / PA-Tn5_barcodes / Nextera_index: confirms the adapter and PCR families previously read in preprint Table 2.
- 🟢 Original supplemental PDF (`mmc1.pdf`) contains Figures S1–S4 and captions. Original Table S2 (`mmc3.xlsx`) records antibodies, library-preparation labels and sequencing statistics; Table S3 (`mmc4.xlsx`) compares CoBATCH/scChIC-seq performance, and Table S4 (`mmc5.xlsx`) gives tissue-digestion recipes. These are not additional sequencing-primer tables.
- 🔴 The final main article remains inaccessible. Its STAR Methods and all relevant predecessor-source implications have not yet been fully checked twice; no missing-information conclusion is accepted.

## Molecular path

🟢 Antibody-bound Protein-A–Tn5 deposits barcoded adapters at target chromatin in a first plate. The cells are pooled and redistributed; i5/i7 PCR indexes encode the second plate. A cell is identified by the combined T5, T7, i5 and i7 values.

🟢 The first-plate oligos are complete transposon top strands, not bare barcode windows. T5 has S5 + `TCCACGC` + 8-nt T5 barcode + `GCGATCGAGGACGGC` + ME. T7 has S7 + `CTGTCCCTGTCC` + 8-nt T7 barcode + `CACCGTCTCCGCCTC` + ME. Both anneal only over their terminal ME to the 5′-phosphorylated common bottom `CTGTCTCTTATACACATCT` before pA–Tn5 loading.

🟢 The source workbook contains 24 T5 and 25 T7 barcode oligos. The schematic draws the exact T5-1 and T7-2 representatives and attaches whitelist feature metadata; it does not render repeated placeholder bases as if they were a published barcode.

🟢 After targeted transfer and gap fill, the first-plate codes sit between their outer and inner connectors at both fragment ends. Indexed PCR adds P5–i5–S5 and S7–i7–P7, making the four-part cellular identity explicit in the completed duplex.

## Sequencing evidence correction

🟢 Original published Table S1 prints the T5/T7 transfer oligos and i5/i7 PCR
primers. Its original sheet rows agree with the representative complete oligos
assembled on this page. The author preprint Methods, lines423–452, specify
paired-end 150-base reads for the Nextera library strategy.

🟡 The former page declared inner-connector-plus-ME sequences as physically used
custom genomic read primers, and showed 50-base reads. Those were unsupported
repository assignments. The displayed adapter subsequences alone cannot
establish the loaded sequencing primers. They have been removed; read length
is qualified to the author preprint, and exact sequencing-entry/index-read
geometry stays unresolved pending the final main article and linked sources.
This correction is not attributed to a publication error.
