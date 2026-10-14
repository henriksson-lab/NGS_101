# TIRTL-seq

**Defining paper:** Pogorelyy et al., *Nature Methods* (2026), doi:10.1038/s41592-025-02907-9. The open article and its methods were read on the publisher site.

## Molecular boundary

🟢 Cells are distributed over replicate wells and lysed during miniaturized RT. A first multiplex PCR uses V- and C-segment primers carrying partial Nextera handles; C-segment reverse primers also carry a plate barcode.

🟢 A second PCR adds full Illumina adapters and one well-specific UDI pair. Alpha/beta pairing is inferred from clonotypes' shared occurrence patterns across replicate wells, not from two chains joined in one molecule.

## Evidence boundary

🟢 The published Supplementary Data 1 (`41592_2025_2907_MOESM3_ESM.xls.txt`) gives the full ordered panel in four sheets: 88 V-segment PCR-I forward primers, 24 plate-barcoded C-segment PCR-I reverse primers (12 alpha and 12 beta), 384 well-specific i5 PCR-II forward primers and 384 matching i7 PCR-II reverse primers. The V primers share the printed `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` partial Nextera arm; the C primers begin `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` before a 6- or 8-base plate identifier and the C-specific arm. The i5 and i7 sheets print complete P5/P7-tailed UDI primer sequences. There are no modification annotations in these four sheets.

🟢 The main article methods (`s41592-025-02907-9_PMC12791010.xml.txt`, line 218) also specify commercial random hexamer and oligo(dT) RT primers. Their synthesized sequences or oligo(dT) length are not printed and are not members of the 880-oligo PCR panel.

🟡 The schematic currently depicts the primer, plate barcode and well UDI **roles**; it does not encode all 880 panel members as segments or derive primer-binding and read-cycle spans from a selected panel row. The panel is exact for the published implementation, while target choice and index combinations can vary. No bases were inferred from raw reads because Supplementary Data 1 supplies the sequences directly.
