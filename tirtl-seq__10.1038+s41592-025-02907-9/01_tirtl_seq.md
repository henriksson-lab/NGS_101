# TIRTL-seq

**Defining paper:** Pogorelyy et al., *Nature Methods* (2026), doi:10.1038/s41592-025-02907-9. The open article and its methods were read on the publisher site.

## Molecular boundary

🟢 Cells are distributed over replicate wells and lysed during miniaturized RT. A first multiplex PCR uses V- and C-segment primers carrying partial Nextera handles; C-segment reverse primers also carry a plate barcode.

🟢 A second PCR adds full Illumina adapters and one well-specific UDI pair. Alpha/beta pairing is inferred from clonotypes' shared occurrence patterns across replicate wells, not from two chains joined in one molecule.

## Evidence boundary

🟢 The paper publishes the oligo panel in Supplementary Data 1. Panel sequences are target- and implementation-specific, so the schematic represents their conserved molecular roles rather than treating one panel as universal.
