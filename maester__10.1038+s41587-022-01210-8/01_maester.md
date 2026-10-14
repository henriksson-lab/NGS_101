# MAESTER — research notes

Defining source: Miller *et al.*, [10.1038/s41587-022-01210-8](https://doi.org/10.1038/s41587-022-01210-8).

## Evidence used

- 🟢 MAESTER enriches mitochondrial transcripts from existing 3′ single-cell RNA-seq cDNA and increases mitochondrial coverage while retaining the parent expression library.
- 🟢 The supplementary oligo tables provide separate implementations for Seq-Well and 10x Genomics 3′ v3 cDNA.
- 🟢 The 10x workbook (Supplementary Table 3) prints the complete 10x-specific mitochondrial target-primer inventory, common `CACCCGAGAATTCCA` handle, 24 i7 and 24 i5 sample-index primer sequences, optional `PartialRead1`/`PartialTSO`, and the two custom index-sequencing spike-in primers. The Seq-Well workbook (Supplementary Table 2) separately prints 5′-biotinylated mitochondrial primers, `SMART-AC`, `CR1P`, Seq-Well index sequencing primer and its indexed PCR primers. These full tables are the authoritative inventories; the two implementations are not interchangeable.
- 🟢 In the 10x implementation, a representative mitochondrial primer is `PvG1218_MT-ND1_4009` = `CACCCGAGAATTCCAAACACCCTCACCACTACAATCT`. In the Seq-Well implementation, the corresponding `TyM61-MT-ND1_4009` is printed with `/5Biosg/` and a Nextera-bearing handle. The 10x branch omits the Seq-Well biotin enrichment step.
- 🟢 Multiplex target amplification preserves cell barcode and UMI identity, allowing mitochondrial variants to be assigned to expression-profiled cells.

## Protocol boundary

The schematic shows the 10x-compatible enrichment branch; the original parent 10x library chemistry remains on its own page.
