# ReDeeM — research notes

Defining source: Weng *et al.*, [10.1038/s41586-024-07066-z](https://doi.org/10.1038/s41586-024-07066-z).

## Evidence used

- 🟢 ReDeeM modifies droplet-based 10x Multiome using whole cells and conditions that retain mitochondrial DNA while preserving RNA and ATAC libraries.
- 🟢 Tiled mtDNA probes hybrid-capture a third, deeply sequenced library carrying the same cell barcodes.
- 🟢 Cell barcode plus mtDNA fragment endpoints act as endogenous molecular identifiers for consensus error correction.

## Protocol boundary

ReDeeM captures mtDNA from Multiome material; MAESTER instead enriches mitochondrial transcripts from scRNA-seq cDNA.
