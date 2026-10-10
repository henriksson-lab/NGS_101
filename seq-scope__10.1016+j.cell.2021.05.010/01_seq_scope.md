# Seq-Scope — research notes

Defining source: Cho *et al.*, [10.1016/j.cell.2021.05.010](https://doi.org/10.1016/j.cell.2021.05.010).

## Evidence used

- 🟢 Seq-Scope uses solid-phase amplification of randomly barcoded HDMI oligos on an Illumina flow cell.
- 🟢 A first sequencing/imaging run maps every HDMI cluster to a physical coordinate.
- 🟢 Tissue RNA is captured through oligo(dT), reverse-transcribed on the decoded surface, and recovered with spatial barcode and UMI for transcript sequencing.

## Protocol boundary

The defining two-pass use of the Illumina flow cell is retained rather than reducing Seq-Scope to a generic spatial oligo(dT) array.
