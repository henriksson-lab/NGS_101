# Seq-Scope — research notes

Defining source: Cho *et al.*, [10.1016/j.cell.2021.05.010](https://doi.org/10.1016/j.cell.2021.05.010).

## Evidence used

- 🟢 Seq-Scope uses solid-phase amplification of randomly barcoded HDMI oligos on an Illumina flow cell.
- 🟢 A first sequencing/imaging run maps every HDMI cluster to a physical coordinate.
- 🟢 Tissue RNA is captured through oligo(dT), reverse-transcribed on the decoded surface, and recovered with spatial barcode and UMI for transcript sequencing.
- 🟢 The open defining paper prints the complete HDMI-DraI and HDMI32-DraI seed oligos, including P5/P7, TruSeq Read 1, constrained-random HDMI, HR1, oligo(dT), DraI site and DraI-adapter regions.
- 🟢 It prints the full `Read1-DraI` custom sequencing primer used for the coordinate-decoding run and the TruSeq Read2-conjugated 9-nt random primer used for second-strand synthesis; the latter random tract serves as the UMI.
- 🟢 Both amplification rounds are source-exact, including their 3′ phosphorothioate bonds and the 8-mer sample-index window in the final reverse primer. The completed library has the standard TruSeq sequencing-primer sites.
- 🟢 The paper also prints the experimental HDMI-TruEcoRI seed, Read1-EcoRI custom sequencing primer, and the 30-dA/8-nt-UMI oligo used to test a UMI-encoded array. These are retained as a distinct alternative rather than blended into the main DraI workflow.

## Protocol boundary

The defining two-pass use of the Illumina flow cell is retained rather than reducing Seq-Scope to a generic spatial oligo(dT) array.

All paper-published oligos are now listed on the schematic page. The main final-library view follows the DraI implementation used for the liver and colon data; the EcoRI set is explicitly an experimental alternative.
