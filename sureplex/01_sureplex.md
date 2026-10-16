# SurePlex with VeriSeq PGS

## Defining sources

🟢 Illumina, *SurePlex Summary Protocol Reference Guide* 15053626 v02, and *VeriSeq PGS Library Prep Reference Guide* 15052877 v05.

## Molecular workflow

🟢 SurePlex accepts a single cell in 2.5 µl PBS (guide p. 6), performs extraction (p. 8), a 12-cycle preamplification with Pre-amp Buffer and Enzyme (p. 9), then a separate 14-cycle amplification with Amplification Buffer and Enzyme (p. 10). The inspected guide does not print the preamplification oligo bases or formulation; whether another relevant vendor document supplies them remains open.

🟢 The associated VeriSeq PGS guide v05, pp. 9–11, uses 1 ng SurePlex dsDNA (5 µl at 0.2 ng/µl), an engineered transposome, NT neutralization, and 12-cycle PCR with S503/S504 and N701–N712 primers. The guide supports both single-index 12-plex and dual-index 24-plex libraries; the schematic shows the dual-index branch. The current official Illumina Adapter Sequences table v22 identifies the S503/S504 and N7xx index bases, but study-specific member selection is unspecified.

🟢 The intended MiSeq acquisition is single-read: Read 1 is analysed, index reads demultiplex samples, and Read 2 is not acquired.

🟢 The [official VeriSeq PGS Guidance Note v00](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/veriseq-pgs/veriseq-pgs-guidance-note-1000000138837-00.pdf), p. 2, reports possible PCR artifacts in parallel no-template controls; these may produce unmapped reads that the system filters. It does not identify the molecules or their origin. The graph renders this as discarded material with a warning, separate from the sample library. The note also calls the assay a 36-cycle short-read recipe, but does not specify the index/read cycle split there.

🟢 The [BlueFuse Workflow Manager quick card v01](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bluefuse-workflow-manager/bwm-qrc-veriseq-pgs-15056206-01.pdf), p. 1, has the operator choose i7 and i5 entries for each sample. The schematic shows the permitted family, not one chosen well's indexes.

🔴 SurePlex and PicoPLEX share a commercial lineage, but they are kept as separate named protocols because PicoPLEX DNA-seq integrates dual indexing whereas SurePlex hands its WGA product to the separate VeriSeq PGS library kit.
