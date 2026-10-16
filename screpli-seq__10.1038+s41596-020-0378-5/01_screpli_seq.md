# scRepli-seq — research notes

## Citation

Miura H et al. *Mapping replication timing domains genome wide in single mammalian cells with single-cell DNA replication sequencing.* Nature Protocols 15, 4058–4100 (2020). [doi:10.1038/s41596-020-0378-5](https://doi.org/10.1038/s41596-020-0378-5)

## Sources read

- 🟢 The defining protocol record, abstract, figures and reagent/library description.
- 🟢 The open method review by Miura and colleagues summarizing the defining workflow.
- 🟢 The independent updated scRepli-seq protocol in [PMC8812816](https://pmc.ncbi.nlm.nih.gov/articles/PMC8812816/), used only to corroborate adapter-library geometry and not to replace the defining SeqPlex workflow.
- 🟢 The publisher's Supplementary Software 1/2 (MOESM3/4) gives the read-trimming patterns for the original SeqPlex implementation. The source bundle and oligo scraper were checked against the original supplementary software.

## 1. What it is

🟢 scRepli-seq does not label nascent DNA. It infers replication from the twofold copy-number difference between genomic domains already replicated and domains not yet replicated in one S-phase cell.

## 2. Construct-changing steps

1. 🟢 Sort individual early-, mid- or late-S cells by DNA content; mid-S cells most directly expose both copy-number states.
2. 🟢 Whole-genome amplify each cell using the SeqPlex Enhanced DNA Amplification Kit, a degenerate oligonucleotide-primed PCR method.
3. 🟢 Remove/trim SeqPlex amplification sequence and prepare an individually indexed Illumina library from each amplified cell. The defining protocol permits KAPA LTP or NEBNext Ultra II preparation.
4. 🟢 Process G1 cells in parallel as controls for amplification and mappability bias.

## 3. Final library and reads

🟢 The final product is an ordinary indexed paired-end Illumina genomic library. Its sample index assigns it to one cell; there is no inline cell barcode or molecule UMI in the defining protocol.

🟢 Mapped coverage, rather than a label-derived mismatch, distinguishes replicated from unreplicated domains.

🟢 Supplementary Software 1 gives an example Illumina single-index trimming pattern `GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG`; the six-N window is the index. It says to infer the SEQXE adapter from the reads and prints an example predicted pattern `TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN` (software lines 383–424). Supplementary Software 2 repeats both patterns. These are processing patterns, not an orderable SeqPlex kit oligo sequence.

🟡 The page models standard Illumina-compatible primer sites and keeps the proprietary SeqPlex priming region as a named role. The example SEQXE pattern is a source-reported prediction from FASTQ reads; its exact relationship to the synthesized degenerate primer and molecular end remains unresolved.

## 4. Boundary

🟢 The later DOP-PCR/barcoded pooling protocol is useful context but is molecularly distinct and is not silently substituted here.
