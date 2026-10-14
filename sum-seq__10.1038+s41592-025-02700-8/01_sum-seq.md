# SUM-seq — research notes

## Sources read

- 🟢 Llorens-Rico et al. 2025, DOI `10.1038/s41592-025-02700-8`, PMC article, Methods and Supplementary Table 1.

## Construct path

🟢 Fixed nuclei receive sample identity twice by modality-specific routes: barcoded Tn5 for ATAC and barcoded oligo(dT) reverse transcription for RNA. RNA/cDNA hybrids are tagmented to add the primer site required for microfluidic indexing. 🟢 Sample-indexed nuclei are pooled and overloaded into 10x droplets, where a second droplet barcode is added. Cell identity is the sample-index × droplet-barcode combination; RNA and ATAC products are split only after droplet breakage and pre-amplification.

## Structural oligos in the original supplement

🟢 Supplementary Table 1b–d, in the publisher XLSX (`MOESM3`), prints the general reagents and the 96-member ATAC and RNA sample-index families. Its indices are **11 bases**, not the eight bases formerly drawn in the page. `SUMseq-ATAC-1` is `CATACGGCAGAAGACGAACACCTCCAACTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` (index `ACCTCCAACTA`). `SUMseq-RNA-1` is `CAGACGTGTGCTCTTCCGATCTNNNNNNNNNNAGTGATTAGCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN` (10-base UMI and index `AGTGATTAGCA`). These are examples from separate published index families, not an identical shared oligo.

🟢 The common Tn5 loading oligo is `[Phos]CTGTCTCTTATACACATCT[ddc]`. The hopping blocker is `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-[invdT]`. The table also gives `Nextera Read 1` `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, phosphorothioate-ending Partial-P5 `AATGATACGGCGACCACCGAG*A`, and Partial-P7 `CAAGCAGAAGACGGCATACGA*G`. Their chemical modifications are part of the ordering sequences.

🟡 The current graph still compresses these published oligos into generic branch objects and a representative final shell. Its full RNA/ATAC junctions and PCR binding geometry need source-based reconstruction. 🔴 The proprietary 10x droplet-bead oligo is not printed in Supplementary Table 1.
