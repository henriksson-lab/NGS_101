# EpiDamID — source notes

Defining paper: Rang and de Luca et al., *Molecular Cell* (2022), [doi:10.1016/j.molcel.2022.03.009](https://doi.org/10.1016/j.molcel.2022.03.009).

- 🟢 PMC full text, STAR Methods and key-resources oligonucleotide table were read.
- 🟢 Chromatin reader domains, proteins or mintbodies recruit Dam; Dam methylates adenine in nearby GATC motifs in living cells.
- 🟢 In scDam&T-seq, DpnI digests methylated GATC after cDNA second-strand synthesis, DamID2 adapters are ligated, wells are pooled and T7 IVT runs for 14 h.
- 🟢 The amplified RNA is fragmented and completed with the CEL-Seq2 workflow. Libraries were sequenced on NextSeq 500 (75 bp) or NextSeq 2000 (100 bp).
- 🟢 AdRt, AdRb and AdR_PCR are transcribed from the paper's oligo table.
- 🟢 The same key-resources table also prints RandomhexRT (`GCCTTGGCACCCGAGAATTCCANNNNNN`), RNA PCR primer 1 (`AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A`) and an example indexed RNA PCR primer (`CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCC*A`). The `*` marks the terminal phosphorothioate linkage in each PCR primer. These three reagents were missing from the earlier oligo list.
- 🔴 The paper does not print a single complete final strand containing all CEL-Seq2-derived sequencing sites, so no base-level final primer map is claimed. The RNA primer sequences establish their own terminal arms but do not settle the full path between the DamID adapter and those arms.
