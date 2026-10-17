# EpiDamID — source notes

Defining paper: Rang and de Luca et al., *Molecular Cell* (2022), [doi:10.1016/j.molcel.2022.03.009](https://doi.org/10.1016/j.molcel.2022.03.009).

- 🟢 PMC full text, STAR Methods and key-resources oligonucleotide table were read.
- 🟢 Chromatin reader domains, proteins or mintbodies recruit Dam; Dam methylates adenine in nearby GATC motifs in living cells.
- 🟢 In scDam&T-seq, DpnI digests methylated GATC after cDNA second-strand synthesis, DamID2 adapters are ligated, wells are pooled and T7 IVT runs for 14 h.
- 🟢 The amplified RNA is fragmented and completed with the CEL-Seq2 workflow. Libraries were sequenced on NextSeq 500 (75 bp) or NextSeq 2000 (100 bp).
- 🟢 AdRt, AdRb and AdR_PCR are transcribed from the paper's oligo table (PMC9153956 JATS text lines 574–583). The table attributes them to Vogel et al. 2007 and labels AdR_PCR “for m6A-PCR”; the Methods use m6A-PCR with gel electrophoresis to check Dam methylation during clone screening (lines 821, 837). They are **not** the scDam&T-seq DamID2 adapter. The earlier model drew AdRt/AdRb as the DamID2 adapter; that was a repository error, corrected 2026-10-17.
- 🟢 The same key-resources table also prints RandomhexRT (`GCCTTGGCACCCGAGAATTCCANNNNNN`), RNA PCR primer 1 (`AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A`) and an example indexed RNA PCR primer (`CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCC*A`). The `*` marks the terminal phosphorothioate linkage in each PCR primer. These three reagents were missing from the earlier oligo list.
- 🟢 Key Resources list “Adapters for DamID2, top and bottom oligonucleotides” and “Primers for SORT-seq / CEL-Seq2” only as Markodimitraki et al. 2020 (lines 627–629), and the Methods state that scDam&T-seq “was performed as previously described in detail (Markodimitraki et al., 2020)” with halved volumes (line 853). That establishes inheritance of the DamID2 adapter design from the cited protocol.
- 🟢 Cited protocol (Markodimitraki 2020, PMC7779467 JATS text lines 135, 284–288): example DamID adapter top `GGTGATCCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATCNNNTGCANNNTATGGA`, bottom `/5Phos/TCCATANNNTGCANNNGATCGTCGGACTGTAGAACTCTGAACCCCTATAGTGAGTCGTATTACCGGGAGCTT`: 6-nt fork, T7 promoter, RA5/P5-side handle, split 8-nt barcode alternating with a split 6-nt UMI, blunt GA end. The model draws this example; the barcode set used on each EpiDamID plate is not assigned.
- 🟢 Single-cell steps (line 853): sort into 50 nL barcoded CEL-Seq2 primer; lysis 65 °C 5 min; SuperScript II RT 42 °C 2 h; second strand (E. coli DNA ligase, RNase H, DNA polymerase I) 16 °C 2 h; proteinase K 50 °C 10 h; DpnI 37 °C 6 h; DamID2 adapters to 25 nM and T4 ligase 16 °C 16 h; pool, 0.8-volume bead cleanup; MEGAscript T7 IVT 14 h; aRNA fragmentation 94 °C 90 s; CEL-Seq2 library preparation with 8–10 PCR cycles.
- 🟡 The DamID-branch final strand is derived: RNA PCR primer 1's 3′ 21 nt equal the adapter's RA5 region, and the RNA PCR index primer ends in the 3′ 21 nt of the RandomhexRT tail (the tail's 5′ G differs). The 6-nt index boundary within the example index primer (`CGTGAT`) is inferred from the “Follow Illumina design” attribution; the source does not delimit it.
- 🔴 The paper does not print a single complete final strand containing all CEL-Seq2-derived sequencing sites, so no base-level final primer map is claimed. The RNA primer sequences establish their own terminal arms but do not settle the full path between the DamID adapter and those arms.

## Open questions (2026-10-17 bounded pass)

- Sequencing primers and paired-end configuration for the EpiDamID runs (the cited 2020 protocol uses NextSeq 500 2 × 75; EpiDamID states only 75-bp or 100-bp reads).
- DamID2 and CEL-Seq2 barcode sets per plate (Markodimitraki 2020 Supplementary Tables 1–2 not read in this pass).
- The CEL-Seq2 transcript branch is carried through the graph but not drawn base by base.
- Rooijers 2019 four-base-fork adapter variant and Hashimshony 2016 CEL-Seq2 library-prep details (RT enzyme, PCR polymerase) not re-read.
