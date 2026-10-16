# DM-tRNA-seq / DM-TGIRT-seq

**Defining source:** Zheng *et al.*, *Nature Methods* (2015),
doi:10.1038/nmeth.3478. The PMC author manuscript and supplements were read.

## Evidence used

- 🟢 Full text lines 52–54 establish combined wild-type AlkB and D135S treatment,
  followed by TGIRT template switching rather than RNA-adapter ligation.
- 🟢 Methods lines 80–96 give deacylation, the blocked 41-nt RNA starter, the one-base-T
  DNA-primer overhang, 60 °C TGIRT, CircLigase II, both PCR primers and HiSeq sequencing.
- 🟡 The separate sequencing primer is not printed, so its exact annealing span is not
  asserted.
- 🟢 The source prints the RNA starter as
  `AGAUCGGAAGAGCACACGUCUAGUUCUACAGUCCGACGAUC/3SpC3/`, and prints both complete
  post-circularization PCR primers. The 3′ spacer blocks starter extension.
- 🟡 The complementary DNA starter primer is described as radiolabeled with a
  single-T 3′ overhang, but its full ordered sequence and label position are not
  printed. Complementarity permits a base-level reconstruction only if the exact
  annealing register is established; no synthetic modification can be inferred from it.
- 🟢 A later [InGex TGIRT-III small-RNA/CircLigase protocol](https://fnkprddata.blob.core.windows.net/domestic/data/datasheet/IGX/TGIRT10.pdf) (updated 2017-01-23,
  pp. 1–3) prints this 41-nt starter RNA, its complementary DNA primer and the same
  PCR-primer architecture. It uses a degenerate `N` DNA-primer overhang, whereas
  the defining 2015 paper explicitly uses `T`. This is a later variant, so the
  diagram follows the 2015 `T` statement. The vendor protocol identifies the
  Read 1 binding region as the NEBNext small-RNA primer site; it does not identify
  the exact HiSeq reagent/primer used for the published run.
- 🔴 The 2015 run's HiSeq model, flow-cell chemistry and sequencing kit revision
  are unspecified. Illumina's HiSeq 2500 primers page lists several possible
  primer mixes, so the distinct run-primer oligo and its exact landing site remain
  unresolved rather than confirmed absent from all linked materials.

## Boundary

This is the 2015 demethylase/TGIRT protocol. It is not mim-tRNAseq, which preserves
modification-induced mismatches, and not YAMAT-seq, which ligates a Y adapter to NCCA.
