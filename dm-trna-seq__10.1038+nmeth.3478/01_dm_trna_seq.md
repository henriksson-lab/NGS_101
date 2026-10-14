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

## Boundary

This is the 2015 demethylase/TGIRT protocol. It is not mim-tRNAseq, which preserves
modification-induced mismatches, and not YAMAT-seq, which ligates a Y adapter to NCCA.
