# Methyl-HiC — source notes

Defining paper: Li et al., *Nature Methods* (2019), [doi:10.1038/s41592-019-0502-z](https://doi.org/10.1038/s41592-019-0502-z).

## Sources read

- 🟢 PMC full text, main description and Methods.
- 🟢 Original Supplementary Protocol, including the reagent list, bulk on-bead
  adapter-ligation procedure and complete scMethyl-HiC oligo table.

## Bulk Methyl-HiC

- 🟢 Fixed nuclei are digested overnight with DpnII. Ends are filled with biotinylated nucleotides and ligated in situ.
- 🟢 After crosslink reversal the DNA is sonicated to about 400 bp and pulled down with streptavidin beads.
- 🟢 Adapters are ligated on beads, then the molecules undergo EZ DNA Methylation-Gold bisulfite conversion and KAPA HiFi HotStart Uracil+ amplification.
- 🟢 Libraries were sequenced paired-end on HiSeq 2500 or 4000 with balanced libraries or PhiX because of bisulfite base imbalance.
- 🟢 The bulk supplementary protocol explicitly lists `TruSeq PCR Primers`, ligates an
  indexed `Illumina Adapter` on beads, and then performs uracil-tolerant PCR. Standard
  TruSeq Read 1, Index 1 and Read 2 entry sites can therefore be located on the compatible
  final shell.
- 🟡 The exact ligated adapter strands and selected index member are not printed. The page
  marks the canonical TruSeq-compatible bases as inferred rather than pretending they were
  transcribed from the paper; the ligation-ready end remains unresolved.

## Related single-cell branch

- 🟢 scMethyl-HiC sorts nuclei after proximity ligation, bisulfite-converts each nucleus,
  uses indexed P5 random priming, Exonuclease I/SAP cleanup, and Adaptase addition of P7
  before dual-index PCR. The supplement prints all four P5-tailed random-primer families,
  eight custom P5 indexing primers and twelve custom P7 indexing primers.
- 🟡 This page models the defining bulk Methyl-HiC library. The single-cell branch is recorded here rather than conflated into the bulk molecule.
