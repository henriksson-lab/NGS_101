# Methyl-HiC — source notes

Defining paper: Li et al., *Nature Methods* (2019), [doi:10.1038/s41592-019-0502-z](https://doi.org/10.1038/s41592-019-0502-z).

## Sources read

- 🟢 PMC full text, main description and Methods.

## Bulk Methyl-HiC

- 🟢 Fixed nuclei are digested overnight with DpnII. Ends are filled with biotinylated nucleotides and ligated in situ.
- 🟢 After crosslink reversal the DNA is sonicated to about 400 bp and pulled down with streptavidin beads.
- 🟢 Adapters are ligated on beads, then the molecules undergo EZ DNA Methylation-Gold bisulfite conversion and KAPA HiFi HotStart Uracil+ amplification.
- 🟢 Libraries were sequenced paired-end on HiSeq 2500 or 4000 with balanced libraries or PhiX because of bisulfite base imbalance.
- 🔴 The exact bulk adapter oligonucleotide sequences are not printed. The public page therefore leaves the arms unresolved and does not assert a sequencing-primer landing site.

## Related single-cell branch

- 🟢 scMethyl-HiC sorts nuclei after proximity ligation, bisulfite-converts each nucleus, uses indexed P5 random priming, Exonuclease I/SAP cleanup, and Adaptase addition of P7 before dual-index PCR.
- 🟡 This page models the defining bulk Methyl-HiC library. The single-cell branch is recorded here rather than conflated into the bulk molecule.
