# YAMAT-seq

**Defining source:** Shigematsu *et al.*, *Nucleic Acids Research* (2017),
doi:10.1093/nar/gkx005. The PMC full text and supplement were read.

## Evidence used

- 🟢 Methods lines 52–60 publish all five Y-adapter oligos, their DNA/RNA composition,
  end modifications, Rnl2 annealing/ligation, TruSeq Small RNA RT-PCR and 100-nt
  single-read sequencing.
- 🟢 Lines 80–88 establish deacylation, recognition of the mature-tRNA 5′-NCCA-3′
  overhang and sealing of both acceptor-stem nicks.
- 🟡 The final flow-cell ends are shown by molecular role because the paper names the
  Illumina kit but does not reproduce its index or sequencing-primer oligos.

## Construct logic

The four Y-5′ variants differ only in their terminal ribodinucleotide and collectively
pair with all discriminator-base/NCCA ends. Y-3′-AD is 5′ phosphorylated and 3′ blocked.
Annealing creates two Rnl2 substrates; amplification is possible only after both relevant
adapter–tRNA boundaries have been sealed.
