# Structure-seq2

**Defining source:** Ritchey *et al.*, *Nucleic Acids Research* (2017),
doi:10.1093/nar/gkx533. The PMC text, main figure and nucleotide-explicit supplementary
Figure S1 were read.

## Evidence used

- 🟢 The paper specifies in-vivo DMS probing, random-hexamer RT and matched −DMS controls.
- 🟢 Methods publish the phosphorylated/blocked six-degenerate-base hairpin donor, T4 DNA
  ligation, both PCR primers, 220–600-nt selection and 150-cycle MiSeq endpoint.
- 🟢 Supplementary Figure S1 gives the custom primer sequence used to skip the common
  22-nt prefix; it is the terminal 35 nt of the published forward primer.
- 🟢 PAGE and biotin-dCTP/streptavidin branches are alternative purification routes.

## Construct logic

Random priming samples transcript positions, DMS establishes the stop, and the degenerate
hairpin creates a ligatable duplex at diverse cDNA ends. The custom primer begins at the
unique boundary rather than consuming cycles on invariant adapter sequence.
