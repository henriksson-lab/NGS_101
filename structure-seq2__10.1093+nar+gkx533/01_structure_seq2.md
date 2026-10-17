# Structure-seq2

**Defining source:** Ritchey *et al.*, *Nucleic Acids Research* (2017),
doi:10.1093/nar/gkx533. The PMC text, main figure and nucleotide-explicit supplementary
Figure S1 were read.

## Evidence used

- 🟢 The paper specifies in-vivo DMS probing, random-hexamer RT and matched −DMS controls.
- 🟢 Methods publish the phosphorylated/blocked six-degenerate-base hairpin donor, T4 DNA
  ligation, both PCR primers, 220–600-nt selection and 150-cycle MiSeq endpoint.
- 🟢 Nucleotide-explicit Supplementary Figure S1 gives the 27-nt RT primer and the
  custom primer used to skip the common 22-nt prefix; the custom primer is the terminal
  33 nt of the published forward primer.
- 🟢 PAGE and biotin-dCTP/streptavidin branches are alternative purification routes.

## Construct logic

Random priming samples transcript positions, DMS establishes the stop, and the degenerate
hairpin creates a ligatable duplex at diverse cDNA ends. The custom primer begins at the
unique boundary rather than consuming cycles on invariant adapter sequence.

The PCR endpoint is represented as a six-base-index member of the paper's stated
6–8-base index family. This makes both the custom 150-cycle Read 1 and the multiplexing
index-primer landing sites exact without pretending that every library used the same
index length or sequence.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5737731 JATS63–105, DMS treatment, RT, purification variants, T4DNA ligation and sequencing.

Added RNaseH digestion after SuperScriptIII55°C50min RT. Explicitly selected the PAGE branch and separated biotin/streptavidin and original-Structure-seq comparison variants. Confirmed published hairpin donor and single-read150.

Deferred: Published SupplementaryFigS1 and full primer geometry, two PAGE selection states, selected index width and ligation temperature wording; no new publication error asserted. This pass does not establish a complete source review or publication-wide absence.
