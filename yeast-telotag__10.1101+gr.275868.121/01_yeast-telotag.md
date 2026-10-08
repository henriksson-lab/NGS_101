# Yeast TeloTag nanopore sequencing

## 1. What it is

Sholes TeloTag marks the physical end of each *S. cerevisiae* telomere with poly(A) and a recognizable oligo before Oxford Nanopore whole-genome library construction.

## 2. Sources read

- 🟢 Sholes et al. (2022), [doi:10.1101/gr.275868.121](https://doi.org/10.1101/gr.275868.121), full text and methods.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), comparison with later long-read methods.

## 3. Construct and reaction order

1. 🟢 Terminal transferase adds poly(A) to native telomeric 3′ ends.
2. 🟢 One of five oligo(dT)-TeloTag primers anneals to that tail.
3. 🟢 Sulfolobus DNA polymerase IV fills the recessed strand without strand displacement; T4 ligase seals the nick.
4. 🟢 Oxford Nanopore adapters are added with the manufacturer workflow.

The exact five TeloTag sequences are in Supplemental Table S6, which the publisher endpoint did not provide automatically. The model therefore records their shared structure without invented bases.

## 4. Final library and sequencing

The TeloTag is a terminal landmark in an otherwise long genomic molecule. Nanopore sequencing uses a motor-loaded adapter, not a sequencing primer.
