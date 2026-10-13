# Quantum Barcoding (QBC) — source notes

**Defining source:** Roy et al., *Communications Biology* (2020),
[doi:10.1038/s42003-020-0896-2](https://doi.org/10.1038/s42003-020-0896-2).

## What the supplement establishes

- 🟢 Fixed cells can carry targeted RNA cDNA and antibody-derived reporters.
- 🟢 Four successive 60-well split–anneal–ligate–wash–pool rounds add SC3, SC2,
  SCN and SC1. Each contributes a 7-nt cell subcode; SC3 also contributes a 10-nt UMI.
- 🟢 The published Splint6 molecule holds the four junctions on one scaffold and
  separates its binding regions with internal C3 spacers.
- 🟢 Antibody reporters use a 9-nt antibody-handle code (AHCA); targeted RNA instead
  carries target-specific cDNA after the common anchor.
- 🟢 The final antibody amplicon alignment in Supplementary Data 2 is 204 nt, of which
  the marked sequenced block is 133 nt. NEBNext indexed primers add Illumina ends.
- 🟡 Run lengths varied (151–258 bases, single- or paired-end), so the page describes
  stable molecular content rather than inventing one mandatory cycle recipe.

The source prints the oligos in several orientations. The model stores every ordered
oligo 5′→3′ and derives sequencing placement from the completed library.
