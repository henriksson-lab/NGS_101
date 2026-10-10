# Fiber-seq — research notes

## Sources read

- 🟢 Stergachis et al. 2020, DOI [10.1126/science.aaz1646](https://doi.org/10.1126/science.aaz1646), defining paper.
- 🟢 Fiber-seq bench protocol, linked by the authors' Fiber-seq documentation; current Hia5 and PacBio workflow.
- 🟢 PacBio SMRTbell prep kit 3.0 procedure, for the current PCR-free dumbbell endpoint.

## Molecular order

🟢 Hia5 methylates solvent-accessible adenines in intact nuclei using SAM. Protein-bound
DNA is protected, so the resulting m6dA pattern records accessibility and footprints on
individual native molecules. High-molecular-weight DNA is extracted without PCR. The
current workflow shears to long inserts, repairs and dA-tails the ends, ligates PacBio
hairpins, removes open molecules enzymatically, then anneals sequencing primer and binds
polymerase for circular-consensus sequencing.

## Interpretation boundary

🟡 The schematic shows representative accessible and protected intervals, not a fixed
nucleosome spacing. 🟡 The defining 2020 experiment predates SMRTbell prep kit 3.0; the
endpoint shown is the current standard Fiber-seq library workflow, linked as an additional
protocol source. Hairpin sequences are proprietary.
