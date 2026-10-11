# sci-Space — research notes

Srivatsan et al., *Science* (2021),
[doi:10.1126/science.abb9536](https://doi.org/10.1126/science.abb9536).

## Sources read

- 🟢 Defining article and supplementary methods (PMCID: PMC9118175).
- 🟢 Cao et al. two-level sci-RNA-seq chemistry,
  [doi:10.1126/science.aam8940](https://doi.org/10.1126/science.aam8940).

## What is spatially encoded

🟢 The agarose-coated slides carry a printed 18 mm × 18 mm grid of 7,056 uniquely
barcoded spots. Hash oligos are single-stranded DNA with poly(A) tails. About 5% of the
spots also contain SYBR Green, forming fluorescent fiducials used to register the ideal
grid to the tissue image. The reported mean spot radius was 73.2 µm and center spacing
222 µm.

🟢 A second polyadenylated hash identifies the tissue section. It permits sections to be
pooled, but is not itself the x–y coordinate.

## Construct-changing steps

1. 🟢 Cryosection fresh-frozen tissue, permeabilize it with the section-hash oligo and
   juxtapose it with the printed agarose grid. Grid hashes diffuse into nearby nuclei.
2. 🟢 Image the sandwich during transfer so fluorescent fiducials preserve the transform
   between barcode coordinates and tissue anatomy.
3. 🟢 Extract nuclei, fix them and perform sci-RNA-seq. The barcoded oligo(dT) RT primer
   captures endogenous poly(A) RNA, grid hashes and section hashes in the same nucleus.
4. 🟢 Subsequent pooling, redistribution, second-strand synthesis, Tn5 fragmentation and
   indexed PCR create separate transcript and hash molecules carrying the same nucleus
   identity.

## Library boundary

🟢 A spatial hash is not covalently fused to a cell’s transcripts. Association is by the
shared sci-RNA-seq combinatorial indices: transcript reads define expression and hash
reads define the most enriched printed coordinate and section. The schematic therefore
shows a representative hash molecule and separately states the transcript branch rather
than drawing a chemically false combined molecule.

🟡 The public diagram uses the established two-level sci-RNA-seq outer primer geometry.
Individual barcode identities are represented as semantic regions because the experiment
uses large panels; they are not synthesized as strings of repeated placeholder bases.
