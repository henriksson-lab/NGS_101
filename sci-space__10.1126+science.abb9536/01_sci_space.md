# sci-Space — research notes

Srivatsan et al., *Science* (2021),
[doi:10.1126/science.abb9536](https://doi.org/10.1126/science.abb9536).

## Sources read

- 🟢 Defining article (PMCID: PMC9118175); the 77.5 MB supplementary-methods PDF is listed by PMC but its download is gated in the source cache.
- 🟢 Author-linked protocols.io guides for [grid preparation](https://doi.org/10.17504/protocols.io.bm64k9gw), [spatial transfer and imaging](https://doi.org/10.17504/protocols.io.btw5npg6), and [two-level sci-RNA-seq library preparation](https://doi.org/10.17504/protocols.io.bquamwse). Their fetched PDFs were checked with the oligo scraper and read at the cited steps.
- 🟢 Cao et al. two-level sci-RNA-seq chemistry,
  [doi:10.1126/science.aam8940](https://doi.org/10.1126/science.aam8940).

## What is spatially encoded

🟢 The agarose-coated slides carry a printed 18 mm × 18 mm grid of 7,056 uniquely
barcoded spots. Hash oligos are single-stranded DNA with poly(A) tails. About 5% of the
spots also contain SYBR Green, forming fluorescent fiducials used to register the ideal
grid to the tissue image. The reported mean spot radius was 73.2 µm and center spacing
222 µm.

🟢 A separate hash identifies the tissue section. It permits sections to be pooled,
but is not itself the x–y coordinate. 🟡 The linked transfer guide calls for a unique
slide-specific oligo in lysis buffer (step 1.3) without printing its sequence or tail;
its poly(A) structure cannot yet be verified from the readable sources.

🟢 The grid-preparation guide, step 7.1, prints the spot/sector hash ordering structure
verbatim as `5′-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG[10bp-barcode sequence]-BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-3′`,
where `B` is C, G or T. The 10-base barcode differs among ordered oligos. The guide
points to `all_plates_for_ordering.xlsx` for the complete panel, but that attachment is
not included in the fetched PDF. The printed structure is also shown on the public page;
the individual plate sequences remain an open source item, not a confirmed absence.

## Construct-changing steps

1. 🟢 Cryosection fresh-frozen tissue, permeabilize it with the section-hash oligo and
   juxtapose it with the printed agarose grid. Grid hashes diffuse into nearby nuclei.
2. 🟢 Image the sandwich during transfer so fluorescent fiducials preserve the transform
   between barcode coordinates and tissue anatomy.
3. 🟢 Extract nuclei, fix them and perform sci-RNA-seq. The barcoded oligo(dT) RT primer
   captures endogenous poly(A) RNA and the printed hashes in the same nucleus. The
   section hash is also recovered with the nucleus, but its exact capture-ready tail
   structure is not established by the readable guide.
4. 🟢 Subsequent pooling, redistribution, second-strand synthesis, Tn5 fragmentation and
   indexed PCR create separate transcript and hash molecules carrying the same nucleus
   identity.

🟢 The author-linked two-level sci-RNA-seq guide prints the terminal amplification
primers as `5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGA-3′`
and `5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′` (PCR section,
extracted lines 579–602). It describes an n7-loaded Tn5 or dual-loaded Illumina Tn5
before PCR (steps 37–39). These are linked predecessor-guide statements; the main
sci-Space paper does not itself spell out every library-primer sequence.

## Library boundary

🟢 A spatial hash is not covalently fused to a cell’s transcripts. Association is by the
shared sci-RNA-seq combinatorial indices: transcript reads define expression and hash
reads define the most enriched printed coordinate and section. The schematic therefore
shows a representative hash molecule and separately states the transcript branch rather
than drawing a chemically false combined molecule.

🟡 The public diagram uses the established two-level sci-RNA-seq outer primer geometry.
Individual barcode identities are represented as semantic regions because the experiment
uses large panels; they are not synthesized as strings of repeated placeholder bases.
