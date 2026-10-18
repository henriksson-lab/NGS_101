# YAMAT-seq

**Defining source:** Shigematsu *et al.*, *Nucleic Acids Research* (2017),
doi:10.1093/nar/gkx005. The PMC full text and supplement were read.

## Evidence used

- 🟢 Methods lines 52–60 publish all five Y-adapter oligos, their DNA/RNA composition,
  end modifications, Rnl2 annealing/ligation, TruSeq Small RNA RT-PCR and 100-nt
  single-read sequencing.
- 🟡 **Source discrepancy:** Methods, “Y-shaped adapter design” (PMC JATS line 52),
  prints Y-3′-AD as `5′-/5Phos/GTATCCAGTTGGAATTCTCGGGTGCCAAGG/3ddC/-3′`,
  with a blocked 3′ end. Results, “Y-AD” description (line 88), calls the
  same strand 3′-OH. The model follows the explicit ordered oligo `/3ddC/`;
  the physical oligo used cannot be verified from these conflicting statements.
- 🟢 Lines 80–88 establish deacylation, recognition of the mature-tRNA 5′-NCCA-3′
  overhang and sealing of both acceptor-stem nicks.
- 🟢 The defining Methods names the TruSeq Small RNA RTP, PML and kit PCR primers. The current source pass has not assigned an applicable historical kit member or its index width. An earlier note quoted a later vendor sequence sheet; that reference is historical and was not independently re-read in this capped pass. It does not fix the native 2017 kit or run.
- 🟡 The kit-derived library ends and index remain inferred roles, without a six-base default. The native RNA/DNA adapter sequences are retained independently of this kit boundary.
- 🟢 The additional nucleotide sequences in the paper's methods (lines 64–68 of the
  PMC full-text extraction) are a synthetic tRNA substrate and RT/qPCR primers for
  **ligation-efficiency experiments**. They are not the library amplification primers
  used for the sequencing runs; that step explicitly uses primers supplied with the
  TruSeq Small RNA kit (line 60). This distinction prevents the assay-control primers
  from being substituted for the unresolved sequencing-primer landing site.

## Construct logic

The four Y-5′ variants differ only in their terminal ribodinucleotide and collectively
pair with all discriminator-base/NCCA ends. Y-3′-AD is 5′ phosphorylated and 3′ blocked.
Annealing creates two Rnl2 substrates; the source describes selective amplification of molecules with the two relevant adapter–tRNA boundaries sealed; this does not establish a physical discard fate for all other molecules.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5605243 JATS46–60, deacylation, adapter design, ligation, RT and sequencing.

Specified deacylation at pH9/37°C/40min, 37°C1h plus overnight4°C ligation, SuperScriptIII55°C60min and 11-cycle library PCR; recorded single-read100. Printed /3ddC/ and the existing Methods/Results discrepancy notice retained.

Deferred: Original Results discrepancy reread, selected kit primer geometry, exact acceptor-stem annealing and gel-selection intermediates; not a new publication error finding. This pass does not establish a complete source review or publication-wide absence.

Corrected the adapter-set state to five separate oligo strands instead of one concatenation, fixing a real source-guard failure.

## Native graph migration, 2026-10-18 batch07

Two originals inspected: published PMC5605243 Methods46–68 and Results80–95 (deacylation, five oligos/modifications, annealing/ligation, RT/PCR, PAGE and single-read100); published supplement opening FigS1–S3 headings/captions and controls. The mandatory bundle oligo scraper was run and relevant five-adapter hits checked against the original Methods. No additional vendor original was read in this pass.

The protocol now owns the entire graph: deacylation followed by ethanol recovery, separate adapter molecules, noncovalent terminal-pairing state, both Rnl2 junctions on a covalent RNA/DNA strand, SSIII RT on that actual strand, eleven-cycle PCR on the first-strand population, and 8% native PAGE. Y5-U retains physical terminal RNA U rather than DNA T. The ordered `/3ddC/` contributes an additional physical terminal C, drawn separately and marked dideoxy; it is not just an inert note on the preceding G. The Methods5P/3ddC Y3 ordering specification remains on the public page beside the conflicting Results3OH statement.

The pairing view computes the source Y-adapter stem and the A-discriminator/Y5-U representative terminal interaction. It does not claim to position the complete folded acceptor stem or both nick geometries. All four terminal classes are retained as the reagent population; the representative is not a chosen biological tRNA sequence. The fully ligated strand retains terminal ddC, consumed5P junctions and its original RNA through RT. Residual original RNA/first-strand template recovery after PCR and gel selection stays unresolved. Bioinformatic CCA/length filtering does not discard physical molecules.

Deferred: complete acceptor-fold layout; applicable original kit primer/member definitions and native index width; precise first-strand boundaries, sequencing landing and class-specific recovery. This capped migration does not establish confirmed missing information or alter historical review statuses.
