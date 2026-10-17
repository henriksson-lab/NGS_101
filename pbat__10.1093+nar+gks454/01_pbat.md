# PBAT — research notes

## 1. What it is

Post-bisulfite adaptor tagging reverses conventional WGBS order: destructive conversion
comes first, and two random-primer extensions install the sequencing handles afterward.

## 2. Sources read

- 🟢 Miura et al. 2012, DOI [10.1093/nar/gks454](https://doi.org/10.1093/nar/gks454),
  complete article and supplement route (PMC3458524).

## 3. Molecular order

🟢 Bisulfite-treated single-stranded DNA is copied from a 5′-biotinylated adapter primer
ending in random N4. Two AMPure purifications remove residual primers and primer-dimers before
streptavidin capture; alkaline washes then expose the immobilized first strand. 🟢 A second
adapter primer, also ending in N4, makes the strand released for sequencing and size
selection. The original demonstration supports PCR-free library production.

## 4. Evidence boundary

🟢 N4 architecture, biotin capture, extension order and directionality come from the
paper. 🟢 The main paper prints BioPEA2N4 as
`5′-biotin-ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNN-3′`, PE-reverse-N4 as
`5′-CAAGCAGAAGACGGCATACGAGATNNNN-3′`, and Primer-3 as
`5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`.
The original readout is unindexed single-end sequencing of the generally G-poor strand
complementary to bisulfite-converted DNA.
🟢 The article also prints PE-forward `AATGATACGGCGACCACCGAGATCTACAC` and
PE-reverse `CAAGCAGAAGACGGCATACGAGAT` for library quantification. These are
measurement primers, not a required global PCR step in the amplification-free
library path. The supplement additionally gives locus-specific qPCR controls,
which do not become part of the sequencing construct.

## 5. Bounded source pass — 2026-10-17

Read the published PMC3458524 JATS Methods lines 28–60, strand interpretation
64–66, method overview 86–90 and mouse readout 138; inspected published
File006 supplement text, Supplementary Figures S1–S3 captions lines 1–35,
qPCR context 100–108 and Table S1 beginning 110–138. Source manifest identifies
these as the published article and supplement. Ran the bundle oligo scraper
and checked the complete three library primers and the BioPEA2N4 biotin
against article lines 40, 48 and 52. Scraper fragments from PDF extraction are
not independent oligo identities.

🟢 First-strand extension uses Klenow fragment (3′→5′ exo−). Two AMPure
purifications precede Dynabeads M280 streptavidin capture and alkaline washes.
Second-strand priming also uses Klenow; after its heat inactivation, Bst DNA
polymerase large fragment completes second-strand synthesis. Primer-3 extension
uses Phusion Hot Start High-Fidelity DNA Polymerase. The graph previously
collapsed these catalysts into a single generic polymerase and omitted Bst.
They are now separate, source-named stages. Primer inputs have source-derived
strand states; oligo ordering text preserves the biotin modification.

🟢 Following Primer-3 extension, the supernatant is transferred, treated with
Exonuclease I and purified twice with AMPure XP. The source reports removal
of fragments smaller than 200 bp. These steps now precede the graph's sequencing
terminal. The model does not claim complete representation of the Exonuclease I
hydrolysis products, bead-associated leftover strands or recovery efficiencies.

🟢 Article lines 64–66 report occasional opposite-strand C-poor reads; File006
Figure S3's caption proposes genomic self-priming and displacement as a plausible
mechanism. The public caveat retains this uncertainty; a proposed pathway is
not established reaction chemistry. No global library PCR is introduced by
these changes; qPCR remains a measurement operation outside the selected path.

Deferred: rendered S3 geometry and remaining supplement tables/figures;
historical kit/instrument originals and DNA preparation predecessor; explicit
bisulfite purification/desulfonation and bead-release states; molecular
opposite-strand branch and a dataset-specific acquired cycle declaration.
The paper reports multiple single-end datasets; primer geometry alone does
not select their lengths. This capped pass does not establish publication
absence or complete the two-reader gate.
