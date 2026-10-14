# SHAPE-MaP

## 1. What it is

SHAPE-MaP measures RNA structural flexibility by chemically modifying flexible nucleotides and encoding those adducts as mutations during reverse transcription.

## 2. Sources and evidence

- 🟢 Siegfried NA et al. *RNA motif discovery by SHAPE and mutational profiling (SHAPE-MaP).* Nature Methods (2014). [doi:10.1038/nmeth.3029](https://doi.org/10.1038/nmeth.3029).
- 🟢 The defining paper's Online Methods (PMC4259394), sections “SHAPE-MaP using targeted gene-specific primers” and “Sequencing library preparation,” give the exact fixed portions of both targeted PCR rounds and the example small-RNA RT primer.
- 🟢 Busan S & Weeks KM. *Accurate detection of chemical modifications in RNA by mutational profiling (MaP) with ShapeMapper 2.* RNA (2018), open full text at https://pmc.ncbi.nlm.nih.gov/articles/PMC5769742/.
- 🟢 SHAPE electrophiles acylate RNA 2′-hydroxyls in a structure-dependent manner; MaP RT reads through adducts while introducing sequence changes.
- 🟢 Modified, untreated and denatured controls separate reagent reactivity from background and sequence-dependent RT error.
- 🟡 Library adapters depend on the chosen targeted or transcriptome implementation; the schematic declares an inferred standard Illumina endpoint.

## 3. Oligos recovered from the defining paper

- 🟢 The small-RNA 3′-structure-cassette RT primer is `GAACCGGACCGAAGCCCG`. HIV-1 targeted reactions use a different gene-specific RT primer; the paper does not give one universal HIV-1 sequence.
- 🟢 Targeted PCR round 1 forward: `GACTGGAGTTCAGACGTGTGCTCTTCCGATC NNNNN [gene-specific primer]`.
- 🟢 Targeted PCR round 1 reverse: `CCCTACACGACGCTCTTCCGATCT NNNNN [gene-specific primer]`.
- 🟢 Targeted PCR round 2 forward: `CAAGCAGAAGACGGCATACGAGAT [Barcode] GTGACTGGAGTTCAGAC`.
- 🟢 Targeted PCR round 2 reverse: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG`.
- 🟢 The paper also describes a separate fragmented-cDNA route: second-strand synthesis, end repair, dA tailing, ligation of Illumina-compatible forked TruSeq adapters and emulsion PCR. It names the adapter family, but does not print an ordered adapter pair or all sample barcode identities. The schematic currently represents this route's endpoint and therefore must not imply that the targeted PCR primers are its reagents.

## 4. Why SHAPE-MaP

This is the foundational general-purpose mutational-profiling workflow with a defining method paper and an authoritative complete protocol. DMS-MaPseq is chemically narrower and should be a separate future page, not folded into this one.
