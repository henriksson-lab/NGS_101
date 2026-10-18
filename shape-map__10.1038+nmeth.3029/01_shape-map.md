# SHAPE-MaP

## 1. What it is

SHAPE-MaP measures RNA structural flexibility by chemically modifying flexible nucleotides and encoding those adducts as mutations during reverse transcription.

## 2. Sources and evidence

- 🟢 Siegfried NA et al. *RNA motif discovery by SHAPE and mutational profiling (SHAPE-MaP).* Nature Methods (2014). [doi:10.1038/nmeth.3029](https://doi.org/10.1038/nmeth.3029).
- 🟢 The defining paper's Online Methods (PMC4259394), sections “SHAPE-MaP using targeted gene-specific primers” and “Sequencing library preparation,” give the exact fixed portions of both targeted PCR rounds and the example small-RNA RT primer.
- 🟢 Busan S & Weeks KM. *Accurate detection of chemical modifications in RNA by mutational profiling (MaP) with ShapeMapper 2.* RNA (2018), open full text at https://pmc.ncbi.nlm.nih.gov/articles/PMC5769742/.
- 🟢 SHAPE electrophiles acylate RNA 2′-hydroxyls in a structure-dependent manner; MaP RT reads through adducts while introducing sequence changes.
- 🟢 Modified, untreated and denatured controls separate reagent reactivity from background and sequence-dependent RT error.
- 🟡 Library adapters depend on the chosen targeted or transcriptome implementation; the schematic now keeps the native fragmented-route molecular mixtures and leaves unassigned adapter/primer/read geometry explicit.

## 3. Oligos recovered from the defining paper

- 🟢 The small-RNA 3′-structure-cassette RT primer is `GAACCGGACCGAAGCCCG`. HIV-1 targeted reactions use a different gene-specific RT primer; the paper does not give one universal HIV-1 sequence.
- 🟢 Targeted PCR round 1 forward: `GACTGGAGTTCAGACGTGTGCTCTTCCGATC NNNNN [gene-specific primer]`.
- 🟢 Targeted PCR round 1 reverse: `CCCTACACGACGCTCTTCCGATCT NNNNN [gene-specific primer]`.
- 🟢 Targeted PCR round 2 forward: `CAAGCAGAAGACGGCATACGAGAT [Barcode] GTGACTGGAGTTCAGAC`.
- 🟢 Targeted PCR round 2 reverse: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG`.
- 🟢 The paper also describes a separate fragmented-cDNA route: second-strand synthesis, end repair, dA tailing, ligation of Illumina-compatible forked TruSeq adapters and emulsion PCR. It names the adapter family, but does not print an ordered adapter pair or all sample barcode identities. The schematic represents this native route and therefore must not imply that the targeted PCR primers are its reagents.

## 4. Why SHAPE-MaP

This is the foundational general-purpose mutational-profiling workflow with a defining method paper and an authoritative complete protocol. DMS-MaPseq is chemically narrower and should be a separate future page, not folded into this one.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC4259394 JATS103–113, fragmented-sample and targeted-primer SHAPE-MaP library routes.

Corrected graph/endpoint mismatch by choosing fragmented RNA with random nonamers and adding second-strand synthesis, end repair, dA tailing and TruSeq-compatible adapter ligation before30-cycle Q5 emulsion PCR. Targeted two-round PCR is no longer implied by this graph.

Deferred: Full neutral-RNA construct geometry, selected forked adapters, cleanup intermediates and run acquisition; targeted and viral-target-specific branches are outside this pass. This pass does not establish a complete source review or publication-wide absence.

## Native graph migration (2026-10-18, batch 12)

Inspected one original: published author-manuscript PDF `nmeth.3029_PMC4259394.pdf.txt`, Online Methods pp. 8–10, lines 420–512 (modification and fragmented-sample preparation), 513–542 (targeted route boundary and source-verified scraper hits), and 550–575 (ShapeMapper and quality trimming). The bundle oligo scraper was run and its source pointers reviewed; this does not establish source completeness.

The selected branch is fragmented ribosome/group II intron/HCV IRES RNA with random nonamer DNA primers. Original RNA 2′-O-adducts remain chemical substituents on RNA; independent primer-extended DNA encodes mutations without copying those chemical groups. Second-strand DNA, possible residual original RNA/oligos, and forked-adapter strands are independent molecular populations. Source-named G-50, PureLink and AMPure physical cleanup fractions are explicit; no class-specific recovery or purity is assumed.

Removed the queued shared TruSeq endpoint and its canonical sequencing/index defaults. The actual ligated mixture now supplies PCRGeometry, while ordered adapter strands, PCR primer sites, strand assembly and run cycles remain unassigned. PCR daughters and possible carried original templates, adapters and primers stay distinct. Quality trimming and demultiplexing are computational fates, not physical removal. Source Qubit/Bioanalyzer QC and pooling are recorded as a terminal limitation because actual aliquots, allocation volumes and companion libraries have not been reconstructed; no invented split or pool is drawn.

Deferred: historical ordered TruSeq adapter/primer originals; detailed NEB module internals and strand-specific end geometry; exact physical QC/pool allocations and run acquisition; HIV LNA and targeted-PCR alternatives; comprehensive source/linked-material absence review. This bounded migration changes no existing absence disposition.
