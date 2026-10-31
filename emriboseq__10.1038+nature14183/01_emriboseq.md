# emRiboSeq: capture a newly formed3′OH and sequence on Ion Torrent

Reijns et al. (2015), [Lagging-strand replication shapes the mutational landscape of the genome](https://doi.org/10.1038/nature14183).

## Capped sources

One reader,28October2026; article-associated accepted manuscript PMC4374164, complete library preparation/design and alignment Methods; original supplementary PDF tables1–3 (data sources, runs/strains). Full-bundle mandatory scraper run; printed trP1/A oligos and modifications verified in the article. Publisher extended-data PPT files have incomplete text conversions and remain unrendered; no chemistry absence is inferred. Later detailed emRiboSeq/EndoSeq protocols are deferred, not substituted for the original.

The original readout is **Ion Torrent PGM/Proton with OneTouch emulsionPCR**, not an Illumina library. The source explicitly reports no pre-sequencing amplification; PhusionFLASH is used for second-strand synthesis, not an invented library-PCR cycle series.

## Ordered original preparation

1. Sonicate yeast genomic DNA to about400bp, ethanol precipitate and1.2×AMPure select. NEBNext end-repair/dA-tail, AMPure purification; ligate trP1 adapter14–18h16°C. Subsequent bead capture uses NaCl/PEG conditions stated in the source.
2. TdT/ddATP blocks remaining free3′OH ends2h37°C, then AMPure purification. **Only subsequently** does recombinant RNaseH2 create new3′OH ends5′of embedded rNMPs; source Nb.BtsI cleavage is a validation alternative. Do not use alkaline cleavage or assume the rNMP itself becomes the captured terminal nucleotide.
3. Heat80°C20min,1.8×AMPure purify. Shrimp alkaline phosphatase removes5′phosphates at37°C, then65°C15min andAMPure. Heat95°C5min and snap-cool to expose single strands.
4. A-top's5′phosphate ligates to a newly available DNA3′OH using a six-random-base splint on A-bottom. A-top is3′dideoxy-blocked; A-bottom has5′biotin and3′C3 block. Ligate14–18h16°C.
5. StreptavidinM-280 captures biotinylated A-bottom-containing complexes. **0.15MNaOH releases the nonbiotinylated DNA strand** while biotinylated splint remains with beads. The biotin is on the splint, not covalently assigned to the retained captured DNA. Ethanol-precipitate the released single strand.
6. PrimerA/PhusionFLASH synthesizes the complement. Select200–300bp fragments using2%E-GelEX; BioanalyzerQC, OneTouch emulsionPCR, Ion Torrent PGM orProton sequencing.

## Source oligos,5′→3′

| Name | Printed sequence/modification |
|---|---|
| trP1-top | CCTCTCTATGGGCAGTCGGTGAT-phosphorothioate-T |
| trP1-bottom | Phosphate-ATCACCGACTGCCCATAGAGAGGC-dideoxy |
| A-top | Phosphate-CTGAGTCGGAGACACGCAGGGATGAGATGG-dideoxy |
| A-bottom | Biotin-CCATCTCATCCCTGCGTGTCTCCGACTCAGNNNNNN-C3phosphoramidite |
| PrimerA | CCATCTCATCCCTGCGTGTCTCCGAC |

Random6 on A-bottom is hybridization degeneracy, not a UMI. The printed trP1 pair has a blocked terminal base outside its maximal exact duplex; that base is retained as printed, not silently removed to force a full duplex. No sample-index sequence is supplied by resemblance to another Ion or Illumina protocol.

## Interpretation and deferreds

The source assigns a ribonucleotide site **one nucleotide upstream and on the opposite strand** to the aligned read5′end. That computational coordinate is distinct from the physical captured upstream fragment's3′OH. The captured fragment does not retain the downstream rNMP at its3′end. Because no pre-sequencing amplification was used, the source did not deduplicate genomic reads; this is not a UMI claim.

Exact complete original Ion library/strand assembly, emulsion-PCR and acquired run-primer/flow configuration remain unassembled in the partial model. Source-linked RNaseH2 purification and later detailed protocols, full extended-data graphics and a second independent reader are deferred. No missing-information assertion is accepted.

## Publication preparation review

The source-qualified public page now exposes actual NanoDrop/Bioanalyzer allocations and selected DNA stocks, the molecular ddATP contribution, and components of the same actual Phusion FLASH prepared master mix. Original captured-template3′dideoxy chemistry is retained through release and cleanup; it is not assigned to newly copied DNA. OneTouch supplies remain functional classes with unassembled whole primer polymers. Accepted-manuscript scope and inferred strand assembly remain explicit.
