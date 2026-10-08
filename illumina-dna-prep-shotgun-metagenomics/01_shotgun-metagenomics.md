# Illumina DNA Prep shotgun metagenomics

## 1. What it is

Illumina DNA Prep applies bead-linked transposomes to total DNA extracted from a microbial community. On-bead tagmentation fragments many genomes without locus selection, and a reduced-cycle PCR adds dual indexes and complete flow-cell arms. The resulting shotgun library samples both taxonomic and functional sequence content.

## 2. Authoritative sources

🟢 Illumina, *DNA Prep Reference Guide*, document 1000000025416 v12, defines bead-linked tagmentation, post-tagmentation wash and index PCR.

🟢 Illumina, *Shotgun metagenomics with Illumina DNA Prep on the NovaSeq 6000 System*, M-NA-00010 v2.0, identifies DNA extraction followed by Illumina DNA Prep as the demonstrated shotgun workflow.

## 3. Construct logic

1. Extract total DNA without selecting a marker locus.
2. Bead-linked Tn5 complexes fragment and tag the DNA. Immobilization makes the reaction saturable and provides a solid-phase wash.
3. Index PCR completes the P5/P7 arms and introduces the i5/i7 sample indexes.
4. Standard Nextera paired-end and index primers sequence randomly sampled community fragments.

## 4. What distinguishes it from amplicon profiling

There is no locus-specific primer and therefore no 16S or ITS boundary in the molecule. “Shotgun metagenomics” describes the biological sampling strategy; the library molecule itself is an Illumina DNA Prep/Nextera-style tagmentation product.

## 5. Proprietary regions

🟡 Current reagent mixes are commercial, but the final Nextera-compatible primer architecture is public. The schematic models the public adapter/read structure and treats the bead linkage as reaction geometry rather than DNA sequence.
