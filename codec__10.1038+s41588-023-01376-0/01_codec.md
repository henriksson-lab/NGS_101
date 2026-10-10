# CODEC — research notes

Defining source: Ding *et al.*, [10.1038/s41588-023-01376-0](https://doi.org/10.1038/s41588-023-01376-0).

## Evidence used

- 🟢 The article defines CODEC as “Concatenating Original Duplex for Error Correction” and states that an adapter quadruplex plus strand-displacing extension physically links information from both original DNA strands in one Illumina read pair.
- 🟢 The read design recovers UMI, sample index and insert while retaining both complementary copies for a single-duplex consensus.
- 🔴 Exact oligo sequences were not available from the accessible article text. The schematic therefore shows adapter roles, not invented bases.

## Protocol boundary

This page covers CODEC library construction. Conventional Duplex Sequencing and NanoSeq remain separate protocols because they recover the two strands as independently amplified read families rather than concatenating them physically.
