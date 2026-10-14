# CODEC — research notes

Defining source: Bae *et al.*, [10.1038/s41588-023-01376-0](https://doi.org/10.1038/s41588-023-01376-0).

## Evidence used

- 🟢 The article defines CODEC as “Concatenating Original Duplex for Error Correction” and states that an adapter quadruplex plus strand-displacing extension physically links information from both original DNA strands in one Illumina read pair.
- 🟢 The read design recovers UMI, sample index and insert while retaining both complementary copies for a single-duplex consensus.
- 🟢 Supplementary Table 5 prints all four CODEC adapter oligos for eight matched
  sample-index sets. It defines `*` as a phosphorothioate linkage and `P-` as a 5′
  phosphate. The schematic uses set 1 as a representative molecule and preserves both
  modifications in its ordering list.
- 🟢 The inventors' [US20240052342A1 sequence listing](https://patents.google.com/patent/US20240052342A1/en)
  independently prints the set-1 `LD4-adap5`, `LD4-adap7`, `LD4-brid5` and
  `LD4-brid7` oligos, with `*` and 5′ `P-` notation matching the model. This is
  useful when the publisher's supplementary PDF cannot be downloaded.
- 🟢 `LD4-brid5` contains the reverse complement of the standard TruSeq Read 1
  primer site; `LD4-brid7` contains the reverse complement of the standard TruSeq Read 2
  primer site. Each is preceded by a 3-nt UMI and an 18-nt inline sample index, consistent
  with the outward Read 1/Read 2 geometry in Fig. 1b and the CODECSuite description.
- 🟡 The paper draws the post-extension linked product by region rather than
  printing it as one contiguous strand. That state is therefore role-scaled; it is not a
  claim that the adapter oligos themselves are unavailable.

## Protocol boundary

This page covers CODEC library construction. Conventional Duplex Sequencing and NanoSeq remain separate protocols because they recover the two strands as independently amplified read families rather than concatenating them physically.
