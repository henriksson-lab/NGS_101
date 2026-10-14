# ChIL-seq — research notes

Authoritative protocol: Handa *et al.*, [10.1038/s41596-020-0375-8](https://doi.org/10.1038/s41596-020-0375-8).

## Evidence used

- 🟢 A ChIL probe is an oligonucleotide-conjugated antibody that marks a chosen DNA-binding protein or chromatin modification in situ.
- 🟢 Tn5 transposition places an adapter in nearby genomic DNA; extension connects target-proximal sequence to the probe-derived promoter architecture.
- 🟢 T7 RNA polymerase provides linear amplification before cDNA sequencing-library construction.
- 🟢 Supplementary Table 4 prints the ChIL probe strands and modifications,
  the MEDS-B transposome oligos, the `Read2+N11` first-strand primer, and the complete
  `Ad1_noMX` and indexed `Ad2_index` library-PCR primers.
- 🟢 The terminal PCR primers create the standard single-index Nextera shell.
  Consequently the final Read 1, Index 1 and Read 2 primer sites are computed from that
  completed library rather than left as unnamed roles.

## Protocol boundary

ChIL-seq is retained separately from CUT&Tag because its antibody oligo and T7 transcription are part of the defining molecular path.
