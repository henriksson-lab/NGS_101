# CELLO-seq — research notes

## Sources read

- 🟢 Berrens et al. 2022, DOI `10.1038/s41587-021-01093-1`, PMC article, Methods and Supplementary Table 5.
- 🟢 Berrens et al. 2025, DOI `10.1038/s41596-025-01203-2`, PMC detailed protocol and supplementary materials.

## Construct path

🟢 Poly(A)-primed, high-temperature RT uses a 3′-amino-blocked TSO. After Exonuclease I cleanup and denaturation, a pre-annealed splint is ligated to first-strand cDNA at 55 °C. It adds the cellular barcode and a 22-base UMI designed as repeated RYN classes to avoid long homopolymers. 🟢 High PCR duplicate counts are deliberate: noisy ONT reads sharing one long UMI are grouped to form a consensus suitable for locus-specific transposable-element assignment.

## Original oligos and version boundary

🟢 Supplementary Table 5 of the defining paper prints the original phosphorylated oligo(dT), amino-blocked TSO, top splint, and all 96 bottom splints. The page models the first listed bottom splint (`splintlig_spc_idx1`) explicitly. It carries, in order, the Nextera R2-derived tail, a 24-base cell barcode, `AGTGGTATC`, a 22-base patterned UMI (`NRY` repeated seven times plus `N`), the reverse complement of the top splint, and a terminal `T14` tract. The top splint has a 3′ inverted-dT block.

🟢 The later 2025 Nature Protocols article changes this design: its bottom splint uses a 50-base UMI and the article explicitly describes that as replacing the original 22-base UMI. Those oligos define a revised implementation, not missing detail that can be copied into the 2022 construct. The public schematic therefore stays with the defining paper's 22-base design and links the revised protocol separately.

🟡 The published oligos determine the internal splint duplex and identifier-bearing product. The terminal Oxford Nanopore motor-loaded adapter is named by the workflow, but its bases and end geometry are not disclosed in the sources read, so that terminal part remains inferred rather than being invented.
