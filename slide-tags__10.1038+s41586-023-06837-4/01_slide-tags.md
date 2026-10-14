# Slide-tags — research notes

Defining source: Russell et al., 2024, [doi:10.1038/s41586-023-06837-4](https://doi.org/10.1038/s41586-023-06837-4). The 2023 preprint DOI is supporting history, not a second protocol.

## Sources read

- 🟢 PMC article and Methods.
- 🟢 Supplementary Methods sections 1.3–1.6: bead synthesis, photocleavage, tissue/nucleus recovery, 10x 3-prime and Multiome branches, spatial-tag PCR and sequencing.
- 🟢 Published article PMC10764288 and Supplementary Methods MOESM1 DOCX, sections 1.3 and 1.5; linked preprint 10.1101/2023.04.01.535228 agrees on the oligos below. The DOCX is used for the ligated-bead oligo because the published JATS text introduces an extra `z` at that location that is absent from the DOCX and preprint.

## Molecular path

🟢 A sequenced puck supplies spatial barcode oligos through a photocleavable linker. Released oligos enter nearby nuclei before dissociation. Standard 10x processing supplies a cell barcode; a feature-barcode-style library connects that cell barcode to one or more spatial barcodes.

🟢 The protocol supports snRNA-seq and joint snRNA/snATAC. The schematic shows RNA, ATAC and spatial-tag outputs while stating that ATAC belongs to the Multiome branch.

## Source-printed oligos

🟢 The published Supplementary Methods specify three puck designs. `PC` denotes the photocleavable linker, `N` a random base and `V` a non-T base under IUPAC notation; the source prints `J` but does not define its chemistry in this passage. The directly synthesized capture oligo is printed as `5′-TTT-PC-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVVGCTTTAAGGCCGGTCCTAGCAA-3′`. Its poly(A) variant ends after `NNNNNNNVV` in `A30`. For the ligated puck design, the source prints a different 5′ arm and the exact extension oligo `GAGCTTTGCTAACGGTCGAGGCTTTAAGGCCGGTCCTAGCAA`, annealed with the splint shown 3′→5′ as `CTGTGTACCCGCCTCGAAACGATTGC`.

🟢 The full ligated-design bead oligo in MOESM1 is `5′-TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVVGCTCGGACACATGGGCG-3′`. The published JATS extraction's `z` immediately after `_PC_` is an extraction discrepancy, not supported by either Supplementary Methods DOCX or the preprint PDF.

🟢 The Multiome branch adds spike-in primer `GTGACTGGAGTTCAGACGT`. Spatial-tag PCR uses `STAG_iP7_a1`, `CAAGCAGAAGACGGCATACGAGATATTTACCGCAGTGACTGGAGTTCAGACGT*G*T`, and `P5-STAG_ip5_a1`, `AATGATACGGCGACCACCGAGATCTACACGACAATAAAGACACTCTTTCCCTACACGACGC*T*C`; asterisks mark the printed terminal phosphorothioate bonds. The 3′ RNA branch instead uses the named `STAG_P701_NEX` and `P5-Truseq Hybrid` primers, whose bases are not printed in these Methods passages.

🟡 The public schematic still uses role tokens for the variable puck barcodes and commercial 10x bead branches. The published fixed puck arms and custom Multiome PCR primers are known; those roles should not be described as wholly unavailable.
