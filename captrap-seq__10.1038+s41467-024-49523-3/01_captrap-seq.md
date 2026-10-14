# CapTrap-seq — research notes

## Sources read

- 🟢 Glinos et al. 2024, DOI `10.1038/s41467-024-49523-3`, accepted article supplements and the author's bioRxiv full text and Table S1 (`10.1101/2023.06.16.543444`).

## Construct path

🟢 Anchored oligo(dT) makes first-strand cDNA from polyadenylated RNA. The intact RNA cap is biotinylated and streptavidin capture retains hybrids whose cDNA reached the 5′ end. RNase treatment releases the cDNA. 🟢 A second full-length selection requires sequential 5′- and 3′-end linker ligation before long-and-accurate PCR. The same selected cDNA pool was adapted for ONT and PacBio; the page does not conflate those platform adapters with the CapTrap selection chemistry.

## Oligos from the author preprint

🟢 Table S1 prints the following ordering sequences. Its `P` denotes phosphate, `U` denotes a uracil residue, `V/N` are degenerate bases, and `N` in the second-strand primer is an 8-base UMI. The 5′/3′ linker rows are separate strands, shown in their printed 5′→3′ orientation; source line wraps are rejoined here.

- OligodT 16VN: `5′-P-TTTTTTTTTTTTTTTTVN-3′`.
- 5′ linker strand 1: `5′-GTGGTAUCAACGCAGAGUACGNNNNN-P-3′`.
- 5′ linker strand 2: `5′-P-CACCATAGTTGCGTCTCATG-P-3′`.
- 3′ linker strand 1: `5′-AAAAAGCAUCGCUGTCTCUTAUACACAUCUCCGAGCCCACGAGAC-P-3′`.
- 3′ linker strand 2: `5′-CGTAGCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG-3′`.
- Second-strand UMI8 primer: `5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNGTGGTATCAACGCAGAGTA-3′`.
- Forward LA-PCR primer: `5′-TCGTCGGCAGCGTC-3′`.
- Reverse LA-PCR primer: `5′-GTCTCGTGGGCTCGG-3′`.

🟢 The Methods specify sequential ligation of these linkers, SAP/USER treatment to convert the duplexes to single-strand handles, second-strand synthesis with UMI8, and LA-PCR with the printed pair. 🟡 The precise linker-to-cDNA junctions still need strand-resolved reconstruction from the source figures; the schematic does not assign an unverified overhang. 🔴 The named ONT motor adapter and PacBio platform hairpin bases are not disclosed by the article.
