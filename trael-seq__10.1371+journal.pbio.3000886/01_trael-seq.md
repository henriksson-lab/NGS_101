# TrAEL-seq — research notes

## Sources read

- 🟢 Kara et al. 2021, DOI `10.1371/journal.pbio.3000886`, PMC article, protocol supplement `pbio.3000886.s009` and oligo section.
- 🟢 NEB E7335 Multiplex Oligos manual, version 8.0 (June 2024), primer table on page 4; the defining paper names E7335S for indexed library PCR.

## Construct path

🟢 A 5′-adenylated hairpin adapter is ligated to native DNA 3′-OH ends with truncated RNA ligase in the absence of ATP. It contains eight random bases, internal biotin and uracils. Bst polymerase fill/displacement stabilizes the junction. After sonication and second-adapter ligation, streptavidin capture and USER cleavage release end-derived molecules for PCR.

🟢 The defining paper prints adapter 1 in full as
`[Phos]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTUGCGCAGGCCATTGGCC[BtndT]GCGCUACACTCTTTCCCTACACGACGCT` and states that this synthesized oligo is subsequently 5′ adenylated. The schematic therefore shows the ligation-ready form as 5′-App and folds the complete 23-bp stem instead of abbreviating its loop. Adapter 2 is likewise printed in full, including four deoxyuridines and its terminal phosphorothioate bond.

🟢 The paper's named E7335S PCR set uses one six-base i7 index and an index-free Universal PCR primer, as confirmed by the NEB primer table. The earlier schematic incorrectly inherited two eight-base indexes from a generic TruSeq shell. The corrected library shows only Index 1 (six cycles). The paper reports 50-base MiSeq single reads or 75-base NextSeq 500 single reads, but does not name the selected E7335 index primer for each library. 🔴 The per-library index identity remains unverified while linked author documents and supplementary originals remain unread.
