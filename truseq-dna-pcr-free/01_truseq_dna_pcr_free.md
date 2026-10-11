# TruSeq DNA PCR-Free

## Sources read

🟢 Illumina, *TruSeq DNA PCR-Free Reference Guide*, document 1000000039279 v00, pp. 5–20 and 25. The guide defines Covaris fragmentation, end repair, bead size selection, 3′ adenylation, indexed-adapter ligation, library QC and pooling.

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, is the central source used by `lib/illumina.py` for TruSeq adapter and sequencing-primer sequences.

## Molecular workflow

1. Covaris shearing makes double-stranded fragments with 3′ or 5′ overhangs. The guide provides separate 350-bp and 550-bp insert workflows.
2. End repair makes blunt fragments, and a two-sided bead selection removes undesired sizes.
3. A-tailing adds one 3′ A to each repaired end.
4. Indexed TruSeq adapters with complementary T overhangs are ligated to both ends.
5. The library is quantified and pooled directly. There is no enrichment PCR, so qPCR is recommended to count cluster-forming molecules rather than all DNA.

## Final construct and readout

🟢 The guide lists single-index, combinatorial-dual-index and TruSeq unique-dual-index components. The schematic shows the current dual-index option.

🟡 Index identities vary with the purchased adapter plate, so the model records i5 and i7 as sample-index features without inventing particular bases. Standard TruSeq Read 1, Index 1, Index 2 and Read 2 primer landing sites are computed from the canonical construct.
