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

🟢 The ligated adaptor is full-length, which is what makes the workflow PCR-free: the guide's ligation step "ligates index adapters to the DNA fragment ends, preparing them for hybridization to a flow cell", and its ligation control is described as joining "3′-T overhang adapters to 3′-A overhang inserts". *Illumina Adapter Sequences* prints both arms of that adaptor — `Index 1 (i7) Adapters GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG` and `Index 2 (i5) Adapters AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT` (single-index kits instead use the one-oligo `TruSeq Universal Adapter`) — so P5, P7 and every sample index are present before ligation and no primer installs them afterwards.

🟢 The same document writes the dA/dT junction into the oligo listing itself: "the additional A base is in parentheses in the i7 adapter", `(A)GATCGGAAGAGC…`. The model therefore draws a ligation junction at both insert ends: the adaptor's own 3′ dT on the left, the insert's 3′ dA on the right.

🟡 Because the adaptor is forked, the sealed molecule's two platform arms are unpaired beyond the 12-bp TruSeq stem until the flow cell copies them. The reagent's fork geometry is modelled in `lib/batch_ngs.py` (`FullLengthAdapter.scene()`), while the library construct keeps the project's linear one-top-strand convention used for every other ligation product.

🟢 The guide lists single-index, combinatorial-dual-index and TruSeq unique-dual-index components. The schematic shows the current dual-index option.

🟡 Index identities vary with the purchased adapter plate, so the model records i5 and i7 as sample-index features without inventing particular bases. Standard TruSeq Read 1, Index 1, Index 2 and Read 2 primer landing sites are computed from the canonical construct.
