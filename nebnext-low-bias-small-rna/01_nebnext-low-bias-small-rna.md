# NEBNext Low-bias Small RNA

## 1. What it is

A gel-free small-RNA library protocol that uses randomized splint ligation to reduce sequence-dependent adaptor ligation bias.

## 2. Sources and evidence

- 🟢 [NEB E3420 manual v2.0 (May 2026)](https://www.neb.com/en/-/media/nebus/files/manuals/manuale3420.pdf), read from `_data/sources/nebnext-low-bias-small-rna/manuale3420.pdf.txt`, and the [product page](https://www.neb.com/en-us/products/e3420-nebnext-low-bias-small-rna-library-prep-kit).
- 🟢 Captured RNAs are under 120 nt and require a 5′ monophosphate and 3′ hydroxyl.
- 🟢 The order is 3′-adaptor ligation, excess-adaptor removal, simultaneous 5′ ligation plus modification of the 3′ adaptor into an RT substrate, RT, bead selection and UDI PCR.
- 🟢 NEB's [E3420 RT-primer FAQ](https://www.neb.com/en-us/faqs/when-is-the-rt-primer-annealed-during-the-nebnext-low-bias-small-rna-library-prep-workflow-neb-e3420) specifies that the bottom strand of the 3′ adaptor becomes the RT primer during 5′ ligation. The kit does not add a separate RT-primer oligo.
- 🟢 The [NEB E3420 product page](https://www.neb.com/en/products/e3420-nebnext-low-bias-small-rna-library-prep-kit) cites [Maguire et al. (2020)](https://doi.org/10.1093/nar/gkaa480) as the research predecessor and explicitly says that the commercial kit has subsequent protocol enhancements. Maguire's Supplementary Table S2 prints prototype oligos 4–7: phosphorylated DNA 3′ adapter, uracil-bearing six-base-randomized 3′ splint with a 3′ inverted dT, RNA 5′ adapter, and a 5′ splint. The article describes preadenylation and splint annealing before 3′ then 5′ ligation. These are source-positive facts for the research method; no checked source establishes that E3420 uses identical oligos or modifications.
- 🔴 The checked E3420 manual names LB 3′ and 5′ adaptors and NEBNext LV Unique Dual Index Primer Pairs, but does not print complete kit adaptor or randomized-splint ordering strands and modifications. The UDI primer pairs are purchased separately as sets E3390/E3392/E3400/E3402/E3404/E3406/E3408; no particular pair belongs to the generic workflow. NEB's patent predecessors and the remaining product FAQs are still under review, so kit-specific sequence absence has not been accepted.

## 3. Molecular path

The protocol directly ligates RNA ends. It therefore does not inherit the dA-tailed DNA junction used by Ultra II DNA workflows. A typical 56-base Read 1 is sufficient for the small-RNA insert.
