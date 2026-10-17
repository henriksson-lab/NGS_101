# TruSeq DNA PCR-Free

## Sources read

🟢 Illumina, *TruSeq DNA PCR-Free Reference Guide*, document 1000000039279 v00 (October 2017; cover line 4, md5 816296f601b8492bdd87bd81d41a3ae5). Sections read in the extracted text: Introduction and DNA input (lines 112–140), Fragment DNA with the Covaris 350-bp and 550-bp tables (lines 604–700), Repair Ends and Select Library Size (lines 817–1126), Adenylate 3′ Ends (lines 1127–1232), Ligate Adapters and Clean Up Ligated Fragments (lines 1235–1470), Check Libraries (lines 1526–1660), Product Contents with the index-adapter catalogue numbers (lines 1890–1925) and Index Adapter Sequences (lines 2480–2483).

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, is the central source used by `lib/illumina.py` for TruSeq adapter and sequencing-primer sequences.

## Molecular workflow

1. Covaris shearing makes double-stranded fragments with 3′ or 5′ overhangs. The guide provides separate 350-bp and 550-bp insert workflows.
2. 🟢 End Repair Mix (ERP 2 or ERP 3, 30 °C for 30 minutes) blunts the ends: "A 3′ to 5′ exonuclease activity removes the 3′ overhangs. A 5′ to 3′ polymerase activity completes the 5′ overhangs" (line 818–819).
3. 🟢 Size selection is two-sided and the two rounds discard different material. "Remove Large DNA Fragments" adds 160 µl *diluted* SPB and transfers 250 µl **supernatant** onward, so the oversized fragments leave on the beads; "Remove Small DNA Fragments" then adds 30 µl *undiluted* SPB and discards the supernatant, so the undersized fragments leave in solution. The graph now declares both rounds with their own waste terminals rather than one generic off-size discard.
4. 🟢 A-tailing (CTA + ATL or ATL 2; 37 °C 30 min, 70 °C 5 min) adds one 3′ A "to prevent them from ligating to each other during adapter ligation reaction. One corresponding thymine (T) nucleotide on the 3′ end of the adapter provides a complementary overhang" (lines 1128–1131).
5. 🟢 Ligation adds CTL, LIG 2 and the separately purchased DNA adapters, 30 °C for 10 minutes, and is then stopped with 5 µl STL (Stop Ligation Buffer) before two rounds of SPB cleanup (lines 1399–1470). STL and the two-round cleanup are declared in the graph; the cleanup's waste terminal carries the unbound adapters and reagents.
6. 🟢 The library is quantified and pooled directly. There is no enrichment PCR, so qPCR is required: "Non-qPCR methods quantify molecules that do not have adapters on both ends and do not form clusters" (lines 1537–1539).

## Final construct and readout

🟢 The ligated adaptor is full-length, which is what makes the workflow PCR-free: the guide's ligation step "ligates index adapters to the DNA fragment ends, preparing them for hybridization to a flow cell", and its ligation control is described as joining "3′-T overhang adapters to 3′-A overhang inserts". *Illumina Adapter Sequences* prints both arms of that adaptor — `Index 1 (i7) Adapters GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG` and `Index 2 (i5) Adapters AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT` (single-index kits instead use the one-oligo `TruSeq Universal Adapter`) — so P5, P7 and every sample index are present before ligation and no primer installs them afterwards.

🟢 The same document writes the dA/dT junction into the oligo listing itself: "the additional A base is in parentheses in the i7 adapter", `(A)GATCGGAAGAGC…`. The model therefore draws a ligation junction at both insert ends: the adaptor's own 3′ dT on the left, the insert's 3′ dA on the right.

🟡 Because the adaptor is forked, the sealed molecule's two platform arms are unpaired beyond the 12-bp TruSeq stem until the flow cell copies them. The reagent's fork geometry is modelled in `lib/batch_ngs.py` (`FullLengthAdapter.scene()`), while the library construct keeps the project's linear one-top-strand convention used for every other ligation product.

🟢 The guide lists separate single-index, combinatorial-dual-index and TruSeq unique-dual-index components (Appendix A, Product Contents, lines 1905–1925: IDT for Illumina–TruSeq DNA UD Indexes 20020590/20022370, TruSeq DNA Combinatorial Dual Indexes 20015949, TruSeq DNA Single Indexes Set A 20015960 and Set B 20015961; "Index adapter components must be purchased separately", line 138). The schematic uses the documented TruSeq DNA combinatorial-dual-index adapter family; it does not represent the other two purchasable families.

🟡 The purchased adapter and its well determine the index identities. The model records i5 and i7 as sample-index features without inventing a sample-specific choice. Standard TruSeq Read 1, Index 1, Index 2 and Read 2 primer landing sites are computed from the canonical construct.

## Open questions

- 🟡 The guide gives insert-size workflows only as Covaris settings and input amounts (350 bp from 1 µg, 550 bp from 2 µg, lines 123–130); the model keeps the insert as a placeholder rather than asserting a drawn length.
- 🔴 The guide never prints adapter bases; it delegates to *Illumina Adapter Sequences* 1000000002694 (lines 2480–2483). The index member and plate well therefore stay unassigned, and the schematic's combinatorial-dual family is representative, not a kit-fixed choice.
- 🟡 The 2017 v00 guide is the revision recorded here. Later TruSeq DNA PCR-Free revisions and the IDT-UD-only repackaging are not compared in this pass.
