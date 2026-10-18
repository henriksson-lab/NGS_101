# TruSeq DNA PCR-Free

## Sources read

🟢 Illumina, *TruSeq DNA PCR-Free Reference Guide*, document 1000000039279 v00 (October 2017; cover line 4, md5 816296f601b8492bdd87bd81d41a3ae5). Sections read in the extracted text: Introduction and DNA input (lines 112–140), Fragment DNA with the Covaris 350-bp and 550-bp tables (lines 604–700), Repair Ends and Select Library Size (lines 817–1126), Adenylate 3′ Ends (lines 1127–1232), Ligate Adapters and Clean Up Ligated Fragments (lines 1235–1470), Check Libraries (lines 1526–1660), Product Contents with the index-adapter catalogue numbers (lines 1890–1925) and Index Adapter Sequences (lines 2480–2483).

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, is the central source used by `lib/illumina.py` for TruSeq adapter and sequencing-primer sequences.

## Molecular workflow

1. Covaris shearing makes double-stranded fragments with 3′ or 5′ overhangs. The guide provides separate 350-bp and 550-bp insert workflows.
2. 🟢 End Repair Mix (ERP 2 or ERP 3, 30 °C for 30 minutes) blunts the ends: "A 3′ to 5′ exonuclease activity removes the 3′ overhangs. A 5′ to 3′ polymerase activity completes the 5′ overhangs" (line 818–819).
3. 🟢 Size selection is two-sided and the two rounds discard different material. "Remove Large DNA Fragments" adds 160 µl *diluted* SPB and transfers 250 µl **supernatant** onward, so the oversized fragments leave on the beads; "Remove Small DNA Fragments" then adds 30 µl *undiluted* SPB and discards the supernatant, so the undersized fragments leave in solution. The graph now declares both rounds with their own waste terminals rather than one generic off-size discard.
4. 🟢 A-tailing (CTA + ATL or ATL 2; 37 °C 30 min, 70 °C 5 min) adds one 3′ A "to prevent them from ligating to each other during adapter ligation reaction. One corresponding thymine (T) nucleotide on the 3′ end of the adapter provides a complementary overhang" (lines 1128–1131).
5. 🟢 Ligation adds CTL, LIG 2 and the separately purchased DNA adapters, 30 °C for 10 minutes, and is then stopped with 5 µl STL (Stop Ligation Buffer) before two rounds of SPB cleanup (lines 1399–1470). STL is a separate quenching step, followed by two sequential SPB binding/wash/elution operations. Their discarded fractions remain separate.
6. 🟢 The library is quantified and pooled directly. There is no enrichment PCR, so qPCR is required: "Non-qPCR methods quantify molecules that do not have adapters on both ends and do not form clusters" (lines 1537–1539).

## Final construct and readout

🟢 The ligated adaptor is full-length, which is what makes the workflow PCR-free: the guide's ligation step "ligates index adapters to the DNA fragment ends, preparing them for hybridization to a flow cell", and its ligation control is described as joining "3′-T overhang adapters to 3′-A overhang inserts". *Illumina Adapter Sequences* prints both arms of that adaptor — `Index 1 (i7) Adapters GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG` and `Index 2 (i5) Adapters AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT` (single-index kits instead use the one-oligo `TruSeq Universal Adapter`) — so P5, P7 and every sample index are present before ligation and no primer installs them afterwards.

🟢 The same document writes the dA/dT junction into the oligo listing itself: "the additional A base is in parentheses in the i7 adapter", `(A)GATCGGAAGAGC…`. The model therefore draws a ligation junction at both insert ends: the adaptor's own 3′ dT on the left, the insert's 3′ dA on the right.

🟡 The adapter stem is computed by antiparallel overlap of the two source arms. Both the reagent and the sealed product retain their forked platform tails. Each physical sealed strand contains a universal arm, a genomic strand and an indexed arm; Scene checks the paired stems, genomic insert and dA/dT junctions. Ligation does not synthesize the complements of the platform tails. Flow-cell copying is outside this preparation graph.

🟢 The guide lists separate single-index, combinatorial-dual-index and TruSeq unique-dual-index components (Appendix A, Product Contents, lines 1905–1925: IDT for Illumina–TruSeq DNA UD Indexes 20020590/20022370, TruSeq DNA Combinatorial Dual Indexes 20015949, TruSeq DNA Single Indexes Set A 20015960 and Set B 20015961; "Index adapter components must be purchased separately", line 138). The schematic uses the documented TruSeq DNA combinatorial-dual-index adapter family; it does not represent the other two purchasable families.

🟡 The purchased adapter and its well determine the index identities. The model records i5 and i7 as sample-index features without inventing a sample-specific choice. Physical eight-base index windows are assigned only to this selected v22 CD reference. The run recipe and primer geometry after flow-cell copying remain unassigned; physical index length does not assign acquired cycles.

## Open questions

- 🟡 The guide gives insert-size workflows only as Covaris settings and input amounts (350 bp from 1 µg, 550 bp from 2 µg, lines 123–130); the model keeps the insert as a placeholder rather than asserting a drawn length.
- 🟡 The inspected guide delegates adapter identity to *Illumina Adapter Sequences* 1000000002694 (lines 2480–2483). The index member and plate well therefore stay unassigned, and the schematic's combinatorial-dual family is representative, not a kit-fixed choice.
- 🟡 The 2017 v00 guide is the revision recorded here. Later TruSeq DNA PCR-Free revisions and the IDT-UD-only repackaging are not compared in this pass.

## Native migration, 2026-10-18 batch07

The local module owns its graph, source arm assembly, forked sealed product and recovery sequence. This pass inspected two originals: v00 title and end repair, two-sided size selection, dA-tailing, ligation, both cleanup rounds and qPCR sections (text817–840,1020–1152,1385–1540); v22 CD adapter arms and D701/D501-family tables (text21090–21230). The v00 end-repair control description (text2123–2168) explicitly assigns ERP2/ERP3 addition of the 5′ phosphates needed for ligation; the substrate 5P marks use that positive statement. The mandatory bundle oligo scraper was run and relevant arm hits checked against those original sections. Earlier broader reading remains historical scope.

The v22 adapter document is an explicitly selected later reference, not proof of a 2017 purchased lot. This pass does not compare kit revisions, select an index well, assign the executed run, or confirm missing information. No source-review status is upgraded. qPCR is analytical quantification, not a library amplification step.

The ligation graph also retains the guide-noted nonclustering molecular fraction lacking two productive adapter ends. Its component-specific recovery through SPB remains unresolved. Analytical qPCR does not physically discard that fraction.
