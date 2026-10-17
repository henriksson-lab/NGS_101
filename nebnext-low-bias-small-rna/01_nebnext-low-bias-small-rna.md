# NEBNext Low-bias Small RNA

## 1. What it is

A gel-free small-RNA library protocol that uses randomized splint ligation to reduce sequence-dependent adaptor ligation bias.

## 2. Sources and evidence

- 🟢 [NEB E3420 manual v2.0 (May 2026)](https://www.neb.com/en/-/media/nebus/files/manuals/manuale3420.pdf), read from `/husky/henriksson/NGS101_papers/sources/nebnext-low-bias-small-rna/manuale3420.pdf.txt`, and the [product page](https://www.neb.com/en-us/products/e3420-nebnext-low-bias-small-rna-library-prep-kit).
- 🟢 Captured RNAs are under 120 nt and require a 5′ monophosphate and 3′ hydroxyl.
- 🟢 The order is 3′-adaptor ligation, excess-adaptor removal, simultaneous 5′ ligation plus modification of the 3′ adaptor into an RT substrate, RT, bead selection and UDI PCR.
- 🟢 NEB's [E3420 RT-primer FAQ](https://www.neb.com/en-us/faqs/when-is-the-rt-primer-annealed-during-the-nebnext-low-bias-small-rna-library-prep-workflow-neb-e3420) specifies that the bottom strand of the 3′ adaptor becomes the RT primer during 5′ ligation. The kit does not add a separate RT-primer oligo.
- 🟢 The [NEB E3420 product page](https://www.neb.com/en/products/e3420-nebnext-low-bias-small-rna-library-prep-kit) cites [Maguire et al. (2020)](https://doi.org/10.1093/nar/gkaa480) as the research predecessor and explicitly says that the commercial kit has subsequent protocol enhancements. Maguire's Supplementary Table S2 prints prototype oligos 4–7: phosphorylated DNA 3′ adapter, uracil-bearing six-base-randomized 3′ splint with a 3′ inverted dT, RNA 5′ adapter, and a 5′ splint. The article describes preadenylation and splint annealing before 3′ then 5′ ligation. These are source-positive facts for the research method; no checked source establishes that E3420 uses identical oligos or modifications.
- 🔴 The checked E3420 manual names LB 3′ and 5′ adaptors and NEBNext LV Unique Dual Index Primer Pairs, but does not print complete kit adaptor or randomized-splint ordering strands and modifications. The UDI primer pairs are purchased separately as sets E3390/E3392/E3400/E3402/E3404/E3406/E3408; no particular pair belongs to the generic workflow. NEB's patent predecessors and the remaining product FAQs are still under review, so kit-specific sequence absence has not been accepted.

- 🟢 Maguire Supplementary Table S2 (supplementary_material.pdf.txt lines 247–270) prints the prototype 3′ adapter `/5Phos/AGATCGGAAGAGCACACGTCT`, 3′ splint `AGACGTGTGCTCTTCCGATCUNNNNNN /3InvdT/` and RNA 5′ adapter `rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC`. 🟡 That predecessor 5′ adapter would give a TruSeq Small RNA Read 2 read-through, not the standard TruSeq Read 2 trim motif printed for E3420 (below), so the prototype 5′ adapter is not the kit adapter.

## 3. Molecular path (E3420 manual v2.0_5/26, steps 1–8)

The protocol directly ligates RNA ends. It therefore does not inherit the dA-tailed DNA junction used by Ultra II DNA workflows.

1. 🟢 Denature 8 µl RNA 2 min at 70 °C, then ice (step 1).
2. 🟢 3′ ligation: LB 3′ Adaptor (diluted 1:10 or 1:30 for low input), LB Ligation Reaction Buffer and LB 3′ Ligase, 30 min at 25 °C (step 2).
3. 🟢 Excess 3′ adaptor removal: LB Adaptor Removal Enzyme Mix, 15 min at 30 °C then 5 min at 75 °C (step 3). The enzyme identity is not disclosed.
4. 🟢 5′ ligation: LB 5′ Adaptor, 10 mM ATP, LB Ligation Enhancer and LB 5′ Ligase, 30 min at 37 °C (step 4). The Figure 1 legend says the ligated 3′ adaptor is enzymatically modified into the RT primer in this same reaction.
5. 🟢 RT: LB RT Reaction Buffer, DTT, LB RT Enzyme Mix; 15 min at 42 °C, 10 min at 70 °C (step 5).
6. 🟢 Bead size selection of the RT reaction (step 6): 1X beads, discard beads carrying large fragments; supernatant plus 1X beads plus 40 µl isopropanol binds cDNA, supernatant discarded. The Figure 1 legend says this removes non-ligated adaptor and larger cDNA. (Model correction this round: the page previously placed the miRNA-versus-broad choice here.)
7. 🟢 PCR with a separately purchased NEBNext LV UDI Primer Pair and NEBNext Ultra II Q5 Master Mix, 98 °C/62 °C/72 °C, 9–22 cycles by input (step 7, Tables 7.1–7.2).
8. 🟢 Post-PCR choice (step 8): 8A 1.2X/0.8X double-sided bead size selection for inserts < 50 nt (miRNA, piRNA, tRF) or 8B 2X cleanup for small RNAs < 120 nt.

## 4. Final library and reads

🟢 The manual (Sequencing and Data Processing) says libraries are trimmed like standard TruSeq libraries (not TruSeq Small RNA): Read 1 `{{= low_bias_small_rna.MANUAL_TRIM["Read 1"] }}`, Read 2 `{{= low_bias_small_rna.MANUAL_TRIM["Read 2"] }}`. 🟡 These motifs equal `A` + the TruSeq Index 1 arm and the reverse complement of the TruSeq Read 1 arm; `final_library()` places the insert directly between those arms and raises if its read-through differs from the printed motifs. The earlier model's undisclosed 8-nt "kit remnant" placeholders on each side of the insert were removed because standard trimming would leave them in the reads. The adaptor and primer oligos themselves are still not printed, so the arms remain inferred.

🟢 A single 56-base read is typically sufficient; longer and paired-end reads can be used. 🟢 The LV UDI manual (v2.0_6/26, Overview) states each well holds premixed 8-base i5 and i7 index primers; the well used for a given library is user-chosen.

## 5. Open questions

- 🔴 Adaptor, splint and UDI primer strands and modifications; adaptor-removal and 3′-adaptor-modification enzymes; ligation geometry.
- Remaining NEB FAQs, supplemental-data page (HTTP 403) and patent predecessors not read in the 2026-10-17 bounded pass.
