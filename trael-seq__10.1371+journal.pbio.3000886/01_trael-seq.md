# TrAEL-seq — research notes

## Sources read

- 🟢 Kara et al. 2021, DOI `10.1371/journal.pbio.3000886`, PMC8021198 JATS (`journal.pbio.3000886_PMC8021198.xml.txt`): Results lines 62–71 (TdT/RNA-ligase rationale, Fig 1A legend, adapter roles) and Methods "TrAEL-seq library preparation and sequencing" lines 199–225; S2 File detailed protocol `pbio.3000886.s009.doc.txt` lines 20–289.
- 🟢 NEB E7335 Multiplex Oligos manual, version 8.0 (June 2024), primer table on page 4; the defining paper names E7335S for indexed library PCR.

## Construct path

### 2026-10-17 bounded pass: corrected reaction order

🟢 Methods line 217 / S2 File: half an agarose plug is equilibrated in TdT buffer, then tailed 2 h at 37 °C with Terminal Transferase (NEB M0315) and 0.4 mM **ATP**; Results line 62 says TdT adds 1–4 adenosines (1–3 nt on >95% of an in-vitro substrate), forming a substrate for RNA-ligase ligation. Adapter 1 is then ligated overnight at 25 °C with T4 RNA ligase 2 truncated KQ and PEG 8000, and plugs are washed extensively. The earlier model omitted this tailing step and joined adapter 1 directly to the native 3′-OH; that was a repository modelling error.

🟢 Line 219: plugs are melted and β-agarase-digested, DNA precipitated, extended with Bst 2.0 WarmStart (65 °C, 30 min), sonicated (Covaris E220, DF10/PIP175/200 cycles) and bound to MyOne C1 streptavidin beads, washed. Line 221: **on beads**, NEBNext Ultra II End Prep (with salmon-sperm carrier), adapter 2 ligation with Ultra II ligation mix/enhancer, washes that remove adapter-2 homodimers (Fig 1A legend), then two USER elutions. The earlier model put end repair and adapter 2 ligation before streptavidin capture; corrected.

🟢 Line 223: test PCR (18 cycles) to pick the cycle number, final PCR with NEBNext Ultra II Q5 master mix, NEBNext Universal and "any NEBNext Index" primer from "a NEBNext Multiplex Oligos set (e.g., NEB E7335S)"; 98 °C 30 s, n × (98 °C 10 s, 65 °C 75 s), 65 °C 5 min; two AMPure XP cleanups (40 µl, then 20 µl beads). E7335S is an example set, so its six-base i7 length is not a confirmed per-library assignment; the page draws the i7 window without a length.

🟡 Read structure (derived): Bst extends adapter 1's 3′ stem arm across `AGATCGGAAG`, the N8, the rA tail and the genomic strand, so the copied strand carries the full TruSeq Read 1 site followed by the UMI, a 1–4 nt T tract and genomic sequence complementary to the captured strand — consistent with Results line 71 ("the TrAEL-seq read is the reverse complement of the original DNA 3′ end"). Line 225: MiSeq 50 bp single read or NextSeq 500 75 bp single end, so Read 2 is not acquired; the index-read length is not printed.

🟡 Adapter 2 end: the Methods and Fig 1A legend speak of "blunting"/"polishing", but the published adapter 2 hairpin leaves one 3′ T (phosphorothioate-linked) outside its 31-bp stem, implying TA ligation to a dA-tailed end. The model draws that junction with the dA labelled inferred and shows this wording/structure tension on the page; the Ultra II End Prep vendor original was not read in this pass.

### Earlier summary

🟢 After TdT rA tailing (see above), a 5′-adenylated hairpin adapter is ligated to the tailed 3′ ends with truncated RNA ligase 2 KQ. It contains eight random bases, internal biotin and uracils. Bst polymerase fill/displacement stabilizes the junction. After sonication and second-adapter ligation, streptavidin capture and USER cleavage release end-derived molecules for PCR.

🟢 The defining paper prints adapter 1 in full as
`[Phos]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTUGCGCAGGCCATTGGCC[BtndT]GCGCUACACTCTTTCCCTACACGACGCT` and states that this synthesized oligo is subsequently 5′ adenylated. The schematic therefore shows the ligation-ready form as 5′-App and folds the complete 23-bp stem instead of abbreviating its loop. Adapter 2 is likewise printed in full, including four deoxyuridines and its terminal phosphorothioate bond.

🟢 The paper's named E7335S PCR set uses one six-base i7 index and an index-free Universal PCR primer, as confirmed by the NEB primer table. The earlier schematic incorrectly inherited two eight-base indexes from a generic TruSeq shell. The corrected library shows only Index 1, drawn without a length because the six-base E7335S set is only the paper's example. The paper reports 50-base MiSeq single reads or 75-base NextSeq 500 single reads, but does not name the selected E7335 index primer for each library. 🔴 The per-library index identity remains unverified while linked author documents and supplementary originals remain unread.
