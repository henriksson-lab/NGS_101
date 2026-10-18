# TrAEL-seq — research notes

## Active native migration, 2026-10-18

Three originals inspected: published PMC8021198 Results69 and Methods199–225; original S2 detailed protocol `pbio.3000886.s009.doc.txt`20–289; NEB E7645S/L,E7103S/L manual **Version6.1_5/20**, title/header and overview79. The last explicitly describes end repair,5′phosphorylation and dA-tailing, qualifying the same named kit's activity. Actual experimental revision/lot remains unpinned. The TA-compatible drawing stays inferred; the paper's “blunting” terminology is not called a publication error. Older statements below that the vendor original was unread describe historical scope.

Native ordered adapter1 has5′[Phos], becomes5′rApp before capture, contains physical dU and BtndT; adapter2 has5′phosphate,31bp stem,4dU loop and terminal phosphorothioate-linked3′T. The graph preserves UMI8 through the Bst-connected original/hairpin/copy polymer, sonicated/end-prepared capture population and actual USER-opened PCR input. Full two-hairpin nick/seal/closure topology and USER footprints remain unassigned. Salmon-sperm carrier is an independent population added at on-bead EndPrep, with possible co-carried original/copy/carrier residuals retained through main PCR and both final cleanups. The test18cycle PCR is a QC aliquot; main cycle n is selected experimentally. E7335S is an example set only; no exact member/index width is assigned. Chemical tail1–4rA differs from computational removal of up to3T. MiSeq50SR/NextSeq50075SE are positive acquisition choices; index acquisition is unassigned.

The mandatory scraper returned27hits/4files; superseded preprint hits were excluded. Remaining deferred targets: actual kit revision/lot, library index set/member, complete native USER/PCR strand geometry and optional barcoded-adapter quantitative variant. This capped pass establishes no accepted absence claim and preserves historical reviews.

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
