# ISSAAC-seq — research notes

## Sources read

- 🟢 Xu et al. 2022 published Supplementary Information (MOESM1), including Supplementary Protocols 1 (droplet) and 2 (plate). The published main article remains inaccessible for direct reading (manual: https://doi.org/10.1038/s41592-022-01601-4). An earlier pass also consulted the bioRxiv v2 preprint; the published supplement takes precedence.
- 🟢 Supplementary Table 4 (`s41592-022-01601-4_41592_2022_1601_MOESM4_ESM.xlsx`), droplet- and plate-workflow oligo tables.
- 🟢 10x Chromium Next GEM Single Cell ATAC v1.1 user guide CG000209 Rev G, Appendix A1, for the gel-bead P5–N16–S5 scaffold. The exact revision used in the ISSAAC experiments is unverified.

## Shared reactions

🟢 `ME_S5` and `ME_S7` are the canonical S5–ME and S7–ME oligos. Each is annealed to `ME_bottom`, written `Phos-C*T*G*T*C*T*C*T*T*A*T*A*C*A*C*A*T*C*/iInvdT/`: a 5′ phosphate, phosphorothioate linkages and a 3′ inverted-dT block are therefore part of the reagent, not merely a generic Nextera end.

🟢 S5/S7 Tn5 first tagments accessible chromatin. In situ RT then uses a 5′-biotin capture oligo containing a partial TruSeq arm, 10 random bases (UMI), 30 thymidines and VN. The second transposome differs by branch: published Supplementary Protocol 1 (droplet; MOESM1 lines 863 and 1087) assembles and uses S5/S5, whereas Protocol 2 (plate; lines 1806 and 2047) uses S7/S7. The plate RNA fragment drawn here retains the capture end and gains S7 at its other end; internal S7/S7 products are not amplified by the asymmetric primer set.

## Plate branch drawn on the public page

🟢 The plate capture oligo is `biotin-CTACACGACGCTCTTCCGATCTNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN`. The sorted well is one cell, and its 8-nt i7 well barcode is therefore the cell identity. An 8-nt i5 index identifies the plate and distinguishes the ATAC and RNA libraries.

🟢 The plate ATAC library has ordinary Nextera geometry. The plate RNA library is hybrid: the capture-derived end is completed to a TruSeq Read 1 arm, while the tagmented end retains Nextera ME–S7 geometry. Consequently it uses TruSeq Read 1 and Index 2 sequencing primers together with Nextera Index 1 and Read 2 primers. These mixed primer assignments are derived from exact landing sites in the final construct.

🟢 Plate order of operations (published MOESM1 Supplementary Protocol 2, lines 1794–2386): Tn5-S5/S7 chromatin tagmentation (steps 19–24); in situ RT with TruseqR1_oligo_dT (25–30); Tn5-S7/S7 hybrid tagmentation (31–36); one in-nucleus "EXO I digestion and gap fill-in" with Maxima H Minus RT, dNTPs and Exonuclease I at 37 °C for 15 min (37–39); DAPI sort of one nucleus per well into a lysis plate pre-loaded with `Truseq_S5_short` and one N7xxx i7 primer (7–10, 43–45; P24 is a no-cell control); 65 °C SDS lysis and Tween-20 quench (46–47); a 72 °C 5-min gap fill-in followed by 4 RNA pre-PCR cycles with Q5 (48–49); addition of `Nextra_S5_short` and 6 ATAC+RNA pre-PCR cycles (51–52); pooling one plate, Zymo purification (53–56); Exonuclease I digestion of residual primers, 80 °C inactivation and 1.2× beads (57–60); one 8-cycle final PCR containing `ATAC_plate_SS5xx`, `RNA_plate_S5xx` and `illumina P7` (`CAAGCAGAAGACGGCATACGAGAT`) (61–62); 1.0× beads (63).

🟡 Model correction (2026-10-17): the earlier graph drew a separate gap fill, an Exonuclease I cleanup, a physical ATAC/RNA split and one dual-index PCR per branch. The graph now follows the order above. The i7 well barcode enters in the well and i5 enters in the pooled final PCR. That final PCR is one tube, drawn as two primer-pair readings so that each product keeps exact geometry. `RNA_plate_S5xx` ends in `…CTACACGACGCTCTTCCG`, so its 3′ 18 nt anneal within the 22-nt capture arm and the primer extends it to the full TruSeq Read 1 site.

🟢 Plate libraries are sequenced 150 + 8 + 8 + 150 cycles (R1 + I1 + I2 + R2) on a NovaSeq 6000 (step 67, line 2386).

🟢 Supplementary Table 4 also prints the full `N7001`–`N7384` i7 well-index family and the RNA/ATAC-specific plate i5 PCR primer families. The page now exposes one source-printed member of each family (`N7001`, `RNA_plate_S502`, `ATAC_plate_SS502`) as ordering examples; the other family members remain in the source workbook, not represented by a falsely universal selected index.

🟢 Exact input states (2026-10-17): the graph's chromatin transposome is the two loaded adapters (ME_S5 and ME_S7, each on 5′-Phos ME_bottom) rather than drawn rows, and `Truseq_S5_short`, `Nextra_S5_short`, `illumina P7`, the `N7xxx` family (P7 + 8-nt i7 well-index window + S7) and the `ATAC_plate_SS5xx` / `RNA_plate_S5xx` families (P5 + 8-nt i5 window + S5 or `ACACTCTTTCCCTACACGACGCTCTTCCG`) carry the Supplementary Table 4 sequences, with the index windows as identifier placeholders.

## Droplet branch

🟢 The droplet workflow uses `biotin-CAGACGTGTGCTCTTCCGATCT` + 10-nt UMI + dT30VN and its second tagmentation uses S5/S5 Tn5, unlike the plate S7/S7 step. It sends tagged material through a 10x Single Cell ATAC v1.1 gel-bead reaction. CG000209 Rev G Appendix A1 prints the gel-bead P5–N16–S5 scaffold, placing the 16-nt cell barcode next to the Nextera S5 handle. The RNA branch is completed from its printed RNA-specific P7/index primer family. The public page names this branch but uses the plate workflow for its exact reaction and final-library drawings so the two architectures are not collapsed into a fictitious common construct. The revision actually used by Xu et al. is not established. 🟡
