# ISSAAC-seq — research notes

## Sources read

- 🟢 Xu et al. 2022 author bioRxiv v2 preprint and the published Supplementary Information (MOESM1). The published main article remains inaccessible for direct reading.
- 🟢 Supplementary Table 4 (`s41592-022-01601-4_41592_2022_1601_MOESM4_ESM.xlsx`), droplet- and plate-workflow oligo tables.
- 🟢 10x Chromium Next GEM Single Cell ATAC v1.1 user guide CG000209 Rev G, Appendix A1, for the gel-bead P5–N16–S5 scaffold. The exact revision used in the ISSAAC experiments is unverified.

## Shared reactions

🟢 `ME_S5` and `ME_S7` are the canonical S5–ME and S7–ME oligos. Each is annealed to `ME_bottom`, written `Phos-C*T*G*T*C*T*C*T*T*A*T*A*C*A*C*A*T*C*/iInvdT/`: a 5′ phosphate, phosphorothioate linkages and a 3′ inverted-dT block are therefore part of the reagent, not merely a generic Nextera end.

🟢 S5/S7 Tn5 first tagments accessible chromatin. In situ RT then uses a 5′-biotin capture oligo containing a partial TruSeq arm, 10 random bases (UMI), 32 thymidines and VN. The second transposome differs by branch: published Supplementary Protocol 1 (droplet; MOESM1 lines 863 and 1087) assembles and uses S5/S5, whereas Protocol 2 (plate; lines 1806 and 2047) uses S7/S7. The plate RNA fragment drawn here retains the capture end and gains S7 at its other end; internal S7/S7 products are not amplified by the asymmetric primer set.

## Plate branch drawn on the public page

🟢 The plate capture oligo is `biotin-CTACACGACGCTCTTCCGATCTNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN`. The sorted well is one cell, and its 8-nt i7 well barcode is therefore the cell identity. An 8-nt i5 index identifies the plate and distinguishes the ATAC and RNA libraries.

🟢 The plate ATAC library has ordinary Nextera geometry. The plate RNA library is hybrid: the capture-derived end is completed to a TruSeq Read 1 arm, while the tagmented end retains Nextera ME–S7 geometry. Consequently it uses TruSeq Read 1 and Index 2 sequencing primers together with Nextera Index 1 and Read 2 primers. These mixed primer assignments are derived from exact landing sites in the final construct.

🟢 Plate libraries are sequenced 150 + 8 + 8 + 150 cycles (R1 + I1 + I2 + R2).

🟢 Supplementary Table 4 also prints the full `N7001`–`N7384` i7 well-index family and the RNA/ATAC-specific plate i5 PCR primer families. The page now exposes one source-printed member of each family (`N7001`, `RNA_plate_S502`, `ATAC_plate_SS502`) as ordering examples; the other family members remain in the source workbook, not represented by a falsely universal selected index.

## Droplet branch

🟢 The droplet workflow uses `biotin-CAGACGTGTGCTCTTCCGATCT` + 10-nt UMI + dT32VN and its second tagmentation uses S5/S5 Tn5, unlike the plate S7/S7 step. It sends tagged material through a 10x Single Cell ATAC v1.1 gel-bead reaction. CG000209 Rev G Appendix A1 prints the gel-bead P5–N16–S5 scaffold, placing the 16-nt cell barcode next to the Nextera S5 handle. The RNA branch is completed from its printed RNA-specific P7/index primer family. The public page names this branch but uses the plate workflow for its exact reaction and final-library drawings so the two architectures are not collapsed into a fictitious common construct. The revision actually used by Xu et al. is not established. 🟡
