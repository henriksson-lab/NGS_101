# ISSAAC-seq — research notes

## Sources read

- 🟢 Xu et al. 2022, DOI `10.1038/s41592-022-01601-4`, article and complete protocol supplement.
- 🟢 Supplementary Table 4 (`s41592-022-01601-4_41592_2022_1601_MOESM4_ESM.xlsx`), droplet- and plate-workflow oligo tables.
- 🟢 Teichlab ISSAAC-seq page, checked as a secondary reconstruction against the source table and protocol.

## Shared reactions

🟢 `ME_S5` and `ME_S7` are the canonical S5–ME and S7–ME oligos. Each is annealed to `ME_bottom`, written `Phos-C*T*G*T*C*T*C*T*T*A*T*A*C*A*C*A*T*C*/iInvdT/`: a 5′ phosphate, phosphorothioate linkages and a 3′ inverted-dT block are therefore part of the reagent, not merely a generic Nextera end.

🟢 S5/S7 Tn5 first tagments accessible chromatin. In situ RT then uses a 5′-biotin capture oligo containing a partial TruSeq arm, 10 random bases (UMI), 32 thymidines and VN. A second, S7/S7 transposome tagments the RNA/cDNA hybrid. The productive RNA fragment retains the capture end and gains S7 at its other end; internal S7/S7 products are not amplifiable by the asymmetric primer set.

## Plate branch drawn on the public page

🟢 The plate capture oligo is `biotin-CTACACGACGCTCTTCCGATCTNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN`. The sorted well is one cell, and its 8-nt i7 well barcode is therefore the cell identity. An 8-nt i5 index identifies the plate and distinguishes the ATAC and RNA libraries.

🟢 The plate ATAC library has ordinary Nextera geometry. The plate RNA library is hybrid: the capture-derived end is completed to a TruSeq Read 1 arm, while the tagmented end retains Nextera ME–S7 geometry. Consequently it uses TruSeq Read 1 and Index 2 sequencing primers together with Nextera Index 1 and Read 2 primers. These mixed primer assignments are derived from exact landing sites in the final construct.

🟢 Plate libraries are sequenced 150 + 8 + 8 + 150 cycles (R1 + I1 + I2 + R2).

## Droplet branch

🟢 The droplet workflow uses `biotin-CAGACGTGTGCTCTTCCGATCT` + 10-nt UMI + dT32VN and sends the tagged material through a 10x Single Cell ATAC gel-bead reaction. Its 16-nt cell barcode is carried next to the Nextera S5 handle. The RNA branch is completed from its printed RNA-specific P7/index primer family. The public page names this branch but uses the plate workflow for its exact reaction and final-library drawings so the two architectures are not collapsed into a fictitious common construct.
