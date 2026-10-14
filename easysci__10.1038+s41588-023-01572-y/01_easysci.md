# EasySci-RNA

**Defining source:** doi:10.1038/s41588-023-01572-y.

## Molecular boundary

🟢 Nuclei receive barcoded oligo(dT) and random-hexamer primers in RT wells, a ligated second barcode after redistribution, and a final indexed-PCR barcode after tagmentation.

## Evidence boundary

🟢 The authors' pipeline confirms separate P5 barcode and P7 index reads and merges matched random-hexamer and oligo(dT) products.

## Oligo source audit

🟢 The published supplementary workbook `s41588-023-01572-y_supp_41588_2023_1572_MOESM3_ESM.xlsx` contains explicit EasySci-RNA oligo lists: Supplementary Table 16 (short-dT RT), Table 17 (random-N RT), Table 18 (ligation) and Table 19 (P7 PCR/index). The detailed protocol in MOESM1 explicitly points to the separate Excel file for these sequences.

🟢 Both RT families carry a 5′ phosphate. For example, the A1 short-dT RT oligo is `/5Phos/ACGACGCTCTTCCGATCTNNNNNNNNTTCTCGCATGTTTTTTTTTTTTTTT`; the A1 random-N RT oligo is `/5Phos/ACGACGCTCTTCCGATCTNNNNNNNNCGGTCAAGAANNNNNN`. The A1 ligation primer is `AATGATACGGCGACCACCGAGATCTACACCCGCGGCTCAACACTCTTTCCCTAC`, and Table 19 gives the A1 P7 primer `CAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG`. Preserve the lowercase index text as printed; its table barcode column is the reverse complement.

🟡 These A1 rows are examples of larger plate-indexed sets. Their constant handles and chemical modifications are source-known, but a specific final library requires selecting matching well barcodes and the tagmentation/indexing configuration. The present role-level page can be made more precise without treating A1 as a universal reagent.
