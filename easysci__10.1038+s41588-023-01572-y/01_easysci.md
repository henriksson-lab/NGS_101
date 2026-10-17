# EasySci-RNA

**Defining source:** Sziraki et al., *Nature Genetics* (2023), doi:10.1038/s41588-023-01572-y (PMC10703679). The article Methods defer library construction to Supplementary Protocol 1 (publisher MOESM1, text lines 21–514); oligos are in MOESM3 Supplementary Tables 16–19. The bioRxiv preprint (10.1101/2022.09.28.509825) is superseded and not used here.

## Molecular boundary

🟢 Nuclei receive barcoded oligo(dT) and random-hexamer primers in RT wells, a ligated second barcode after redistribution, and a final indexed-PCR barcode after tagmentation.

## Step order (Supplementary Protocol 1, bounded pass 2026-10-17)

- 🟢 Steps 13–16 (text 249–262): each barcoded ligation primer is annealed 1:1 with the Ligation Adaptor Primer (95 °C 2 min, −1 °C/min to 20 °C) and diluted to 3.125 µM.
- 🟢 Step 17 (264–283): in-house Tn5 (Hennig et al. 2018) is loaded with a single adapter, Tn5-ME-B `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` annealed to `/5Phos/CTGTCTCTTATACACATCT` (Tn5-ME-rev).
- 🟢 Steps 30–33 (377–403): ≤20,000 nuclei per well of 2×96 plates; 1 µL each of 50 µM short-dT and randomN primer, 55 °C 5 min, ice; Maxima H Minus RT with a 4→50 °C ramp and 55 °C 15 min.
- 🟢 Steps 34–40 (406–425): pool, wash in NBB, redistribute into four 96-well plates, add 1 µL annealed ligation complex, T4 DNA ligase 30 min at room temperature, EDTA stop.
- 🟢 Steps 41–55 (426–469): pool, filter, count, 10,000 nuclei per well; NEBNext Ultra II Non-Directional second-strand module 16 °C 1 h; in-well DNA binding buffer + 0.8× AMPure.
- 🟢 Steps 56–62 (470–496): 1:100 Tn5 in tagmentation buffer, 55 °C 5 min; SDS/BSA with Universal P5 primer 55 °C 15 min; Tween-20; indexed P7 primer; NEBNext HF 2× mix, 72 °C 5 min, 98 °C 30 s, 12–15 cycles.
- 🟢 Steps 63–71 (497–515): pool, 0.8× AMPure, quantify; NovaSeq Read 1 100 bp, Read 2 100 bp, Index 1 10 bp, Index 2 10 bp.

🟡 Derived from the sequences, not drawn by the source: the Ligation Adaptor Primer `AGATCGGAAGAGCGTCGT|GTAGGGAAAGAGTGT` pairs with both the RT-primer head `ACGACGCTCTTCCGATCT` and the ligation-primer arm `ACACTCTTTCCCTAC`, so it bridges the ligation-primer 3′-OH to the /5Phos/ RT-primer end of the cDNA (splinted ligation; the model's `Scene` checks this pairing). The joined arm plus head is the TruSeq Read 1 sequence. The 72 °C 5 min PCR start is taken as Tn5 gap fill. Only fragments carrying the ligated P5 end and one Tn5-ME-B end can give P5/P7 products; ME-B/ME-B internal fragments are a candidate off-product with no described separation.

## Evidence boundary

🟢 The [authors' current analysis README](https://github.com/JunyueCaoLab/EasySci) labels demultiplexed files **R1 = genomic Read 1, R2 = P5 barcode, R3 = genomic Read 2, I1 = P7 barcode**. It demultiplexes on P7 alone, then retains the P5 barcode separately, and merges matched random-hexamer and oligo(dT) products. The README explicitly calls this an updated pipeline and links its [publication-era archive](https://zenodo.org/records/8395492). The original archive contents and [GEO GSE212606](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606) run metadata remain uninspected here. Supplementary Protocol 1 step 71 gives the run lengths (Read 1 100, Read 2 100, Index 1 10, Index 2 10 bp); which sequencing-primer reagents were used is not stated, so read starts remain derived.

## Oligo source audit

🟢 The published supplementary workbook `s41588-023-01572-y_supp_41588_2023_1572_MOESM3_ESM.xlsx` contains explicit EasySci-RNA oligo lists: Supplementary Table 16 (short-dT RT), Table 17 (random-N RT), Table 18 (ligation) and Table 19 (P7 PCR/index). The detailed protocol in MOESM1 explicitly points to the separate Excel file for these sequences.

🟢 Published discrepancy (shown on the page): Supplementary Protocol 1 materials (MOESM1 text 98–99) write the short-dT primer with 16 T after the barcode, whereas all 192 rows of Supplementary Table 16 have 15 T and the Fig. 1d legend says "a stretch of 15 thymine nucleotides". The model uses T15. The Fig. 1d legend also calls randomN oligos "composed of 8 random nucleotides" while Table 17 and the protocol give an N8 UMI and a 3′ N6; the model uses Table 17.

🟢 Both RT families carry a 5′ phosphate. For example, the A1 short-dT RT oligo is `/5Phos/ACGACGCTCTTCCGATCTNNNNNNNNTTCTCGCATGTTTTTTTTTTTTTTT`; the A1 random-N RT oligo is `/5Phos/ACGACGCTCTTCCGATCTNNNNNNNNCGGTCAAGAANNNNNN`. The A1 ligation primer is `AATGATACGGCGACCACCGAGATCTACACCCGCGGCTCAACACTCTTTCCCTAC`, and Table 19 gives the A1 P7 primer `CAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG`. Preserve the lowercase index text as printed; its table barcode column is the reverse complement.

🟡 These A1 rows are examples of larger plate-indexed sets. Their constant handles and chemical modifications are source-known, but a specific final library requires selecting matching well barcodes and the tagmentation/indexing configuration. The page now draws the constant handles with placeholder barcodes; A1 is not treated as a universal reagent.
