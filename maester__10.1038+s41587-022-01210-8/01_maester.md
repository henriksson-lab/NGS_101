# MAESTER — research notes

Defining source: Miller *et al.*, “Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations”, *Nature Biotechnology* 40:1030–1034 (2022), [10.1038/s41587-022-01210-8](https://doi.org/10.1038/s41587-022-01210-8). The bioRxiv preprint (10.1101/2021.03.08.434450) is superseded and not used.

## Sources read (bounded pass 2026-10-17)

- 🟢 Published article, PMC9288977 author-manuscript JATS (`s41587-022-01210-8_PMC9288977.xml.txt`): Methods 10x library and sequencing (line 133), Seq-Well MAESTER (lines 139–153), MAESTER for 10x 3′ v3 (lines 157–171), MAESTER data processing (line 193).
- 🟢 Supplementary Table 3 (MOESM5 workbook, sheet “MAESTER 10x”): universal primers (lines 4–8), 12 mitochondrial primer mixes (lines 11–77), 24 GoT-P7-i7 and 24 GoT-P5-i5 sample-index primers and the footnotes (lines 78–125). Supplementary Table 2 (MOESM4, Seq-Well) was checked for its universal rows and first primer rows only (for example line 9, TyM61-MT-ND1_4009 = `/5Biosg/GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAACACCCTCACCACTACAATCT`, same targeting sequence as the 10x PvG1218 primer).
- 🟢 Supplementary Information PDF (MOESM1), Supplementary Figure 3, inspected as a rendered page because the text layer is font-encoded: the 10x final product, P7 orientation.
- Oligo scraper hits were checked against those tables.

## What is new

🟢 MAESTER enriches all 15 mitochondrial transcripts from the full-length cDNA of an existing 3′ scRNA-seq run (a separate fraction from the expression library), keeping cell barcode and UMI so that mtDNA variants can be assigned to expression-profiled cells. The Seq-Well/Drop-seq branch (Supplementary Table 2) uses 5′-biotinylated primers with a Nextera handle, `SMART-AC`, streptavidin capture and `CR1P`. The 10x branch uses different primers and no biotin enrichment (Methods line 157). The schematic shows the 10x branch only.

## Oligos (Supplementary Table 3, 5′→3′, no modifications printed)

- Mitochondrial primers: common handle `CACCCGAGAATTCCA` + transcript-specific sequence; for example PvG1218_MT-ND1_4009 = `CACCCGAGAATTCCAAACACCCTCACCACTACAATCT`. There are 65 primers in mixes 1–8 and R1–R4.
- GoT-P5-i5-BCXX = P5 + 8-nt i5 + `TCTTTCCCTACACGACGCTC` (PCR1).
- P5-generic = `AATGATACGGCGACCACCGAGATCTACAC`; GoT-P7-i7-BCXX = P7 + 8-nt i7 + `GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` (PCR2).
- 10x-Ci7P = `TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC`, 10x-Ci5P = `AGATCGGAAGAGCGTCGTGTAGGGAAAGA`: index spike-ins for the NovaSeq reverse-complement workflow only.
- PartialRead1 / PartialTSO: optional primers for full-length cDNA amplification (also used in the TCR-seq add-on).

🟡 GoT-P5-i5 anneals to the 10x Read 1 end through its 3′ 13 nt (overlap with PartialRead1; the model asserts it when the module is imported), rebuilding the Read 1 arm `TCTTTCCCTACACGACGCTCTTCCGATCT`.

## Step by step (10x, Methods lines 161–171)

1. 🟢 PCR1: 20 ng pooled 10x cDNA in 16 µl, 4 µl primer mix (mitochondrial primers 1 µM each, plus GoT-P5-i5-BCXX 10 µM, 10× stocks) and 20 µl KAPA HiFi HotStart ReadyMix. 95 °C 3 min, then 6 cycles of 98 °C 20 s / 65 °C 15 s / 72 °C 3 min, then 72 °C 5 min. There are 12 reactions, one per mix.
2. 🟢 Pool with empirically determined ratios (Supplementary Fig. 4D), purify with 1× AMPure XP and elute in 20 µl H2O.
3. 🟢 PCR2: 18 µl eluate, 2 µl of a 5 µM P5-generic + GoT-P7-i7 mix and 20 µl KAPA. 95 °C 3 min, then 6 cycles of 98 °C 20 s / 60 °C 30 s / 72 °C 3 min, then 72 °C 5 min. Purify with 0.8× AMPure XP and elute in 20 µl TE. The resulting libraries are 300–1,500 bp.
4. 🟢 NovaSeq SP 300-cycle kit: Read 1 28, Read 2 256, 2 × 8 index. No custom primers on the forward-strand workflow; 10x-Ci7P/10x-Ci5P on the reverse-complement workflow. i7 is demultiplexed with its reverse complement, and i5 according to the workflow (Table 3 footnotes).

## Final library and reads

🟢 Supplementary Figure 3 (read from the P7 end): P7–i7–`GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`–mitochondrial target…poly(A)–UMI–CB–`AGATCGGAAGAGCGTCGTGTAGGGAAAGA`–i5–P5. The figure labels the P7 arm “10x_R2”, but its sequence is the GoT-P7 small-RNA-type arm plus handle, not the TruSeq Read 2 arm of the parent 10x expression library. 🟢 Read 1 reads CB + UMI. 🟢 Read 2 starts at the target-primer sequence; the first 24 bp are trimmed before alignment (Methods line 193).

Model correction 2026-10-17: the previous page used a shared generic dual-index TruSeq shell, which put a TruSeq Read 2 arm where the GoT-P7 arm actually is. The library is now assembled from the printed primers.

## Discrepancy

🟢 The Methods (line 169) and Supplementary Figure 3 name the index primers XV-P7-i7-BCXX (and, in the figure, XV-P5-i5-BCXX). Supplementary Table 3 prints only the GoT-P7-i7 and GoT-P5-i5 series, and its PCR2 footnote uses GoT-P7-i7. The sequences drawn in the figure match the GoT series, so the model treats XV as an alias. The page shows this discrepancy.

## Open questions

- The Read 2 sequencing primer that anneals to the GoT-P7 arm is not printed. The Methods say only that the forward-strand workflow needs no custom primer.
- 🟡 The CB (16) and UMI (12) lengths are read from the Supplementary Figure 3 drawing and the 28-cycle Read 1. The poly(dT) length is not modelled.
- Deferred: 10x 3′ v3 user guide (parent cDNA ends), Nam/Landau 2019 GoT Table S5 (origin of P5-generic), Seq-Well branch chemistry beyond the Table 2 header.
