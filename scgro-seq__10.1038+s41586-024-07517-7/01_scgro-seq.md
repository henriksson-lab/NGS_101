# scGRO-seq — research notes

## Sources read

- 🟢 Mahat et al. 2024, DOI `10.1038/s41586-024-07517-7`, PMC article, Methods and Supplementary Table 7 (`41586_2024_7517_MOESM1_ESM.pdf`).
- 🟢 Round 2026-10-17-breadth-first-16 re-read the published Methods in `s41586-024-07517-7_PMC11222150.xml.txt`: propargyl nucleotides and adaptor design lines 185–195; cell line/culture/sample preparation lines 197–215; nuclear run-on line 217; single-cell sorting line 219; CuAAC lines 221–227; RT and pre-amplification lines 229–233; DASH lines 235–237; final PCR and NGS lines 241–243; read trimming and demultiplexing line 247. Supplementary Information MOESM1 lines 131–150 (Table 7) and 167–240 (Tables 8 and 9).

## Construct path

🟢 Engaged polymerases extend nascent RNA during a nuclear run-on with 3′-O-propargyl NTPs, which both terminate extension and leave a clickable alkyne. 🟢 An azide-bearing oligo is joined by click chemistry and supplies the recovery handle plus molecular/cellular identifiers. Reverse transcription reads across the triazole junction before indexed Illumina amplification.

## Published azide adaptors

🟢 Supplementary Table 7 prints three alternative 5′-AzScBc designs. Here `J` denotes a well-specific cell-barcode base and `N` a random UMI base. The reported 43/86-nt lengths include the modified azide-dT at the 5′ end; the paper does not enumerate all 96 actual `J` sequences per design.

- Design 1, 10-nt cell barcode and 6-nt UMI: `5′-[azideT]JJJJJJJJJJNNNNNNGATCGTCGGACTGTAGAACTCTGAAC-3′`.
- Design 2, 12-nt cell barcode and 4-nt UMI: `5′-[azideT]NNNNJJJJJJJJJJJJGATCGTCGGACTGTAGAACTCTGCGGCCGTGCTCGTTTTCGAGCACGGCCGCAGAGTTCTACAGTCCGA-3′`. This 86-nt oligo folds into an EagI-site hairpin; its 3′ end primes reverse transcription, and EagI cuts the site before PCR.
- Design 3, 10-nt cell barcode and 5-nt UMI: `5′-[azideT]ACAGGNNNNNJJJJJJJJJJAGATCGGAAGAGCGTCGTGTAG[SpC3]-3′`. The terminal C3 spacer blocks unwanted extension. The 5′ `ACAGG` supports DASH depletion of empty-adaptor amplicons.

## Reaction conditions read in this pass

🟢 Nuclei are isolated by douncing in Triton X-100 buffer or cells are permeabilized in Tween-20/NP-40 buffer, then stored. The run-on uses 2× buffer with 400 mM KCl, 0.05% Sarkosyl and 50 µM of each 3′-(O-propargyl)-ATP/CTP/GTP/UTP (Jena NU-945 to NU-948) for 5 min at 37 °C, after which the nuclei are washed three times. 🟢 Single cells or nuclei are FACS-sorted into 96-well plates holding 2.5 µl of 8 M urea, one barcoded 5′-AzScBc oligo per well. 🟢 CuAAC runs at 30 nM adaptor, 800 mM sodium ascorbate, 15% PEG 8000, 1 mM CuSO4, 5 mM BTTAA and 2.66 M urea in 7.5 µl for 2 h at 50 °C, is quenched with 5 mM EDTA, pooled, freed of PEG with TRIzol and of the copper reagents through a 3 kDa Amicon filter, then fragmented with 10 mM ZnCl2 for 5 min at 65 °C. 🟢 RT uses RNase H− M-MuLV reverse transcriptase (Thermo EP0751) with 16% PEG 8000, 1 µM RT primer (omitted for the hairpin design) and 1 µM TSO for 2 h at 50 °C.

🟢 The cDNA is then size-selected on 10% denaturing PAGE away from the unclicked 5′-AzScBc DNA and the empty cDNA formed between adaptor and TSO — a physical separation that precedes DASH; the graph now draws it, and DASH remains a cleavage that blocks amplification rather than a physical removal. 🟢 Pre-amplification and the final PCR both use NEBNext Ultra II Q5 2× master mix (NEB M0544) with 0.5 µM unique dual-index primers.

🟡 Supplementary Tables 8 and 9 print the same six-cycle 98 °C / 60 °C / 72 °C program although the Methods call the final one a "two-step" program; the model therefore states six cycles only for the pre-amplification and makes no claim about the final cycle count. This is shown on the page.

🟢 For the two 43-nt designs, a separate RT primer is added; the 86-nt hairpin is itself the RT primer. The Methods specify a template-switching oligo, dual-index PCR primers, and a 64/43/8/8-cycle NovaSeq run, but do not print those ordering sequences. After the two-hour click reaction, the authors fragment the clicked RNA with ZnCl2, reverse-transcribe through the junction, size-select the cDNA on 10% denaturing PAGE, pre-amplify six cycles, cleave empty-adaptor amplicons with Cas9/guide RNA (DASH), then perform final indexed PCR (article Methods, lines 222–245). The undesired amplicons arise from RT of unclicked azide adaptor with TSO; they are an explicit discarded branch in the graph. The Methods print Read 1 trim motif `CCCCTGTCTCTTATACACAT` and Read 2 trim motif `AGATCGGAAGAGCGTCGTGT`, and extract cell barcodes and UMIs from the start of Read 1 (line 247).

🟡 The schematic uses design 3 as a representative 10/5 barcode/UMI construct; its final Illumina junction is inferred. 🔴 The exact 96-member cell-barcode inventories and the TSO, RT, PCR, and DASH guide ordering sequences remain unresolved from the source material reviewed here. The original DASH paper and its linked Table S1 explain Cas9 depletion but list guides for unrelated targets, so those guide sequences cannot be imported. The cited Miller–Kool 2004 azide-synthesis paper is purchase-gated; the defining article itself reports the iodinated-to-azide substitution conditions.
