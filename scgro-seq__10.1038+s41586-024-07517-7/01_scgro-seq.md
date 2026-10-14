# scGRO-seq — research notes

## Sources read

- 🟢 Mahat et al. 2024, DOI `10.1038/s41586-024-07517-7`, PMC article, Methods and Supplementary Table 7 (`41586_2024_7517_MOESM1_ESM.pdf`).

## Construct path

🟢 Engaged polymerases extend nascent RNA during a nuclear run-on with 3′-O-propargyl NTPs, which both terminate extension and leave a clickable alkyne. 🟢 An azide-bearing oligo is joined by click chemistry and supplies the recovery handle plus molecular/cellular identifiers. Reverse transcription reads across the triazole junction before indexed Illumina amplification.

## Published azide adaptors

🟢 Supplementary Table 7 prints three alternative 5′-AzScBc designs. Here `J` denotes a well-specific cell-barcode base and `N` a random UMI base. The reported 43/86-nt lengths include the modified azide-dT at the 5′ end; the paper does not enumerate all 96 actual `J` sequences per design.

- Design 1, 10-nt cell barcode and 6-nt UMI: `5′-[azideT]JJJJJJJJJJNNNNNNGATCGTCGGACTGTAGAACTCTGAAC-3′`.
- Design 2, 12-nt cell barcode and 4-nt UMI: `5′-[azideT]NNNNJJJJJJJJJJJJGATCGTCGGACTGTAGAACTCTGCGGCCGTGCTCGTTTTCGAGCACGGCCGCAGAGTTCTACAGTCCGA-3′`. This 86-nt oligo folds into an EagI-site hairpin; its 3′ end primes reverse transcription, and EagI cuts the site before PCR.
- Design 3, 10-nt cell barcode and 5-nt UMI: `5′-[azideT]ACAGGNNNNNJJJJJJJJJJAGATCGGAAGAGCGTCGTGTAG[SpC3]-3′`. The terminal C3 spacer blocks unwanted extension. The 5′ `ACAGG` supports DASH depletion of empty-adaptor amplicons.

🟢 For the two 43-nt designs, a separate RT primer is added; the 86-nt hairpin is itself the RT primer. The Methods specify a template-switching oligo, dual-index PCR primers, and a 64/43/8/8-cycle NovaSeq run, but do not print those ordering sequences. 🟡 The schematic uses design 3 as a representative 10/5 barcode/UMI construct; its final Illumina junction is inferred. 🔴 The exact 96-member cell-barcode inventories and the TSO, RT, PCR, and DASH guide ordering sequences remain unresolved from the source material reviewed here.
