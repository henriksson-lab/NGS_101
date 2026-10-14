# 10x Visium Spatial Gene Expression

## 1. What it is

Visium captures polyadenylated RNA from a permeabilized tissue section on a glass array whose oligo-dT primers carry known spatial coordinates and UMIs.

## 2. Sources and evidence

- 🟢 10x Genomics, *Visium Spatial Gene Expression Reagent Kits User Guide*, CG000239 Rev F.
- 🟢 The guide defines imaging followed by permeabilization, on-slide reverse transcription, second-strand synthesis, cDNA denaturation and amplification, fragmentation, end repair/dA-tailing, adapter ligation and sample-index PCR.
- 🟢 The sequencing layout is Read 1 28 cycles, i7 10 cycles and Read 2 90 cycles; Read 1 contains a 16-nt spatial barcode and 12-nt UMI.
- 🟢 The **Oligonucleotide Sequences appendix** (CG000239 Rev F, pp. 70–71) prints the slide capture primer as `5′-CTACACGACGCTCTTCCGATCT-N16-N12-(T)30VN-3′`, with spatial barcode `N16` and UMI `N12`. It prints the TSO `5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′`, preserving three terminal RNA guanosines.
- 🟢 The same appendix prints second-strand primer `AAGCAGTGGTATCAACGCAGAG`, cDNA PCR forward `CTACACGACGCTCTTCCGATCT`, and cDNA PCR reverse `AAGCAGTGGTATCAACGCAGAG` (all 5′→3′). The reverse primer is the same sequence as the second-strand primer.
- 🟢 The ligation adapter is shown as top `5′-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3′` and bottom `3′-TCTAGCCTTCTCG-5′`, leaving a 3′ T overhang on the bottom strand for the dA-tailed insert. The appendix also gives the fixed flanks of both 10-nt sample-index PCR primers; the variable index values depend on the chosen Dual Index TT well.
- 🟢 [10x Space Ranger spatial-output documentation](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs) states that the spot barcode inclusion list and slide-coordinate mapping are bundled with Space Ranger and written to `tissue_positions.csv`; it lists 4,992 barcodes for the 6.5-mm Visium array. The barcode set is therefore accessible as analysis reference data, although individual spot bases are not recopied into this schematic.
- 🔴 The guide does not identify the slide attachment chemistry or a single fixed spatial-barcode/UMI species; the appendix gives the complete *sequence pattern* and the random/spot-specific windows. The spatial-barcode whitelist and chosen sample-index well are separate from the structural oligo pattern.

## 3. Molecular path

The spatial coordinate is installed during reverse transcription, before cDNA leaves the slide. Subsequent library construction is bulk, but the inline barcode preserves the originating capture spot.
