# 10x Visium Spatial Gene Expression

## 1. What it is

Visium captures polyadenylated RNA from a permeabilized tissue section on a glass array whose oligo-dT primers carry known spatial coordinates and UMIs.

## 2. Sources and evidence

- 🟢 10x Genomics, *Visium Spatial Gene Expression Reagent Kits User Guide*, CG000239 Rev F.
- 🟢 The guide defines imaging followed by permeabilization, on-slide reverse transcription, second-strand synthesis, cDNA denaturation and amplification, fragmentation, end repair/dA-tailing, adapter ligation and sample-index PCR.
- 🟢 The sequencing layout is Read 1 28 cycles, i7 10 cycles and Read 2 90 cycles; Read 1 contains a 16-nt spatial barcode and 12-nt UMI.
- 🔴 The complete capture-oligo sequence and barcode whitelist are commercial and are not assigned invented bases here.

## 3. Molecular path

The spatial coordinate is installed during reverse transcription, before cDNA leaves the slide. Subsequent library construction is bulk, but the inline barcode preserves the originating capture spot.
