# Slide-seq (v1)

**Defining paper:** Rodriques et al., *Science* (2019), doi:10.1126/science.aaw1219. The open author manuscript was read from PMCID PMC6927209, together with the bioRxiv supplementary methods and Table S1.

## Molecular boundary

🟢 Uniquely DNA-barcoded 10-µm beads are packed into a puck and decoded spatially by SOLiD sequencing-by-ligation. A fresh-frozen tissue section is transferred to that mapped surface.

🟢 Released mRNA is captured on bead oligos and converted into a Drop-seq-style 3′ RNA library. The bead barcode supplies position and the UMI identifies the captured molecule.

🟢 The supplement prints two bead-oligo designs. Sequence 2 contains the TruSeq handle, an 8-base split-pool barcode block, `TCTTCAGCGTTCCCGAGA`, a 7-base barcode block, an 8-base random UMI and oligo(dT)30. Sequence 1 additionally carries a photocleavable thymidine and T7-promoter region.

🟢 Table S1 prints the 58-nt TruSeq5 tagmentation primer, the 23-nt SMART PCR primer and the 22-nt TruSeq PCR-handle primer. Libraries were sequenced as 42 cycles Read 1, 8 cycles i7 and 50 cycles Read 2 on a NovaSeq S2 flow cell. The page therefore models the source-defined mixed endpoint: TruSeq Read 1 reads bead barcode and UMI, while Nextera Read 2 enters the transcript insert from the tagmented end.

## Evidence boundary

🟢 This page describes the original Slide-seq workflow, not Slide-seqV2, whose template-switching and second-strand changes have their own page.
