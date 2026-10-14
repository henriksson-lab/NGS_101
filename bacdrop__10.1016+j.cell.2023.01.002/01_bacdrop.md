# BacDrop — research notes

Defining source: Ma et al., 2023, [doi:10.1016/j.cell.2023.01.002](https://doi.org/10.1016/j.cell.2023.01.002).

## Sources read

- 🟢 PMC full text, STAR Methods and Key Resources oligo table: fixation/permeabilization, RNase-H rRNA depletion, 384-well RT indexing, 10x droplet indexing, TdT poly(A) tailing and second-strand synthesis.
- 🔴 Supplementary Table S1 (the 384 RT primer sequences) is listed by the article but its PMC download is gated in the current source fetch. The article's other numbered supplements are likewise manual-retrieval entries.

## Construct

🟢 CB1 is introduced on the random RT primer and CB2 by the 10x bead in droplets; CB1+CB2 identifies a cell. The RT primer also carries a UMI.

🟢 The paper's sequencing specification allocates 16 cycles to Index 2 (i5), consistent with the 16-nt 10x droplet barcode, while the 13-nt plate barcode and 8-nt UMI are on the Read-1 side of the cDNA. The prior schematic compressed CB1 to eight bases and put CB2 in the insert; both positions have now been corrected.

🟢 TdT poly(A)-tails bacterial cDNA because bacterial transcript ends do not support the mammalian template-switching route efficiently.

🟡 Proprietary 10x bead bases remain placeholders.

## Source-printed structural oligos

🟢 The article's Key Resources table prints `SMRT_dT` (second-strand primer after TdT tailing) as `AAGCAGTGGTATCAACGCAGAGTTTTTTTTTTTTTTTVN`, `SMRT_PCR` as `AAGCAGTGGTATCAACGCAGAGT`, and the partial P5 primer as `AATGATACGGCGACCACCGAGA` (spaces removed from the printed table). `SMRT_dT` provides the synthetic priming site later read by `SMRT_PCR`; these are not 10x bead sequences. The final Nextera XT library PCR pairs P5 with the Illumina N7** Index 1 primer. The paper specifies Read 1 60 bp, Read 2 39 bp, i7 8 bp and i5 16 bp.

🔴 The 384 well-specific RT primer sequences are assigned to inaccessible Table S1. The main paper supplies their structure—partial Illumina primer binding site, 8-nt UMI, 13-nt CB1 and 6-nt random 3′ primer—but not every barcode member. The commercial 10x bead barcode (CB2) member sequences are also undisclosed in the accessible paper.
