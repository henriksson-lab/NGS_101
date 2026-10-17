# BacDrop — research notes

Defining source: Ma et al., 2023, [doi:10.1016/j.cell.2023.01.002](https://doi.org/10.1016/j.cell.2023.01.002).

## Sources read

- 🟢 PMC full text, STAR Methods and Key Resources oligo table: fixation/permeabilization, RNase-H rRNA depletion, 384-well RT indexing, 10x droplet indexing, TdT poly(A) tailing and second-strand synthesis.
- 🟢 10x Genomics, *Chromium Next GEM Single Cell ATAC Reagent Kits v1.1 User Guide*, CG000209 Rev G, Appendix p. 51: PN-2000210 gel-bead oligo constant arms and 16-base barcode position. This is the kit named in BacDrop STAR Methods, not a BacDrop-specific barcode list.
- 🔴 Supplementary Table S1 (the 384 RT primer sequences) is listed by the article but its PMC download is gated in the current source fetch. The article's other numbered supplements are likewise manual-retrieval entries.

## Construct

🟢 CB1 is introduced on the random RT primer and CB2 by the 10x bead in droplets; CB1+CB2 identifies a cell. The RT primer carries a partial primer binding site (PBS1), an 8-base UMI, a 13-base plate barcode and a 6-base random priming tract (Figure 1B legend, JATS line 827). The earlier panel's 8-base plate barcode and 8-base random tract were transcription errors and have been corrected.

🟢 Figure 1B image (PMC `nihms-1869349-f0002.jpg`, panel "For Sequencing") draws the final library as `P5 – CB2 – PBS1 – UMI – CB1 – cDNA – PBS2 – Index1 – P7`, with the Index 2 primer reading CB2, the Read 1 primer on PBS1 and the Read 2 primer on PBS2. **Model correction (2026-10-17):** the model had CB1 before the UMI; Read 1 now reads UMI cycles 1–8 then CB1 cycles 9–21 (derived by `seqprimers.feature_spans`).

🟡 The final library was previously drawn as a TruSeq ligation library with a dA junction. STAR Methods build it with Nextera XT tagmentation and P5/N7 PCR, and the bead primer ends in the Nextera Read 1 S5 arm `TCGTCGGCAGCGTC`. The model now uses a local Nextera endpoint: P5 and S5 from the bead oligo, an **inferred** remainder of the Read 1 site (`AGATGTGTATAAGAGACAG`, the rest of the standard Nextera Read 1 primer) on the RT primer PBS1 handle, and the Nextera XT mosaic-end/S7 tagmentation end (kit sequences, not printed by the article). Table S1 is needed to confirm PBS1.

🟢 The paper's sequencing specification allocates 16 cycles to Index 2 (i5), consistent with the 16-nt 10x droplet barcode, while the 13-nt plate barcode and 8-nt UMI are on the Read-1 side of the cDNA. The prior schematic compressed CB1 to eight bases and put CB2 in the insert; both positions have now been corrected.

🟢 Graph steps now follow STAR Methods lines 177–217: Tween-20/lysozyme permeabilization; NEBNext bacterial rRNA probes + thermostable RNase H; **DNase I gDNA removal** (previously missing); Maxima H Minus in-cell RT; TdT + dATP tailing; Chip H/ATAC gel-bead partitioning with KAPA HiFi and SMRT_dT; **SMRT_dT second-strand synthesis** (39 °C 5 min, 65 °C 10 min) as a step separate from the four bead-capture cycles; Recovery Agent/Dynabeads SILANE/0.6X AMPure cleanup; P5/SMRT_PCR enrichment; Nextera XT tagmentation; P5/N7 PCR; 0.6X AMPure.

🟢 After plate RT and pooling, TdT poly(A)-tails cDNA **before** droplet generation (STAR Methods, "In-cell cDNA 3′ poly-A tailing" then "Droplet generation"). The earlier panel reversed this order; it has been corrected.

🟢 The 10x guide prints the PN-2000210 bead-primer structure `AATGATACGGCGACCACCGAGATCTACAC-[16-nt barcode]-TCGTCGGCAGCGTC`. The barcode is a variable 16-base member; the version-matched 10x software whitelist has not yet been checked for this audit.

## Source-printed structural oligos

🟢 The article's Key Resources table prints `SMRT_dT` (second-strand primer after TdT tailing) as `AAGCAGTGGTATCAACGCAGAGTTTTTTTTTTTTTTTVN`, `SMRT_PCR` as `AAGCAGTGGTATCAACGCAGAGT`, and the partial P5 primer as `AATGATACGGCGACCACCGAGA` (spaces removed from the printed table). `SMRT_dT` provides the synthetic priming site later read by `SMRT_PCR`; these are not 10x bead sequences. The final Nextera XT library PCR pairs P5 with the Illumina N7** Index 1 primer. The paper specifies Read 1 60 bp, Read 2 39 bp, i7 8 bp and i5 16 bp.

🔴 The 384 well-specific RT primer sequences are assigned to inaccessible Table S1. The main paper supplies their structure—partial Illumina primer binding site, 8-nt UMI, 13-nt CB1 and 6-nt random 3′ primer—but not every barcode member. The 10x guide supplies the bead-primer constant arms and barcode length; its software whitelist still needs a version-matched source check. These are blocked information gaps pending a second source review and recovery of Table S1, not confirmed publication absences.

## Bounded pass 2026-10-17 (round breadth-first-07)

Read: PMC10014032 JATS STAR Methods lines 177–217, Figure 1 legend 825–827, Key Resources Table 1018–1083; Figure 1 image f0002.jpg (panel B); 10x CG000209 Rev G (now fetched) bead-oligo line 3766 and ATAC Read 1N/2N lines 3737–3781. Supplements 1–14 remain behind the PMC download gate (manual URLs in the source manifest); Table S1, the Nextera XT/N7 kit documents and the version-matched 10x barcode whitelist are deferred.
