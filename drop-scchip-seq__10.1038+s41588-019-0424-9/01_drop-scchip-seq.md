# Drop-scChIP-seq — research notes

## Sources read

- 🟢 Grosselin et al. 2019, DOI `10.1038/s41588-019-0424-9`, Supplementary Figures 3–4 and Supplementary Tables 4–5.

## Construct path

🟢 MNase fragments chromatin, droplet barcode beads assign cellular identity, and immunoprecipitation enriches the chosen histone-mark fraction. The target is defined by the antibody used for a run rather than an inline target barcode, unlike Multi-CUT&Tag.

🟢 Supplementary Figure 3 shows that a bead barcode is assembled by three split-pool ligations. Each round has 96 double-stranded index reagents and 4-nt 5′ cohesive ends. Supplementary Table 4 gives both strands of all three plates. The schematic uses the A01 representative from each plate; the construction code derives rather than types their overhangs:

- round 1: `TTCG` incoming end, 20-bp paired identity, `GTCA` on the other strand;
- round 2: `TGAC` incoming end, 20-bp paired identity, `GGGA` on the other strand;
- round 3: `TCCC` incoming end, 20-bp paired identity, `GTTG` on the other strand.

🟡 Adjacent ends are reverse complements (`CGAA`/`TTCG`, `GTCA`/`TGAC`, `GGGA`/`TCCC`, `GTTG`/`CAAC`). The shared cohesive-duplex constructor rejects the assembly if that stops being true.

🟢 SEQ1 supplies the photocleavable bead flank, half PacI site, T7 promoter and Read 2 primer site. Its bottom strand ends in a C3 spacer. SEQ2 supplies the common terminal flank after round 3. Incomplete products are capped by SEQ3 oligos carrying a 5′ C3 spacer and 3′ inverted dT.

🟢 Supplementary Figure 4 defines the amplified product as P7 — Read 2 primer — 76-bp cell-barcode region — 12-bp linker — 150–300-bp nucleosomal DNA — Read 1 primer — P5. It assigns 100 cycles to the barcode-side read and 50 cycles to genomic Read 1.

🟢 Enriched barcoded nucleosomes undergo T7 in-vitro transcription. The RNA is reverse-transcribed with SEQ4, whose terminal six random bases are a priming tract, not a UMI. PCR uses the exact P5 primer and an indexed P7/Read-2 primer.

🟡 From the final construct and standard primer binding, Read 2 covers the three index oligos plus the terminal `CAAC` boundary in cycles 1–76, the 12-bp linker in cycles 77–88, and reaches genomic DNA at cycle 89. Read 1 starts directly in nucleosomal DNA. Index 1 reads the PCR sample index.

## Remaining limitation

🔴 The supplied figure fixes the final barcode and linker dimensions but does not expose, base by base, which residues around the circularized PacI half-sites remain after processing. The page therefore draws the source-resolved 12-bp linker and states this limitation instead of inventing a complete PacI junction.
