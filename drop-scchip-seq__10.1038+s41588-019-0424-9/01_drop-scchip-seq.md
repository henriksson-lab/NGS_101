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

🟢 The original article Methods (local source lines 1948–1951) state that the first four Read 2 cycles were dark. Molecule coordinates below therefore do not imply that every depicted Read 2 base was called; the exact relation of these dark cycles to the 100 reported barcode-side cycles is not resolved by the diagram.

🟢 Enriched barcoded nucleosomes undergo T7 in-vitro transcription. The RNA is reverse-transcribed with SEQ4, whose terminal six random bases are a priming tract, not a UMI. PCR uses the exact P5 primer and an indexed P7/Read-2 primer.

🟡 From the final construct and standard primer binding, Read 2 points across the three index oligos plus the terminal `CAAC` boundary at molecular positions 1–76, the 12-bp linker at positions 77–88, and toward genomic DNA after position 88. These are construct coordinates, not base-called cycle positions because of the four initial dark cycles. Read 1 starts directly in nucleosomal DNA. Index 1 reads the PCR sample index.

🟢 The original GEO [GSE117309 barcode FASTA](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309) contains 884,736 80-base entries, the 96³ combinations described in article Methods lines 1964–1980. Each entry has three variable 20-mers (96 possibilities each) followed by the invariant 20-mer `ACAGATATCCACAGATATCC`. For the A01 example, the variable segments are the reverse complements of the Table 4 round 3, 2 and 1 paired cores in that order. Thus the file is a read-processing whitelist with a fixed flank, not 80 variable barcode bases; the physical synthesis oligos and modifications are in the publisher's Supplementary Table 4.

## Remaining limitation

🔴 The supplied figure fixes the final barcode and linker dimensions but does not expose, base by base, which residues around the circularized PacI half-sites remain after processing. The page therefore draws the source-resolved 12-bp linker and states this limitation instead of inventing a complete PacI junction.

🟢 The published article's Methods refer to SEQ1 and SEQ2 as “Supplementary Table 4” at local lines 1880–1890, and to SEQ4/SEQ5 the same way at lines 1943–1947. They cite “Supplementary Table 5” for the barcode structure at line 1969. In the original supplementary workbook, Table 4 (lines 4172–4749) holds the A/B/C index plates, whereas Table 5 (lines 4750–4762) holds SEQ1–5. The model follows the workbook for these oligos. The reason for the article's reversed table references is unknown.

## Bounded curation round 2 (2026-10-16)

Inspected: Original publisher paper PDF text1860–1900,1938–1954,2500–2511.

Added explicit RNA-digestion boundary between random RT and PCR, with warning because exact reagent/products are unassigned. Source confirmsT7 amplification,50-base genomic R1 and100-cycle barcode R2 with first4dark; existing table-number discrepancy retained.

Deferred: Original Supplementary Tables4/5 and PacI processing details; digestion reagent and exact products; biochemical binding and cleanup transitions. This capped pass does not establish a complete source review or publication-wide absence.
