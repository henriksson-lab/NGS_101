# SUM-seq — research notes

## Sources read

- 🟢 Llorens-Rico et al. 2025, DOI `10.1038/s41592-025-02700-8`, PMC article, Methods and Supplementary Table 1.
- 2026-10-17 bounded pass: PMC12165849 JATS text lines 224 (overview), 398–427 (Transposome generation, DNA transposition, In situ RT, cDNA/mRNA hybrid tagmentation, Gap fill and ExoI, GEM generation, Pre- and final library amplification, Sequencing); MOESM3 Table 1a–d (sheet lines 42–262 of the workbook text).

## Construct path

🟢 Step order (Methods lines 400–425): barcoded-Tn5 DNA transposition of each fixed-nuclei aliquot (30 °C, 30 min; 2× EDTA stop buffer; three washes) → in situ RT with the sample's barcoded oligo(dT) (Maxima H Minus, PEG 8000 in most experiments, 50 °C then three 8–50 °C cycles) → pooling → tagmentation of cDNA–mRNA hybrids with a Read1N-only Tn5 (37 °C, 30 min) → gap fill with Maxima H Minus plus Thermolabile Exonuclease I (37 °C, 15 min) → 10x scATAC GEM generation with ~500 nM blocking invdT oligo (omitted for macrophages), 72 °C 5 min and 4–12 linear barcoding cycles → recovery, MyOne Silane and 1.4× SPRI → three-primer pre-amplification (Partial-P5-PTO, TruseqR2, ATAC-spacer; 5 cycles; 1.2× SPRI) → split in two → modality PCR (Partial-P5-PTO with P7-i7-spacer for ATAC or P7_i7_TruSeqR2 for RNA; 5 cycles plus qPCR-determined cycles) → 1.2× (ATAC) or 0.72× (RNA) SPRI. The earlier graph split nuclei into separate RNA/ATAC branches before RT and omitted gap fill/ExoI, the blocker, pre-amplification and cleanups; it has been corrected.

🟢 Fixed nuclei receive sample identity twice by modality-specific routes: barcoded Tn5 for ATAC and barcoded oligo(dT) reverse transcription for RNA. RNA/cDNA hybrids are tagmented to add the primer site required for microfluidic indexing. 🟢 Sample-indexed nuclei are pooled and overloaded into 10x droplets, where a second droplet barcode is added. Cell identity is the sample-index × droplet-barcode combination; RNA and ATAC products are split only after droplet breakage and pre-amplification.

## Structural oligos in the original supplement

🟢 Supplementary Table 1b–d, in the publisher XLSX (`MOESM3`), prints the general reagents and the 96-member ATAC and RNA sample-index families. Its indices are **11 bases**, not the eight bases formerly drawn in the page. `SUMseq-ATAC-1` is `CATACGGCAGAAGACGAACACCTCCAACTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` (index `ACCTCCAACTA`). `SUMseq-RNA-1` is `CAGACGTGTGCTCTTCCGATCTNNNNNNNNNNAGTGATTAGCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN` (10-base UMI and index `AGTGATTAGCA`). These are examples from separate published index families, not an identical shared oligo.

🟢 The common Tn5 loading oligo is `[Phos]CTGTCTCTTATACACATCT[ddc]`. The hopping blocker is `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-[invdT]`. Table 1b also prints `Nextera Read 1` `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, phosphorothioate-ending Partial-P5 `AATGATACGGCGACCACCGAG*A` and Partial-P7 `CAAGCAGAAGACGGCATACGA*G`, `TruseqR2` `CAGACGTGTGCTCTTCCGATCT`, `ATAC-spacer` `CATACGGCAGAAGACGAAC`, and the two modality-specific P7-index PCR oligos displayed on the page. `[library_index]` in those PCR oligos is the source's placeholder, not a specified index sequence. The chemical modifications are part of the ordering sequences.

## Final libraries and reads

🟡 The final constructs are now assembled from the Table 1 oligos (derived): ATAC = P5 · droplet barcode (16) · Nextera Read 1 · insert · Read2N rc · sample index (11) · ATAC-spacer rc · library index · P7 rc; RNA = P5 · droplet barcode · Nextera Read 1 · cDNA (sense) · poly(A) · sample index (11) · UMI (10) · TruSeq Read 2 rc · library index · P7 rc. 🔴 The 10x scATAC gel-bead oligo (P5 remainder and the droplet barcode in front of Read 1N) is not printed in the article or Table 1; those regions are drawn inferred. 🔴 `[library_index]` has no printed length; an 8-nt inferred placeholder is used. PCR primer landing is not yet modelled as separate geometry.

🟢 Sequencing (Methods line 427): NovaSeq6000/NextSeq2000 — ATAC Read 1 55, Index 1 11, Index 2 16, Read 2 55; RNA Read 1 95, Index 1 6, Index 2 16, Read 2 21. Aviti (T cell, CRISPRi/a screen) — ATAC 79/11/16/79, RNA unchanged. 🟡 Derived from the assembled constructs: RNA Read 2 covers UMI then sample index (10 + 11 = 21 cycles); ATAC Index 1 from the Nextera Index 1 site reads the 11-nt sample index (Table 1a: "either this sequence or the reverse complement is used for sample demultiplexing"); Index 2 reads the droplet barcode. The ATAC library index lies 19 nt beyond the 11 acquired Index 1 cycles; how (or whether) it is read is not stated in the inspected Methods.

## Open questions

- 10x scATAC gel-bead oligo and kit revision (10x user guide not inspected in this pass).
- `[library_index]` sequences/length and ATAC library-index readout.
- Datlinger 2021 (ref. 6) two-stage indexing and Hennig 2018 (ref. 84) Tn5 preparation, and the author SUMseq repository (manual; git.embl.org), were not followed.
