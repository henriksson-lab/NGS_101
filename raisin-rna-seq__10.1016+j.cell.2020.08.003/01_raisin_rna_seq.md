# RAISIN RNA-seq — intact nuclei with ribosome-bound RNA

## Sources read

- 🟢 Drokhlyansky *et al.*, “The human and mouse enteric nervous system at single cell resolution,” *Cell* (2020), DOI `10.1016/j.cell.2020.08.003`, PMCID `PMC8358727` (defining paper and STAR Methods).
- 🟢 Shekhar *et al.*, *Cell* (2016), DOI `10.1016/j.cell.2016.07.054`, is the plate-library procedure explicitly cited by the defining paper's STAR Methods. Picelli *et al.*, SMART-seq2, DOI `10.1038/nmeth.2639`, is the earlier chemistry family; it is not the procedure directly cited at this point.

- Bounded pass 2026-10-17 (round breadth-first-13): published PMC8358727 JATS STAR Methods lines 1032–1039 (CST/TST composition), 1041–1054 (FACS into TCL/1% β-ME), 1056–1072 (whole-transcriptome amplification, library construction, sequencing), Key Resources 2915–2972 (Maxima H Minus EP0753, KAPA HiFi KK2602, RNAClean, AMPure XP, Nextera XT FC-131-1096, NextSeq 500 high output v2 75-cycle FC-404-2005) and 3030–3041 (Oligonucleotides); cited Shekhar *et al.* 2016 PMC5003425 JATS “Plate-based RNA-seq experiments” lines 288–294 (2.2× RNAClean, 21 cycles, 0.075 ng cDNA, ¼ Nextera XT volumes). Scraper hits at PMC8358727 lines 3225–3246 are unrelated Key-Resources template placeholders (siRNA, morpholino, GFP primer), and the Shekhar table at lines 276–284 holds Drop-seq oligos, not plate-library oligos. Figures S1–S7 and Tables S1–S6 remain PMC download-gated (unread).

## Molecular interpretation

- 🟢 RAISIN means **Ribosomes And Intact Single Nucleus**. Its CST extraction buffer contains 0.49% CHAPS, 146 mM NaCl, 1 mM CaCl2, 21 mM MgCl2 and 10 mM Tris pH 8.0; the isolated nucleus retains rough ER, attached ribosomes and ribosome-bound RNA.
- 🟢 In the plate implementation, one RAISIN is sorted into TCL with 1% β-mercaptoethanol. RNA is purified with RNAClean XP, reverse-transcribed from oligo-dT with Maxima reverse transcriptase (Key Resources: Maxima H Minus, EP0753) and a locked TSO, amplified for 21 cycles with KAPA HiFi HotStart ReadyMix, and the cDNA is purified twice with AMPure XP.
- 🟡 The cited Shekhar *et al.* plate procedure (the base the paper modifies) used 2.2× RNAClean SPRI beads without final elution, 0.075 ng cDNA per cell and ¼ of the standard Nextera XT volume for tagmentation and final PCR. The RAISIN paper does not restate these values, so they are inherited by citation, not confirmed for RAISIN libraries.
- 🟢 Amplified cDNA is converted with Nextera XT plus custom barcode adapters, pooled 384 wells at a time and sequenced on NextSeq 500 (Key Resources: high output kit v2, 75 cycles). 🔴 The acquired read configuration is not printed in the paper; no cycle counts are modelled.
- 🟢 The paper also implements droplet RAISIN RNA-seq on GemCode v2. The separate MIRACL-seq page shows that endpoint because deliberate droplet overloading is MIRACL-seq's defining change.
- 🟢 The paper's Key Resources Table prints the reverse-transcription primer as `5′-AAGCAGTGGTATCAACGCAGAGTACT(30)VN-3′`, the locked template-switching oligo as `5′-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3′`, and the PCR primer as `5′-AAGCAGTGGTATCAACGCAGAGT-3′` (PMC8358727 full-text XML, Oligonucleotides table). Preserve `rG` and `+G` as written: they denote the ribonucleotide and locked-nucleotide chemistry at the TSO 3′ end.
- 🔴 The custom Nextera XT barcode-adapter sequences alone are stated to be “available upon request” in the plate-library methods and remain unpublished in the checked paper and preprint. The schematic must not infer their individual indexes.

## Final readout

The plate workflow is full-length SMART-seq2-like transcript capture from both nuclear RNA and ribosome-bound exterior mRNA. Standard Nextera sequencing primers bind the tagmented endpoint; the unpublished custom barcode identities are left symbolic.
