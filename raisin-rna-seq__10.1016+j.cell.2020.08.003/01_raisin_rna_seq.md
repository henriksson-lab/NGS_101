# RAISIN RNA-seq — intact nuclei with ribosome-bound RNA

## Sources read

- 🟢 Drokhlyansky *et al.*, “The human and mouse enteric nervous system at single cell resolution,” *Cell* (2020), DOI `10.1016/j.cell.2020.08.003`, PMCID `PMC8358727` (defining paper and STAR Methods).
- 🟢 Shekhar *et al.*, *Cell* (2016), DOI `10.1016/j.cell.2016.07.054`, is the plate-library procedure explicitly cited by the defining paper's STAR Methods. Picelli *et al.*, SMART-seq2, DOI `10.1038/nmeth.2639`, is the earlier chemistry family; it is not the procedure directly cited at this point.

- Bounded pass 2026-10-17 (round breadth-first-13): published PMC8358727 JATS STAR Methods lines 1032–1039 (CST/TST composition), 1041–1054 (FACS into TCL/1% β-ME), 1056–1072 (whole-transcriptome amplification, library construction, sequencing), Key Resources 2915–2972 (Maxima H Minus EP0753, KAPA HiFi KK2602, RNAClean, AMPure XP, Nextera XT FC-131-1096, NextSeq 500 high output v2 75-cycle FC-404-2005) and 3030–3041 (Oligonucleotides); cited Shekhar *et al.* 2016 PMC5003425 JATS “Plate-based RNA-seq experiments” lines 288–294 (2.2× RNAClean, 21 cycles, 0.075 ng cDNA, ¼ Nextera XT volumes). Scraper hits at PMC8358727 lines 3225–3246 are unrelated Key-Resources template placeholders (siRNA, morpholino, GFP primer), and the Shekhar table at lines 276–284 holds Drop-seq oligos, not plate-library oligos. Figures S1–S7 and Tables S1–S6 remain PMC download-gated (unread).

## Molecular interpretation

- 🟢 RAISIN means **Ribosomes And Intact Single Nucleus**. Its CST extraction buffer contains 0.49% CHAPS, 146 mM NaCl, 1 mM CaCl2, 21 mM MgCl2 and 10 mM Tris pH 8.0; the isolated nucleus retains rough ER, attached ribosomes and ribosome-bound RNA.
- 🟢 In the plate implementation, one RAISIN is sorted into TCL with 1% β-mercaptoethanol. RNA is purified with RNAClean XP, reverse-transcribed from oligo-dT with Maxima reverse transcriptase (Key Resources: Maxima H Minus, EP0753) and a locked TSO, amplified for 21 cycles with KAPA HiFi HotStart ReadyMix, and the cDNA is purified twice with AMPure XP.
- 🟡 The cited Shekhar *et al.* plate procedure (the base the paper modifies) used 2.2× RNAClean SPRI beads without final elution, 0.075 ng cDNA per cell and ¼ of the standard Nextera XT volume for tagmentation and final PCR. The RAISIN paper does not restate these values, so they are inherited by citation, not confirmed for RAISIN libraries.
- 🟢 Amplified cDNA is converted with Nextera XT plus custom barcode adapters, pooled 384 wells at a time and sequenced on NextSeq 500 (Key Resources: high output kit v2, 75 cycles). 🟡 Acquired read/index allocation is unassigned in this model; the named kit alone does not assign a mask.
- 🟢 The paper also implements droplet RAISIN RNA-seq on GemCode v2. The separate MIRACL-seq page shows that endpoint because deliberate droplet overloading is MIRACL-seq's defining change.
- 🟢 The paper's Key Resources Table prints the reverse-transcription primer as `5′-AAGCAGTGGTATCAACGCAGAGTACT(30)VN-3′`, the locked template-switching oligo as `5′-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3′`, and the PCR primer as `5′-AAGCAGTGGTATCAACGCAGAGT-3′` (PMC8358727 full-text XML, Oligonucleotides table). Preserve `rG` and `+G` as written: they denote the ribonucleotide and locked-nucleotide chemistry at the TSO 3′ end.
- 🔴 The plate-library Methods state that custom Nextera XT barcode-adapter sequences are “available upon request”. Their exact members remain unassigned here; this statement is not a completed publication-wide disclosure audit.

## Final readout

The plate workflow is full-length SMART-seq2-like transcript capture from both nuclear RNA and ribosome-bound exterior mRNA. The custom Nextera XT endpoints, primer sites and barcode identities remain unassigned; no generic canonical endpoint or acquired read mask is substituted.

## Native graph migration (2026-10-18, batch04)

Inspected: Published2020 PMC8358727 Methods1032–1076 and KeyResources oligos3032–3040; explicitly cited Shekhar2016 PMC5003425 plateMethods288–294, contrasting adjacent Drop-seq reagents/run choices. Mandatory whole-bundle scraper completed; no preprint/vendor recipe substituted.

Removed queued Nextera endpoint/defaults and rare-cell workflow constructors. Native selected plate graph separates nucleus/ribosomes/RNA components, TCL lysis and purified RNA populations. Source modifiedRT/lockedTSO/PCR oligos are retained; a clearly inferred representative switched ssDNA preserves copied handle and known source primer relation, without attaching the physical TSO covalently. RNA templates and remaining TSO are molecular/unresolved. Actual switched ssDNA enters21cycle preamplification, whose source duplex products enter twoAMPure recoveries. Native tagged-strand classes, rather than expected final XT library, enter indexing PCR. Nonamplifying products and unselected particles are unresolved rather than falsely discarded. Custom XT endpoints have unknown sample-index features/widths; no default8nt indices or fixed36nt fragment length. Source75-cycle kit does not establish acquired read/index split. GemCode droplet RAISIN and cited predecessor run masks remain separate.

Deferred: Native switch junction heterogeneity/first-copy proof, exact custom XT adapter/index/primer members and gap geometry; exact acquired read/index allocation and native terminal cleanup scope remain unassigned. Shekhar0.075ng/quarter reactions are explicitly inherited plate scope, not Drop-seq primers or run recipes. This capped pass does not establish complete molecular proof or publication-wide absence.
