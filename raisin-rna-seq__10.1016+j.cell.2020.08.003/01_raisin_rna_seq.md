# RAISIN RNA-seq — intact nuclei with ribosome-bound RNA

## Sources read

- 🟢 Drokhlyansky *et al.*, “The human and mouse enteric nervous system at single cell resolution,” *Cell* (2020), DOI `10.1016/j.cell.2020.08.003`, PMCID `PMC8358727` (defining paper and STAR Methods).
- 🟢 Shekhar *et al.*, *Cell* (2016), DOI `10.1016/j.cell.2016.07.054`, is the plate-library procedure explicitly cited by the defining paper's STAR Methods. Picelli *et al.*, SMART-seq2, DOI `10.1038/nmeth.2639`, is the earlier chemistry family; it is not the procedure directly cited at this point.

## Molecular interpretation

- 🟢 RAISIN means **Ribosomes And Intact Single Nucleus**. Its CST extraction buffer contains 0.49% CHAPS, 146 mM NaCl, 1 mM CaCl2, 21 mM MgCl2 and 10 mM Tris pH 8.0; the isolated nucleus retains rough ER, attached ribosomes and ribosome-bound RNA.
- 🟢 In the plate implementation, one RAISIN is sorted into TCL with 1% β-mercaptoethanol. RNA is purified with RNAClean XP, reverse-transcribed from oligo-dT with Maxima reverse transcriptase and a locked TSO, and amplified for 21 cycles with KAPA HiFi.
- 🟢 Amplified cDNA is converted with Nextera XT plus custom barcode adapters, pooled 384 wells at a time and sequenced on NextSeq 500.
- 🟢 The paper also implements droplet RAISIN RNA-seq on GemCode v2. The separate MIRACL-seq page shows that endpoint because deliberate droplet overloading is MIRACL-seq's defining change.
- 🟢 The paper's Key Resources Table prints the reverse-transcription primer as `5′-AAGCAGTGGTATCAACGCAGAGTACT(30)VN-3′`, the locked template-switching oligo as `5′-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3′`, and the PCR primer as `5′-AAGCAGTGGTATCAACGCAGAGT-3′` (PMC8358727 full-text XML, Oligonucleotides table). Preserve `rG` and `+G` as written: they denote the ribonucleotide and locked-nucleotide chemistry at the TSO 3′ end.
- 🔴 The custom Nextera XT barcode-adapter sequences alone are stated to be “available upon request” in the plate-library methods and remain unpublished in the checked paper and preprint. The schematic must not infer their individual indexes.

## Final readout

The plate workflow is full-length SMART-seq2-like transcript capture from both nuclear RNA and ribosome-bound exterior mRNA. Standard Nextera sequencing primers bind the tagmented endpoint; the unpublished custom barcode identities are left symbolic.
