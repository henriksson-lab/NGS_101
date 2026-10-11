# RAISIN RNA-seq — intact nuclei with ribosome-bound RNA

## Sources read

- 🟢 Drokhlyansky *et al.*, “The human and mouse enteric nervous system at single cell resolution,” *Cell* (2020), DOI `10.1016/j.cell.2020.08.003`, PMCID `PMC8358727` (defining paper and STAR Methods).
- 🟢 Picelli *et al.*, SMART-seq2, DOI `10.1038/nmeth.2639` (the library method explicitly named by the defining paper).

## Molecular interpretation

- 🟢 RAISIN means **Ribosomes And Intact Single Nucleus**. Its CST extraction buffer contains 0.49% CHAPS, 146 mM NaCl, 1 mM CaCl2, 21 mM MgCl2 and 10 mM Tris pH 8.0; the isolated nucleus retains rough ER, attached ribosomes and ribosome-bound RNA.
- 🟢 In the plate implementation, one RAISIN is sorted into TCL with 1% β-mercaptoethanol. RNA is purified with RNAClean XP, reverse-transcribed from oligo-dT with Maxima reverse transcriptase and a locked TSO, and amplified for 21 cycles with KAPA HiFi.
- 🟢 Amplified cDNA is converted with Nextera XT plus custom barcode adapters, pooled 384 wells at a time and sequenced on NextSeq 500.
- 🟢 The paper also implements droplet RAISIN RNA-seq on GemCode v2. The separate MIRACL-seq page shows that endpoint because deliberate droplet overloading is MIRACL-seq's defining change.
- 🔴 The locked TSO and custom barcode-adapter sequences were not printed (“available upon request”), so the schematic does not invent them.

## Final readout

The plate workflow is full-length SMART-seq2-like transcript capture from both nuclear RNA and ribosome-bound exterior mRNA. Standard Nextera sequencing primers bind the tagmented endpoint; the unpublished custom barcode identities are left symbolic.
