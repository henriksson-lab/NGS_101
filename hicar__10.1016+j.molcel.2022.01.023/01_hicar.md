# HiCAR — research notes

Wei et al., *Molecular Cell* (2022),
[doi:10.1016/j.molcel.2022.01.023](https://doi.org/10.1016/j.molcel.2022.01.023).

## Sources read

- 🟢 The defining article, STAR Methods and supplement (PMCID: PMC8934281).
- 🟢 The authors’ expanded protocol, Lu et al., *Current Protocols* (2023),
  [doi:10.1002/cpz1.899](https://doi.org/10.1002/cpz1.899), including its oligo table.

The defining molecular boundary is the 2022 CviQI/NlaIII implementation. The expanded
protocol changes the first digest to MseI; that change is useful operationally but is not
silently projected backward onto this page.

## Public oligos

- 🟢 HiCAR7-i0 is written `/5Phos/taGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGGAGATGTGTATAAGAGACAG`.
  The leading lower-case `ta` is retained because it is the ligatable flank, not a display
  convention; the terminal 19 bases are the Tn5 mosaic end.
- 🟢 The TruseqR1 splint is `ACACGACGCTCTTCCGATCT`.
- 🟢 NEB primer i501 is `AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
- 🟢 Nextera-pcr-i7-10-L is `CAAGCAGAAGACGGCATACGAGATCAGCCTCGGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`.
- 🟢 The source lists eight staggered HiCAR7 variants. The schematic shows i0 rather than
  normalizing variants whose stagger differs, and therefore does not imply that their
  sequences are interchangeable.

## Construct-changing steps

1. 🟢 Crosslink 30,000–100,000 cells and tagment accessible chromatin with Tn5 loaded
   with the engineered HiCAR7 adapter.
2. 🟢 Digest nuclei with 100 U CviQI for 2 h at room temperature in the presence of
   TruseqR1 splint, then perform in-situ T4 ligation for 4 h. The adapter’s single-stranded
   flank is thereby joined to spatially proximal restriction-cut genomic DNA.
3. 🟢 Reverse crosslinks, purify DNA, and use T4 DNA polymerase plus dNTPs for 30 min at
   room temperature to repair the staggered Tn5 gaps.
4. 🟢 Digest with NlaIII for 1 h at 37 °C, purify, dilute DNA to 1 ng/µl and circularize
   by intramolecular T4 ligation for at least 2 h at room temperature.
5. 🟢 Perform 12-cycle Q5 inverse PCR with the splint-side and mosaic-end-side primers.
   This opens selected circles into the final linear flow-cell-compatible library.

## Reads

🟢 Read 1 derives from CviQI-digested interacting genomic DNA. Read 2 derives from the
Tn5-tagmented accessible anchor. Thus a pair is asymmetric by design: R2 supports the
one-dimensional open-chromatin signal and R1 supplies its spatial partner. The i501 and
i7-10-L sequences provide the two sample indexes; no UMI or cell barcode is present.
