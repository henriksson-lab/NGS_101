# HiCAR — research notes

Wei et al., *Molecular Cell* (2022),
[doi:10.1016/j.molcel.2022.01.023](https://doi.org/10.1016/j.molcel.2022.01.023).

## Sources read

- 🟢 The defining article and STAR Methods (PMCID: PMC8934281).
- 🟡 Original Table S1 and other supplements are manual/gated entries in the current bundle; earlier blanket supplement-reading provenance is not reverified by this pass.
- 🟢 The authors’ expanded protocol, Lu et al., *Current Protocols* (2023),
  [doi:10.1002/cpz1.899](https://doi.org/10.1002/cpz1.899), including its oligo table.

The defining molecular boundary is the 2022 CviQI/NlaIII implementation. Exact oligos below are transcribed from the later 2023 protocol; applicability to the original 2022 Table S1 remains unresolved, so the final arms are inferred. The expanded
protocol changes the first digest to MseI; that change is useful operationally but is not
silently projected backward onto this page.

## Later-protocol oligos — defining-source assignment unresolved

- 🟢 HiCAR7-i0 is written `/5Phos/taGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGGAGATGTGTATAAGAGACAG`.
  The leading lower-case `ta` is retained because it is the ligatable flank, not a display
  convention; the terminal 19 bases are the Tn5 mosaic end.
- 🟡 The earlier TruseqR1 splint assignment is not reverified: the defining Table S1 is gated and the inspected later tables do not establish that string. The graph leaves this oligo unresolved.
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
5. 🟢 Treat the purified circle pool with PmeI for 1 h at 37 °C, then perform
   12-cycle Q5 inverse PCR with the splint-side and mosaic-end-side primers
   (defining Methods, “Step 5. PmeI digestion and PCR”). The source does not
   identify every PmeI cleavage product. PCR yields the final linear
   flow-cell-compatible library.

## Reads

🟢 Read 1 derives from CviQI-digested interacting genomic DNA. Read 2 derives from the
Tn5-tagmented accessible anchor. Thus a pair is asymmetric by design: R2 supports the
one-dimensional open-chromatin signal and R1 supplies its spatial partner. The i501 and
i7-10-L sequences provide the two sample indexes; no UMI or cell barcode is present.

## Bounded source check — 2026-10-17

Read the published defining article (PMC8934281 author-manuscript JATS): Design
text line 84, transposome assembly 192, detailed Steps 1–5 198–216 and
Figure S3 caption 284. Read later published author protocol PMC10575683
Tables 1/2 (1388–1446), adapter assembly (398–422), overview 434 and changed
MseI/splint2 digest 682. Ran the bundle oligo scraper and followed HiCAR7-i0
and both PCR-primer hits to those tables. These are two originals, not two
independent audit readers. No chemistry facts were taken from upstream drawings.

Added source-stated EDTA stop/wash/SDS–Triton preparation, T4 polymerase heat
inactivation, 0.9X SPRI recovery after NlaIII, purification after circularization,
and final 300–750-bp gel extraction. Named T4 DNA ligase and Q5. Recovery steps
remain summaries where exact chemical states are unresolved. Candidate Tn5
products are unresolved downstream fates, not physically discarded just for
lacking a selected anchor configuration. Gel exclusion has a represented separation.

The original Table S1 PDF and other supplement workbooks remain manifest manual
entries, so no identity equivalence with the later oligo table is asserted.
HiCAR7-i0 and the two PCR primers are shown with explicit later-source provenance;
the compatible final arms are inferred. The source names original primers i501
and i7-10-L, but matching names alone does not verify exact variant equivalence.
The older 20-nt splint string is retained only as a historical module constant,
removed from the displayed oligo list and its unsupported transcription check.
PCR geometry is unavailable on the actual PmeI-treated template; exact landing
on an illustrative final amplicon does not validate inverse PCR from that template.
Removed assumed 150/8/8/150 acquisition counts. R1/R2 roles are source-supported;
exact acquired cycle counts and index configuration remain open.

Deferred: original 2022 Table S1 and remaining supplements, exact transferred
adapter/splint and junction geometry (including PmeI/NlaIII ends), historical
Tn5/predecessor 4C/vendor originals, complete off-product propagation and selected
study run. Do not substitute the later MseI method or call the unavailable
original oligos confirmed absent. This is a capped pass, not complete validation.
