# Fast-ATAC — research notes

Corces et al., *Nature Genetics* (2016),
[doi:10.1038/ng.3646](https://doi.org/10.1038/ng.3646).

## Sources read

- 🟢 The defining article and Online Methods, available as the PMC author manuscript
  (PMCID: PMC5042844).
- 🟢 Buenrostro et al., *Current Protocols in Molecular Biology* (2015),
  [doi:10.1002/0471142727.mb2129s109](https://doi.org/10.1002/0471142727.mb2129s109),
  for the amplification workflow reused by Fast-ATAC. The defining article
  explicitly cites Buenrostro et al., *Nature* (2015),
  [doi:10.1038/nature14590](https://doi.org/10.1038/nature14590), for its
  **modified primers**. That paper's MOESM36 workbook prints the dual-index
  v2_Ad1 and v2_Ad2 series. The older single-index Ad1_noMX and Ad2 series
  appear in Buenrostro et al., *Nature Methods* (2013),
  [doi:10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688), MOESM200;
  they do not match the Fast-ATAC run layout.

## The one-step change

1. 🟢 Five thousand purified primary blood cells are FACS-sorted, pelleted at 500 r.c.f.
   for 5 min at 4 °C, and the supernatant is removed. There is no separate lysis step.
2. 🟢 The pellet is resuspended directly in 50 µl containing 25 µl 2× TD buffer,
   2.5 µl TDE1, 0.5 µl 1% digitonin and 22 µl nuclease-free water. Digitonin therefore
   permeabilizes the plasma membrane while Tn5 enters and fragments accessible chromatin.
3. 🟢 The one-tube permeabilization/transposition is incubated for 30 min at 37 °C with
   mixing at 300 r.p.m. (TD buffer/TDE1 Illumina FC-121-1030; digitonin Promega G9441),
   then purified with a QIAGEN MinElute Reaction Cleanup kit (Cat# 28204) and eluted in
   10 µl elution buffer (10 mM Tris-HCl, pH 8).
4. 🟢 The article itself says only that transposed fragments "were amplified and
   purified as described previously" (its ref. 48, Buenrostro et al. 2015 *Curr Protoc
   Mol Biol*) "with modified primers" (its ref. 23, Buenrostro et al. 2015 *Nature*), and
   that libraries were quantified by qPCR before sequencing.
5. 🟡 Inherited from ref. 48 (PMC4374986 Basic Protocol, PCR Amplification; read
   2026-10-17), not restated in the Fast-ATAC Methods: 10 µl transposed DNA with 2.5 µl
   each 25 µM primer in NEBNext High-Fidelity 2× PCR Master Mix; 72 °C 5 min gap-fill
   then 98 °C 30 s; five cycles of 98 °C 10 s, 63 °C 30 s, 72 °C 1 min; a 5 µl SYBR
   Green I qPCR side reaction picks the extra cycles N at ¼ of maximum fluorescence; the
   remaining 45 µl is cycled N more times and purified on a MinElute PCR Purification
   column (20 µl elution). The page shows this as the cited workflow, not as an
   independently reported Fast-ATAC parameter set.

Fast-ATAC is therefore a sample-preparation variant, not a new adapter family. It is kept
as its own protocol page because simultaneous digitonin permeabilization and tagmentation
is its defining operation and its intended input is unusually small primary-cell samples.

## Library and reads

🟢 TDE1 installs the standard 19-base Nextera mosaic end. The cited 2015
v2_Ad1 and v2_Ad2 primer families then select heterologous S5/S7 fragments,
fill the staggered transposition gaps, and complete P5 and P7. The illustrative
v2_Ad1.1 primer carries an eight-base i5 (`TAGATCGC`), while v2_Ad2.1 carries
an eight-base i7 ordering sequence (`TCGCCTTA`, sample-sheet Index 1
`TAAGGCGA`). 🟡 These are family examples, not assigned indexes for any
particular Fast-ATAC sample.

🟢 The defining article's Online Methods specify NextSeq paired-end,
dual-index `76×8×8×76` sequencing. Read 1 and Read 2 enter accessible
genomic DNA from opposite Tn5 mosaic ends; the eight-base i7 and i5 reads
identify the selected sample indexes. No molecular UMI or cell barcode is present.

## Open questions

- 🔴 The per-sample v2_Ad1/v2_Ad2 index pairs are not assigned in the inspected article
  text or supplements; both indexes stay placeholders.
- Whether every Fast-ATAC library followed the ref. 48 cycling exactly (for example
  the post-transposition MinElute Reaction Cleanup kit here versus the MinElute PCR
  Purification kit in ref. 48) is not stated; these are different protocol texts, not
  a published discrepancy.

## Review trail

- 2026-10-17 bounded pass (round 2026-10-17-breadth-first-10): re-read PMC5042844
  Online Methods "Fast-ATAC sequencing" (JATS text lines 194–198) and references 23/48;
  read ref. 48 PCR Amplification (PMC4374986 text lines 103–169); verified v2_Ad1.1 and
  v2_Ad2.1 against Buenrostro 2015 *Nature* MOESM36 Sheet1 lines 6 and 99. The graph now
  shows the cold pellet as a cleanup (not a separate permeabilization), digitonin in the
  transposition step, the qPCR-chosen extra cycles and the final MinElute cleanup.
