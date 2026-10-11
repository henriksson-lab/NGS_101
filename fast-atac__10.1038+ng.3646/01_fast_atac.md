# Fast-ATAC — research notes

Corces et al., *Nature Genetics* (2016),
[doi:10.1038/ng.3646](https://doi.org/10.1038/ng.3646).

## Sources read

- 🟢 The defining article and Online Methods, available as the PMC author manuscript
  (PMCID: PMC5042844).
- 🟢 Buenrostro et al., *Current Protocols in Molecular Biology* (2015),
  [doi:10.1002/0471142727.mb2129s109](https://doi.org/10.1002/0471142727.mb2129s109),
  for the amplification chemistry reused by Fast-ATAC.

## The one-step change

1. 🟢 Five thousand purified primary blood cells are FACS-sorted, pelleted at 500 r.c.f.
   for 5 min at 4 °C, and the supernatant is removed. There is no separate lysis step.
2. 🟢 The pellet is resuspended directly in 50 µl containing 25 µl 2× TD buffer,
   2.5 µl TDE1, 0.5 µl 1% digitonin and 22 µl nuclease-free water. Digitonin therefore
   permeabilizes the plasma membrane while Tn5 enters and fragments accessible chromatin.
3. 🟢 The one-tube permeabilization/transposition is incubated for 30 min at 37 °C with
   mixing at 300 r.p.m., then stopped by MinElute purification.
4. 🟢 The defining workflow uses five initial amplification cycles and an aliquot qPCR
   curve to select the remaining cycles rather than applying one fixed PCR count.

Fast-ATAC is therefore a sample-preparation variant, not a new adapter family. It is kept
as its own protocol page because simultaneous digitonin permeabilization and tagmentation
is its defining operation and its intended input is unusually small primary-cell samples.

## Library and reads

🟢 TDE1 installs the standard 19-base Nextera mosaic end. The original Ad1_noMX and one
indexed Ad2 primer then select heterologous S5/S7 fragments, fill the nine-base staggered
transposition gaps, and complete P5 and P7. The example Ad2.1 ordering sequence contains
`TCGCCTTA`, corresponding to the sample-sheet Index 1 value `TAAGGCGA`; there is no i5
index in the defining construct.

🟢 The study sequenced 50-base paired-end libraries on an Illumina HiSeq 4000. Read 1 and
Read 2 enter accessible genomic DNA from opposite Tn5 mosaic ends, while the eight-base
Index 1 read identifies the sample. No molecular UMI or cell barcode is present.
