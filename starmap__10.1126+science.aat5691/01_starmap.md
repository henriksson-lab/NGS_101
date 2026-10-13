# STARmap — source notes

**Defining source:** Wang et al., *Science* (2018),
[doi:10.1126/science.aat5691](https://doi.org/10.1126/science.aat5691).

## Evidence used

- 🟢 The defining paper's supplementary methods describe each SNAIL pair as a
  19–25-nt primer probe and 19–25-nt padlock probe binding adjacent RNA regions,
  with 6-nt complementary probe ends and one ligation site.
- 🟢 The primer probe carries the 12-nt template that receives both 6-nt padlock
  ends. It therefore splints the padlock's single closure junction and remains as
  the primer for rolling-circle amplification; it is not incorporated into the circle.
- 🟢 The same methods specify four SNAIL pairs per gene, a 5-base gene code,
  T4 DNA ligase, phi29 rolling-circle amplification with aminoallyl-dUTP,
  polyacrylamide embedding, proteinase-K clearing, and six SEDAL cycles.
- 🟢 In SEDAL an 11-nt partly degenerate reading probe and an 8-nt fluorescent
  decoding probe form a 19-bp ligated duplex only at a correct interrogation base.
- 🟡 The schematic uses symbolic target arms and code bases because the defining
  supplement available to this repository describes their design but does not expose
  the per-gene sequence tables in machine-readable form.
- 🟢 A later STARmap-family protocol, doi:10.1038/s41596-025-01248-3, is retained as
  supporting evidence for the named SEDAL reagents, but its oligos are not presented
  as if they were the original 2018 panel.

## Boundary

STARmap produces spatially anchored DNA rolonies read by microscopy and iterative
ligation chemistry. It does not produce an Illumina R1/I1/I2/R2 library.
