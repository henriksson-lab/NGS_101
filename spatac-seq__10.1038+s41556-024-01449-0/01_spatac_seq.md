# SPATAC-seq — source notes

## What it is

SPATAC-seq is the **split-pool ligation-based assay for transposase-accessible chromatin
using sequencing**. Fixed nuclei receive a two-ended first barcode during custom-Tn5
transposition, then one barcode is ligated to each end in two successive 48-well
pool-and-split rounds. Indexed PCR creates sublibraries and completes the Illumina arms.

This method is unrelated to Spatial-ATAC-seq despite the similar spelling.

## Sources read

- 🟢 Sun et al., *Nature Cell Biology* (2024),
  [doi:10.1038/s41556-024-01449-0](https://doi.org/10.1038/s41556-024-01449-0), is the
  defining paper. Figure 1 gives the molecular architecture and the four indexing rounds.
- 🟢 Supplementary Table 1 prints every Tn5, ligation-splint, round-2, round-3, i5 and i7
  oligo and identifies the `5′P` and `3′ddC` modifications.
- 🟢 The authors' [Lan-lab/SPATAC-seq](https://github.com/Lan-lab/SPATAC-seq) repository
  supplies the workflow drawing, the 110,592-entry 32-nt valid-barcode list and executable
  read preprocessing. Its script extracts round-2 at Read 1 cycles 1–8, the forward half
  of round-1 at Read 1 cycles 39–46, round-3 at Read 2 cycles 1–8 and the reverse half of
  round-1 at Read 2 cycles 24–31. Genomic sequence begins at Read 1 cycle 66 and Read 2
  cycle 51.
- 🟢 GEO series [GSE215901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215901)
  independently describes the 48 indexed transposomes, two 48-well ligation rounds and
  indexed sublibrary PCR.

## Molecular reconstruction

🟢 `Tn5_F_D1` and `Tn5_R_D1` each carry an 8-nt round-1 barcode between a 15-nt leader
and the 19-nt mosaic end. The two halves are retained as distinct structured parts of the
same combinatorial cell barcode. The authors use 48 specified forward/reverse pairs rather
than all 64 possible combinations of the eight oligos on each side.

🟢 The round-2 oligo is
`ACGACGCTCTTCCGATCT-[8 nt]-ATCCACGTGCTTGAG`. Its 3′ end ligates to the
5′-phosphorylated forward round-1 adapter. `Tn5_F_LA` is a 23-nt splint: its reverse
complement is exactly the 15-nt round-2 bridge followed by the first eight bases of the
round-1 leader. The page therefore computes and marks that ligation boundary.

🟢 The round-3 oligo is `S7-mosaic end-[8 nt]`. Its 3′ end ligates to the
5′-phosphorylated reverse round-1 adapter. Each `Tn5_R_LA-*` splint contains eight bases
complementary to the round-1 leader, the reverse complement of the corresponding
round-3 barcode and seven bases of the round-3 mosaic end. This second ligation boundary
is likewise anchored to named segments.

🟢 Round-4 primers add P5, i5 and a complete TruSeq Read 1 site on the forward end, and
P7, i7 and S7 on the reverse end. The result intentionally combines a TruSeq Read 1 entry
with Nextera Index 1 and Read 2 entries.

## Read layout

🟢 The authors' processing recipe uses raw Read 1 bases 66–115 and raw Read 2 bases
51–100 as 50-nt genomic mates. The discarded prefixes are molecular information, not
arbitrary trimming: Read 1 is `round2(8) + bridge(30) + round1-forward(8) + ME(19)`, and
Read 2 is `round3(8) + leader(15) + round1-reverse(8) + ME(19)`.

🟡 The page labels the full demonstrated raw layouts as Read 1 = 115 cycles and Read 2 =
100 cycles so the structured-feature renderer can report all barcode coordinates. The
paper's construct, rather than those total cycle counts, is the invariant: longer reads
can continue farther into genomic DNA without changing the barcode positions.
