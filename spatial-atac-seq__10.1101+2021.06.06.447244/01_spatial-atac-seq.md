# Spatial-ATAC-seq — source notes

## What it is

Spatial-ATAC-seq performs Tn5 tagmentation in a fixed tissue section and then ligates two
orthogonal microfluidic barcode sets. The ordered A×B combination records the tissue pixel
for each accessible-chromatin fragment.

## Sources read

- 🟢 Deng et al., *Nature* (2022),
  [doi:10.1038/s41586-022-05094-1](https://doi.org/10.1038/s41586-022-05094-1)
  (PMC9452302), is the defining published paper; the catalogue directory retains its
  earlier defining-preprint DOI
  [10.1101/2021.06.06.447244](https://doi.org/10.1101/2021.06.06.447244), which is
  superseded and not used for active claims. Read: Methods "Preparation of the
  transposome", "DNA oligos, DNA barcode sequences and other key reagents",
  "Spatial-ATAC-seq profiling" and "Data preprocessing" (JATS text lines 213–235), and
  Supplementary Table 3 (MOESM5 workbook: sheets "DNA oligos" and "DNA barcode
  sequences"). Supplementary Table 4 (MOESM6) lists reagents and catalogue numbers.
- 🟢 Farzad et al., *Nature Protocols* (2024),
  [doi:10.1038/s41596-024-01013-y](https://doi.org/10.1038/s41596-024-01013-y), is the
  authors' later detailed protocol for Spatial-ATAC-seq and Spatial-CUT&Tag. It is not
  used to override the 2022 values; its Supplementary Table 1 conflicts with them (below).
- 🟡 The authors' `dyxmvp/Spatial_ATAC-seq` processing repository (cited in the 2022
  Methods) is reported in earlier passes to parse barcode B from raw Read 1 bases 1–8 and
  barcode A from bases 39–46; the model reproduces those coordinates from the oligos.
  The repository itself was not re-read in the 2026-10-17 bounded pass.

## Oligos (Deng et al. 2022, verbatim)

Transposome (Methods): Tn5MErev `/5Phos/CTGTCTCTTATACACATCT`; Tn5ME-A
`/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG`; Tn5ME-B
`GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`. Supplementary Table 3: ligation linker 1
`AGTCGTACGCCGATGCTGTCTCTTATACACATCT`; ligation linker 2
`CGAATGCTCTGGCCTCTCAAGCACGTGGAT`; first barcodes `Tn5ME-A_BC1_1…50`
(`/5Phos/AGGCCAGAGCATTCG`-8-nt barcode-`AGATGTGTATAAGAGACAG`); second barcodes
`BC2_R1_1…50` (custom Read 1 site-8-nt barcode-`ATCCACGTGCTTGAG`); New-P5-SMART PCR
hybrid oligo `AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C`;
Ad2.1–Ad2.8 indexed primers; custom Read 1 primer
`GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`.

## Molecular reconstruction

🟢 The Methods describe bulk tagmentation of the permeabilized section (37 °C, 30 min,
EDTA stop), barcode A1–A50 pre-annealed to ligation linker 1 and ligated with T4 DNA
ligase in the first channels (37 °C, 30 min), the same for barcode B1–B50 with ligation
linker 2 in perpendicular channels, proteinase K lysis and overnight 65 °C crosslink
reversal, Zymo DCC-5 purification, PCR with New-P5-SMART and Ad2 (72 °C 5 min, 5 cycles,
then qPCR-determined additional cycles) and a 1× AMPure XP cleanup.

🟡 The splint geometry is computed from the printed oligos (the module raises if it does
not hold). Ligation linker 1 is the reverse complement of the barcode-A 3′ mosaic end
followed by `CATCGGCGTACGACT`, the 5′ extension of the 5′-phosphorylated Tn5ME-A strand
that Tn5 transfers onto genomic DNA. Linker 1 therefore joins barcode A's 3′ end to
Tn5ME-A, and the sequenced strand keeps both mosaic ends: barcode-A ME, Tn5ME-A
extension, Tn5ME-A ME, genomic DNA. Linker 2 is the reverse complement of the barcode-B
3′ end followed by the barcode-A 5′ end, so it joins barcode B to the phosphorylated
barcode-A oligo.

🟡 Model correction (2026-10-17): the page previously used the shared
`spatial_epigenome` construct, which placed linker 1's reverse complement directly
against genomic DNA and omitted the Tn5ME-A strand, so Read 1 reached genomic DNA one
mosaic-end length too early. The module now builds the construct locally.

🟡 Candidate Tn5ME-A/A and Tn5ME-B/B fragments follow from the mixed transposome; no
physical separation is described, so the graph leaves their fate unresolved.

## Published discrepancy

Deng et al. 2022 Methods print Tn5ME-A as `/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG`.
Farzad et al. 2024 Supplementary Table 1, used by its Procedure Steps 16–24 (transposome
assembly for spatial-ATAC-seq), prints Tn5ME-A as `/5Phos/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`.
Ligation linker 1, identical in both sources, pairs with the 2022 extension and not with
the 2024 one. The model uses the 2022 oligo; whether the 2024 entry belongs to another
experiment or is a table error is not stated.

## Read layout

🟢 HiSeq 4000 or NovaSeq 6000, paired-end 150, with the custom Read 1 primer. 🟡 Raw
Read 1 enters barcode B (cycles 1–8), linker 2, barcode A (cycles 39–46), the two mosaic
ends and the Tn5ME-A extension, and only then reaches genomic DNA (the page computes
the cycle from the construct). The i7 index is read from the S7 side and Read 2 enters
genomic DNA from the Tn5ME-B end. 🟢 The Methods convert raw Read 1 by filtering on the
two linkers and moving the genomic part to a new Read 1 and barcodes to a new Read 2.

## Open questions

- Whether the 2024 Farzad Table 1 Tn5ME-A reflects a different experiment.
- Ad2 index used per sample and acquired index-read length (8-nt i7 assumed from the Ad2
  names; the run's index cycles are not printed).
- The tagmented intermediate is drawn as a filled duplex; the single-stranded Tn5ME-A
  extension and the 9-nt Tn5 gaps before the 72 °C fill are described, not drawn.
