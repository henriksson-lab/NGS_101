# scSPRITE — research notes

## Sources read

- 🟢 Arrastia et al. preprint, DOI [10.1101/2020.08.11.242081](https://doi.org/10.1101/2020.08.11.242081), defining source used by the catalogue.
- 🟢 Arrastia et al. 2022, DOI [10.1038/s41587-021-00998-1](https://doi.org/10.1038/s41587-021-00998-1), PMC11588347, published full text, Methods and supplements.

## Molecular order

🟢 Intact permeabilized nuclei receive DPM, Odd and Even tags in three 96-well
split-pool ligations. These three parts identify a cell. Filtered nuclei are then lysed;
crosslinked chromatin complexes are sonicated, coupled to beads and receive Odd, Even and
Y-even tags in three more rounds. The complete six-part identifier defines a spatial
cluster while its first three parts still define the cell.

🟢 Paired-end sequencing uses at least 120 bp in Read 1 for genomic DNA and DPM, and at
least 95 bp in Read 2 for Odd–Even–Odd–Even–Y-even. Only complete six-tag molecules were
retained by the published analysis.

## Protocol boundary

This is not bulk SPRITE with a cell label added afterward. The first barcoding phase is
performed while nuclei are intact, and the second phase occurs after selecting nuclei and
isolating spatial complexes. It therefore has its own page and source identity.

## Oligo source review

🟢 The preprint's Methods and published article specify the six-tag order as
DPM–Odd–Even–Odd–Even–Y-even. The first three tags identify a cell; all six
identify a spatial complex. The original 2018 SPRITE Supplementary Table S5,
recovered from the authors' Caltech archive, prints the complete DPM, Odd,
Even and Terminal oligos. The [scSPRITE author parser](https://github.com/caltech-bioinformatics-resource-center/Guttman_Ismagilov_Labs)
lists 96 read-barcode strings each for DPM, Odd, Even and Y; its 96 Y strings
exactly match the 2018 Table S5 Terminal Bottom barcode column. 🟡 A match in
read-barcode strings does not prove that scSPRITE used the 2018 full Terminal
duplex for its Y-even ligation. The scSPRITE-specific Y-even ordering oligo,
end geometry and exact library PCR primer pair remain unresolved.

🟢 The parser's Odd read strings are the Table S5 Odd Bottom barcodes without
their first base. Its Even strings are likewise the Table S5 Even Bottom
barcodes without the first base, although the 96 entries are ordered
differently in the two files. Its DPM strings are the reverse complements of
the Table S5 DPM Top eight-base barcodes followed by `TT`. These are read
parsing relationships, not evidence for a complete ligation-ready duplex.
The adjacent author parser configuration uses the same Y, Odd and Even read
strings but prepends a constant 17-base motif to every DPM read string. The
paper does not identify which parser setting applied to every published
library, so the exact raw-read offset of DPM remains open.
