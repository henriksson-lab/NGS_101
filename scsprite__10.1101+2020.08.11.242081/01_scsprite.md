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
DPM–Odd–Even–Odd–Even–Y-even. Barcode assignment uses a known tag-sequence
table, but neither accessible main text prints the ordering oligos, terminal
Y-even sequence, or library PCR primer pair. The first three tags identify a
cell; all six identify a spatial complex. 🟡 Bulk SPRITE's 2018 tag table and
the authors' 2022 detailed SPRITE workbook may explain the shared DPM/Odd/Even
chemistry, but they do not establish the scSPRITE-specific Y-even adapter without
an explicit match to this protocol. Exact end geometry remains open.
