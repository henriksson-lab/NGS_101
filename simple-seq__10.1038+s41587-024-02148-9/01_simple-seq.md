# SIMPLE-seq — source intake

Defining paper: Bai et al., *Nature Biotechnology* (2024), DOI
[10.1038/s41587-024-02148-9](https://doi.org/10.1038/s41587-024-02148-9).

## Sources read and oligos

- 🟢 The publisher's Supplementary Tables 1–8 workbook, sheet
  `Supplementary_Table1_1`, prints the library-building oligos with their chemical
  modifications. Sheets `Supplementary_Table1_2` through `_1_4` print the first,
  second and third rounds of cell-barcoding oligos. The workbook is retrievable from
  the article's Supplementary information and through `tools/get_sources.py`.
- 🟢 The publisher's supplementary PDF is a three-page reporting summary; visual
  inspection shows it contains no chemistry protocol. The publisher workbook also
  contains model oligos with 5mC/5hmC controls in `Supplementary_Table1_5` and DNA
  integrity qPCR primers in `_1_6`, separate from the library-building rows below.
- 🟡 The main Methods article is subscription-only in this source intake. The
  reaction order and placement of every oligo still require original Methods review
  before a base-by-base schematic is built.

The exact named oligos in Supplementary Table 1 are `Linker-R01`
`CTCGAATGCCTCAAGCACGTGGAT`, `Linker-R02`
`GGTCTGAGTTCGCACCGAAACATCGGCCAC`, `Blocker-R01`
`ATCCACGTGCTTGAGAGGCCAGAGCATTCG` and `Blocker-R02`
`GTGGCCGATGTTTCGGTGCGAACTCAGACC`. The modified oligos are
`pMENTs` `/5Phos/CTGTCTCTTATACACATCT/ddC/`, bridge-adapter bottom
`AGATGTGTATAAGAGACAGCTGACGNNNNNN/-NH2`, bridge-adapter top
`/5Phos/CGTCAGCTGTCTCTTATACACATCTGACGCTGCCGACGA`, and the indicator
primer `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG/caC/TGA/caC/G` (`caC` is
5-carboxylcytosine). The final P7 index primer is
`CAAGCAGAAGACGGCATACGAGAT[6nt i7 index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC`;
the N5 universal primer is `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTC`.

🟢 The first-round representative `DNA_#01` is
`GCATTCGAGACGCAAGATGTGTATAAGAGACAG`. The second-round representative
`R01_#01` is `/5Phos/GTGCGAACTCAGACCAAACCGGATCCACGTGCTTGAG`.
The complete barcode inventories are published in the workbook; one example is
insufficient to infer the chosen oligo in any particular indexed well.
The third-round representative `R02_#01` is
`CAGACGTGTGCTCTTCCGATCTNNNNNNNNNNAAACCGGGTGGCCGATGTTTCG`;
its ten random bases must be kept distinct from the three well-barcode rounds.
