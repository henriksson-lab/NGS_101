# Spear-ATAC — source intake

Defining paper: Pierce et al., *Nature Communications* (2021), DOI
[10.1038/s41467-021-23213-w](https://doi.org/10.1038/s41467-021-23213-w).

🟢 The protocol produces chromatin-accessibility and guide-identity libraries that retain
the same droplet cell barcode.

## Exact oligos and library ends

🟢 Supplementary Note 3 and the paper's Methods print the guide-enrichment oligos:
`oSP2053` is 5′-biotinylated
`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGTATCCCTTGGAGAACCACCTTG`,
`oSP1594` is `AATGATACGGCGACCACCGAGA`, and the reverse PCR primer is
`CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTG`.
The custom guide Read 1 primer `oMCB1672` is
`GCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC`.

🟢 The same supplement prints the 10x bead oligo architecture as P5, a 16-nt cell
barcode, then S5; and prints the Nextera Read 1, Read 2, ATAC Index 1, guide Index 1,
and cell-barcode Index 2 primers. These data determine both completed library shells
without inference. The 16-nt cell barcode, 10-nt optional guide barcode, and 8-nt sample
index remain typed variable fields because the protocol uses barcode sets rather than
one fixed member.
