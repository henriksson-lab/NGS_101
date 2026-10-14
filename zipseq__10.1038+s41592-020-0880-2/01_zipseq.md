# ZipSeq — source notes

Defining paper: Hu et al., *Nature Methods* (2020), [doi:10.1038/s41592-020-0880-2](https://doi.org/10.1038/s41592-020-0880-2).

- 🟢 PMC manuscript, Methods and oligonucleotide list were read.
- 🟢 Antibody- or lipid-linked anchors are protected by NPOM-caged strands. Patterned 405-nm light uncages a selected region.
- 🟢 A complementary zipcode with a 30-adenosine tail hybridizes after each illumination round; rounds can be repeated with different cages/zipcodes.
- 🟢 Dissociated cells undergo 10x single-cell RNA-seq, which captures zipcode and endogenous transcript molecules in the same GEM.
- 🟢 The paper's Methods oligonucleotide list gives the antibody anchor as
  `Amine-C6–ATCGTTTTTTTTTTTTGATGTTGACGG`, the O1 blocking strand as
  `CGATCTGTGGTTGCTAC`, and zipcode 1 as
  `GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC` followed by 28 A bases.
  The final eight bases before poly(A), `CATGATGC`, are the region barcode;
  the structured feature is spatial identity, distinct from the droplet cell barcode.
- 🟢 The same list gives the caged O1 strand, second-layer O1′ strand,
  second-layer caged O2 strand, zipcodes 2–4, fluorescent readout oligo,
  lipid-anchor adaptors and the alternative design-2 anchor/zipcode oligos.
- 🟡 Text extraction loses the paper's red highlighting of the precise NPOM-caged
  thymidines. Their positions must be read from the rendered paper before modeling
  caging as residue-level modifications. The 10x droplet capture primers follow
  the kit chemistry rather than an author-specified ordering oligo.
