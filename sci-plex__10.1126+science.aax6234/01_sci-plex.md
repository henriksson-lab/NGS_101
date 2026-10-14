# sci-Plex — source notes

Defining paper: Srivatsan et al., *Science* (2020), [doi:10.1126/science.aax6234](https://doi.org/10.1126/science.aax6234).

- 🟢 PMC full text and supplementary methods were read.
- 🟢 The authors' [sci-Plex reagent repository](https://github.com/cole-trapnell-lab/sci-plex/tree/master/reagents/IDT_hash_oligo_order_forms)
  publishes twelve 96-well IDT order forms (1,152 entries, all downloaded and
  checked for the common 77-character ordering pattern). Plate 1, well A1 (`sciPlex_0`)
  prints the Nextera Read 2 handle, a 10-base well barcode (`TTCTCGCATG`),
  literal `B`, and 32 A bases in that order. The representative complete
  oligo is now copyable on the page. `B` is preserved as printed because the
  order form does not explain whether it denotes a mixed base or another
  ordering convention.
- 🟢 Permeabilized nuclei are labeled with well-specific, unmodified, polyadenylated ssDNA hash oligos and fixed before pooling.
- 🟢 Hash oligos and endogenous poly(A) RNA are indexed by the same sci-RNA-seq rounds, linking sample/perturbation identity to the transcriptome.
- 🟢 The large screen used pairs of hash oligos per treatment population.
- 🟡 The page focuses on the new nuclear-hashing front end; the downstream library endpoint follows the repository's separate sci-RNA-seq model.
- 🟡 The page displays a representative hash only. The other 1,151
  author-published ordering rows can be retrieved from the linked workbooks;
  the paper's supplemental methods do not require inferring them from reads.
