# SiC-seq — hydrogel-retained genomes with droplet barcoding

Lan et al. (2017), *Single-cell genome sequencing at ultra-high-throughput with microfluidic droplet barcoding*, [10.1038/nbt.3880](https://doi.org/10.1038/nbt.3880).

## Source scope

The defining published-version author manuscript was read as
`nbt.3880_PMC5531050.xml.txt`, Results workflow (39–43), Online Methods barcode
emulsions, microgel handling/lysis, tagmentation and barcoding (107–131).
The author-manuscript supplement Figure 1/2 captions were read; the publisher
MOESM21 PDF correction notice identifies the corrected 21 June 2017 revision
(Figure 11a taxonomy display). The molecular workflow does not rely on that
taxonomy panel. Complete rendered-figure review remains deferred. The source-bundle oligo scraper was run and relevant sequences verified.
This is a capped pass, not a complete missing-information audit.

## Physical retention and molecular operations

🟢 Microbial cells are encapsulated in agarose microgels. Following transfer to
aqueous buffer, cell-wall digestion and detergent/protease treatment release
high-molecular-weight genomes while the gel retains them. Washes physically
remove soluble material before tagmentation. The gel is re-encapsulated during
Nextera processing; transposase-associated fragments remain a macromolecular
complex. This is a source-declared qualitative retention mechanism, not a
computed diffusion or perfect genome-retention prediction.

🟢 A separate emulsion amplifies single random 15-base barcode templates.
Tagmented-gel droplets, barcode droplets and PCR reagents are merged. Barcodes
are spliced to genomic fragments via transposase-added handles. After PCR, large
coalesced droplets are **physically removed**, then the remaining emulsion is
broken, DNA column-purified and 200–600-bp products selected with AMPure.
The source synthetic-community run is single-end MiSeq 150 bases with a custom
index primer; exact final molecule/primer-cycle mapping remains unassigned.
The barcode is a cell identifier, not a UMI or sample-sheet index.

## Source oligos

🟢 Barcode template:
`GCAGCTGGCGTAATAGCGAGTACAATCTGCTCTGATGCCGCATAGNNNNNNNNNNNNNNNTAAGCCAGCCCCGACACT`.
FL127: `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTC`.
FL128: `CTGTCTCTTATACACATCTCCGAGCCCACGAGACGTGTCGGGGCTGGCTTA`.
FL129: `CAAGCAGAAGACGGCATACGAGATCAGCTGGCGTAATAGCG`.
Custom index primer: `GCCCACGAGACGTGTCGGGGCTGGCTTA`.
Canonical Illumina/Nextera regions are imported by the local module.

The source uses Nextera DNA Library Prep **FC-121-1030** for tagmentation, and
Nextera XT **FC-131-1024 NT buffer** later during barcoding. These are distinct
reagent roles, not an interchangeable modern tagmentation kit. Full loaded-end
geometry and outcome classes remain unassigned rather than borrowed from a
similar protocol. The graph begins with source-prepared cell-containing gels;
marine collection and gel fabrication are upstream scope items.

## Lineage and deferred originals

The article explicitly uses the barcode-droplet strategy attributed to Lan et al.
(2016), 10.1038/ncomms11784. Amini/CPT-seq is cited for the transposase-retention
mechanism, not as evidence that nuclei combinatorial indexing is used here.
These are component-scoped links requiring predecessor-source qualification.
Deferred: rendered supplement figures, directly cited barcoding/retention
originals, historical kit design, exact splicing/product/read geometry, other
sample branches and independent review. No sequence or reaction is confirmed
absent. Status remains notes.
