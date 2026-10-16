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

🟢 The published article specifies at least 120 bases of Read 1 and 95 of Read 2
for genomic and tag information. These are source-stated minima; a fixed
acquired cycle count is not assigned in this model. Biological tag order
is distinct from the author parser's READ2 pattern, Y–Even–Odd–Even–Odd–DPM. Only complete six-tag molecules were
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
strings but prepends a constant {{= len(scsprite.PARSER_PREFIX) }}-base motif to every DPM read string. The
author pipeline snapshot at commit
`2aeb96e7641fd90f1e7c00bab38f05a418977c7c` selects this second configuration
in `config.yaml`; its `Snakefile` forwards that setting to the barcode matcher.
That positively identifies the configuration in this code snapshot, without
assigning it to every historical library. The Java matcher consumes actual
input reads left to right; its matching order does not establish physical
ligation-adapter junctions or sequencing-primer sites.

## Corrected model boundary

The former endpoint helper added a TruSeq adapter shell, fixed 9/17-base tag
widths and computed read-cycle spans without a source-specific scSPRITE
amplicon declaration. Those were repository modelling assumptions. The endpoint
now shows genomic payload plus the logical six-part identifier with unspecified
physical lengths; the page retains the source-reported Read 1/2 minima and leaves
sequencing-primer sites and exact cycle coordinates unresolved. This does not
claim that a physical oligo or primer is absent from published sources.

The earlier audit text called the alternate parser prefix 17 bases. Direct
comparison of the original configurations shows the prefix
`{{= scsprite.PARSER_PREFIX }}` has {{= len(scsprite.PARSER_PREFIX) }} bases.
This is a repository counting correction, not a publication discrepancy.

## Linked author pipeline and layout figure

🟢 The original article's Code Availability links the author repository. At the
pinned commit above, `config.yaml` selects `config_dpm6_y-stag_scSPRITE2.txt`
and trims Read 1 to 100 bases before alignment. `Snakefile` passes the selected
configuration and paired input reads directly to the barcode-identification
program. `BarcodeIdentification.java`, `checkRead` and `checkReadForTag`,
consumes the selected input read left to right while matching the supplied
patterns. These are analysis declarations rather than ordered oligo designs.

🟢 The embedded `images/tag_layout.png` depicts Read 1 pointing through DPM
into DNA and Read 2 pointing toward a chain labeled DPM reverse complement,
Odd, Even, Odd and Y. This is a five-part identifier illustration within the
reused pipeline; the article instead specifies a six-part scSPRITE identifier
and explicitly assigns DPM to Read 1 and the remaining five tags to Read 2
(Methods, “scSPRITE data generation” and “Sequencing analysis pipeline”).
The page follows that explicit scSPRITE method. The generic embedded image
is not used to infer an extra physical adapter or omit a scSPRITE tag. The
figure's historical experimental scope is not established by its location
in the repository.

The README describes this folder as work in progress. Downstream alignment,
cluster, contact-map and cell-selection settings in `config.yaml` operate on
identified reads and do not establish additional library oligos. The
new original pipeline files and figure have two independently attributed
readings in the source review. A fresh independent linked-inventory closure
is still pending; no publication-absence claim is accepted.
