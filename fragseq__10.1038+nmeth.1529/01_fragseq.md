# FragSeq

**Defining source:** Underwood *et al.*, *Nature Methods* (2010),
doi:10.1038/nmeth.1529. The article and method overview were read.

## Evidence used

- 🟢 Purified nuclear RNA is refolded and partially digested with nuclease P1.
- 🟢 Nuclease P1 preferentially cleaves single-stranded nucleic acid, so mapped fragment
  termini report accessible RNA positions.
- 🟢 Digested and undigested-control fragments are converted with a small-RNA sequencing
  workflow and compared computationally.
- 🟡 Proprietary SOLiD adapter and sequencing-primer sequences are not represented.

## Construct logic

FragSeq is the closest direct RNA analogue of a nuclease-accessibility assay: structure
determines enzyme access, and sequencing locates the resulting molecular ends. Unlike
ATAC-seq, the enzyme cleaves RNA rather than installing adapters in one reaction.
