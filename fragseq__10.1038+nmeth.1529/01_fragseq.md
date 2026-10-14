# FragSeq

**Defining source:** Underwood *et al.*, *Nature Methods* (2010),
doi:10.1038/nmeth.1529. The article and method overview were read.

## Evidence used

- 🟢 Purified nuclear RNA is refolded and partially digested with nuclease P1.
- 🟢 Nuclease P1 preferentially cleaves single-stranded nucleic acid, so mapped fragment
  termini report accessible RNA positions.
- 🟢 Digested and undigested-control fragments are converted with a small-RNA sequencing
  workflow and compared computationally.
- 🟢 The supplementary library methods specify ABI SOLiD Small RNA Expression Kit adapters ligated to fragments with 5′ phosphate and 3′ hydroxyl, then one sample-specific 3′ PCR primer (Supplementary Methods, “Library generation and high-throughput sequencing”). Supplementary Data File 4 is the authors' actual read-filter FASTA. Its named `P1` pattern is `NNNNNCCACTACGCCTCCGTTTCCTCTCTATGGGCAGTCGGTGATNNNNN`; `adaptor_bc1_P2` is `NNNNNGCGGAACCGGCATGTCGTCCCCGAATCTCTTACTCCTTGGGCCCCGTCNNNNN`. It also lists 15 other six-base barcode variants. The flanking `N` characters are filter-context wildcards, not claimed synthesized bases.
- 🟡 Those filter patterns establish adapter-associated sequence used by the read-processing pipeline, but do not establish the complete ordered, modified ligation oligos or the SOLiD sequencing primer. The schematic's unresolved oligo roles remain warranted.

## Construct logic

FragSeq is the closest direct RNA analogue of a nuclease-accessibility assay: structure
determines enzyme access, and sequencing locates the resulting molecular ends. Unlike
ATAC-seq, the enzyme cleaves RNA rather than installing adapters in one reaction.
