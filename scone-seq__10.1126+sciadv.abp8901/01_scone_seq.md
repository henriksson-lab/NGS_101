# scONE-seq

**Defining paper:** Yu et al., *Science Advances* (2023), doi:10.1126/sciadv.abp8901. The open article and supplement were read from PMCID PMC9812385.

## Molecular boundary

🟢 Custom Tn5 adapters label genomic fragments with a six-base DNA-branch barcode and eight-base UMI. Internally primed RT labels total-RNA cDNA with a distinct six-base RNA-branch barcode and eight-base UMI.

🟢 Poly(C) tailing installs the opposite cDNA amplification end. DNA and RNA products are then amplified together in one tube and separated after sequencing by their branch identifiers.

## Evidence boundary

🟢 The source defines the exact adapter and RT oligos in supplementary tables. The schematic records their invariant barcode lengths and molecular roles; individual indexing choices remain in the linked research source.
