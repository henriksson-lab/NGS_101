# MARIO — research notes

Defining source: Nguyen *et al.*, [10.1038/ncomms12023](https://doi.org/10.1038/ncomms12023).

## Evidence used

- 🟢 MARIO maps RNA–RNA contacts mediated by a single protein in unperturbed cells.
- 🟢 UV crosslinking preserves RNPs; a biotinylated RNA linker is ligated so a selected chimera has the form RNA1–linker–RNA2.
- 🟢 Streptavidin purification selects linker-containing chimeras and a random barcode supports PCR-deduplication.
- 🟢 The defining paper's Methods section publishes every custom library oligo, rather than merely naming them: the internally biotinylated RNA linker, phosphorothioate-protected linker-complement oligo, activated RT adapter, three RT-primer variants, Cut_oligo, DP5/DP3 and full Illumina PE1/PE2 PCR primers.
- 🟢 The RNA linker is printed as `rCrUrArG/iBiodT/rArGrCrCrCrArUrGrCrArArUrGrCrGrArGrGrA`; the DNA oligo used to remove terminal linker is `T*C*G*C*ATTGCATGGGCTACTAGCAT`.
- 🟢 The pre-adenylated RT adapter is `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. The three 5′-phosphorylated RT primers contain random bases, the source's four-base sample barcode, a BamHI-bearing body, and the adapter-binding end.
- 🟢 The final amplification uses the fully printed Illumina PE PCR Forward Primer 1.0 and Reverse Primer 2.0, followed by paired-end HiSeq sequencing. Standard read-primer binding can therefore be derived from the completed Illumina shell.

## Protocol boundary

MARIO is protein-contact capture, not a psoralen duplex assay; its shared linker is therefore represented separately from PARIS, SPLASH and LIGR-seq.

The oligo table on the schematic preserves the source's internal modifications and all three experimental RT-primer variants. The diagram uses one generic UMI/sample-index-bearing payload, because those variants change identifiers rather than the reaction topology.
