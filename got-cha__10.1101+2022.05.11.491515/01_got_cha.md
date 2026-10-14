# GoT-ChA

**Defining preprint:** Namikawa et al. (2022), doi:10.1101/2022.05.11.491515. A later detailed protocol, *Mapping Genotypes to Chromatin Accessibility Profiles in Single Cells*, was also used for reaction order and sequencing layout.

## Molecular boundary

🟢 A locus-specific primer is spiked into the single-cell ATAC droplet barcoding reaction, copying the scATAC cell barcode onto an amplicon spanning the genomic variant.

🟢 An aliquot undergoes biotinylated hemi-nested PCR, streptavidin capture and on-bead indexing. The resulting genotype library can be sequenced with the parent accessibility library and linked by cell barcode.

## Evidence boundary

🟢 Three locus-specific primers are designed per target; their genomic binding sequences necessarily change with the variant. The conserved partial Nextera/P5/IS2 handle roles are shown without promoting one locus panel into a protocol constant.

## Oligo source audit

🟢 The defining bioRxiv PDF (`2022.05.11.491515_v1.full.pdf.txt`, Methods “GoT-ChA library construction”) prints the conserved primer handles. `GoT-ChA_R1N` is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-[22bp locus specific]`; `GoT-ChA_IS2` is `AGCAAGTGAGAAGCATCGTGTC-[22bp locus specific]`. The hemi-nested reaction uses P5 `AATGATACGGCGACCACCGAGATCTACAC` and a biotinylated `GoTChA_nested` primer, `/5BiosG/CCTTGGCACCCGAGAATTCCA-[22bp locus specific]`. The on-bead `RPI-X` indexing primer is `CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`, with `XXXXXXXX` defined as the chosen sample index.

🟢 Extended Data Fig. 1b in the same PDF prints the demonstrated TP53 R248Q and JAK2 V617F genomic binding sites as image text; Fig. 1c diagrams the two amplification steps and final library. These locus sequences can be transcribed for a *named example* after visual cross-check, but they are not a universal GoT-ChA primer set.

🟡 The source gives enough invariant handle sequence to model the fixed junctions and custom library read geometry. Its target-specific bases and user-selected i7 index must remain explicit variables unless one demonstrated panel is selected.
