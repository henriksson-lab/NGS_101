# GoT-ChA

**Defining preprint:** Namikawa et al. (2022), doi:10.1101/2022.05.11.491515. A later detailed protocol, *Mapping Genotypes to Chromatin Accessibility Profiles in Single Cells*, was also used for reaction order and sequencing layout.

## Molecular boundary

🟢 A locus-specific primer is spiked into the single-cell ATAC droplet barcoding reaction, copying the scATAC cell barcode onto an amplicon spanning the genomic variant.

🟢 An aliquot undergoes biotinylated hemi-nested PCR, streptavidin capture and on-bead indexing. The resulting genotype library can be sequenced with the parent accessibility library and linked by cell barcode.

## Evidence boundary

🟢 Three locus-specific primers are designed per target; their genomic binding sequences necessarily change with the variant. The conserved partial Nextera/P5/IS2 handle roles are shown without promoting one locus panel into a protocol constant.
