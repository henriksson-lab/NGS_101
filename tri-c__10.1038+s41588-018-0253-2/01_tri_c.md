# Tri-C — research notes

## Citation

Oudelaar AM et al. *Single-allele chromatin interactions identify regulatory hubs in dynamic compartmentalized domains.* Nature Genetics 50, 1744–1751 (2018). [doi:10.1038/s41588-018-0253-2](https://doi.org/10.1038/s41588-018-0253-2)

## Sources read

- 🟢 The defining author manuscript and Methods in [PMC6265079](https://pmc.ncbi.nlm.nih.gov/articles/PMC6265079/).
- 🟢 The later Capture-C protocol's Tri-C module in [PMC7613269](https://pmc.ncbi.nlm.nih.gov/articles/PMC7613269/), used as corroboration.

## 1. What it is

🟢 Tri-C targets multiway rather than only pairwise chromatin contacts. A single physical 3C concatemer must retain the viewpoint and at least two reporter fragments, so its two proximity-ligation junctions represent a concurrent single-allele contact.

## 2. Construct-changing steps

1. 🟢 Fix cells in 2% formaldehyde, digest with NlaIII and proximity-ligate in nuclei.
2. 🟢 NlaIII was selected because it produces small fragments; the defining experiment reports a median mouse-genome fragment size of 131 bp.
3. 🟢 Reverse crosslinks and sonicate to about 450 bp. About half of these fragments can contain multiple ligation junctions.
4. 🟢 Install Illumina TruSeq adapters with NEBNext DNA Library Prep reagents and perform seven cycles of indexing PCR.
5. 🟢 Perform double oligonucleotide capture with a biotinylated probe targeted to the middle of the short viewpoint restriction fragment.

## 3. Final library and reads

🟢 The defining experiments used 250-bp or 150-bp paired-end Illumina sequencing. The molecular unit is viewpoint–reporter 1–reporter 2 (or more), not three independently observed pairwise reads.

🟢 Random sonication ends serve as unique molecular identifiers for PCR duplicate removal.

🟡 Capture-probe bases are locus-specific and therefore remain role tokens on the generic page.
