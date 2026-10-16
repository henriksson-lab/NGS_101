# OAK — research notes

## Sources read

- 🟢 Wu et al. 2024, DOI `10.1038/s41467-024-53227-z`, PMC article, “OAK scRNA-Seq” and “OAK paired snRNA-Seq and snATAC-Seq” Methods. The [2024 bioRxiv v1 preprint](https://doi.org/10.1101/2024.01.23.576918) prints the same example RNA primers. The original extended-data figures, Supplementary Fig. 1, reporting summary, peer-review file and figure source-data workbook were checked separately in the source-review checklist.
- 🟢 The original 10x [CG000183 Rev A 3′ v3 guide](https://raw.githubusercontent.com/Teichlab/scg_lib_structs/master/data/10X-Genomics/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf) prints the 16-base bead barcode, **12-base UMI** and poly(dT) architecture. The later [CG000338 Rev F Multiome guide](https://cdn.10xgenomics.com/image/upload/v1666737555/support-documents/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevF.pdf) prints separate RNA and ATAC bead structures and final-library shells. The exact guide revisions used in OAK have not been established.

## Construct path

🟢 OAK stands for overloading and unpacking. In the RNA-only route, fixed cells enter overloaded Chromium droplets before reverse transcription. In the paired Multiome route, fixed nuclei undergo **bulk Tn5 transposition before GEM loading**; subsequent GEM reactions barcode RNA and tagged chromatin. The emulsion is broken without destroying the fixed carriers; cells or nuclei are mixed and redistributed into aliquots, where PCR introduces an aliquot-specific second index. The two identifiers together resolve cells that shared a droplet. The prior graph placed Tn5 after GEM RT; that was a repository modeling error, now corrected.

## Published RNA-library oligos

🟢 The defining paper's “OAK scRNA-Seq” Methods print the TSO recognition primer `AAGCAGTGGTATCAACGCAGAGT` and an example aliquot-index primer `AATGATACGGCGACCACCGAGATCTACACAACGTGATACACTCTTTCCCTACACGACGCTCTTCCGATCT` (the eight-base example i5 index is `AACGTGAT`). Subsequent library PCR uses a partial P5 primer `AATGATACGGCGACCACCGAGA` and an example i7 primer `CAAGCAGAAGACGGCATACGAGATCGCATGTTACGTGACTGGAGTTCAGACGTGT` (example i7 index `CGCATGTT`). These are published RNA-route reagents, so a blanket claim that the aliquot primers are unavailable is incorrect.

🟢 The original 10x guides print generic 3′ v3 and Multiome bead and final-library structures. 🟡 The page uses their 16-base droplet barcode and 12-base RNA UMI and shows the published OAK RNA runs as 28/10/8/90 cycles and ATAC runs as 50/8/24/49 cycles (R1/i7/i5/R2). The RNA-only and paired Multiome branches have separate graph paths. 🔴 The exact OAK-era 10x kit guide revisions and complete branch-specific PCR primer families remain unverified; no missing-information conclusion is accepted from the source bundle alone.
