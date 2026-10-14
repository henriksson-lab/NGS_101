# OAK — research notes

## Sources read

- 🟢 Darmanis et al. 2024, DOI `10.1038/s41467-024-53227-z`, PMC article, Methods and supplements.

## Construct path

🟢 OAK stands for overloading and unpacking. Fixed cells or nuclei are intentionally overloaded into Chromium droplets. In situ reverse transcription (and, for Multiome, the ATAC route) transfers the droplet bead barcode. The emulsion is broken without destroying the fixed carriers; cells are mixed and redistributed into aliquots, where PCR introduces an aliquot-specific second index. The two identifiers together resolve cells that shared a droplet.

## Published RNA-library oligos

🟢 The defining paper's “OAK scRNA-Seq” Methods print the TSO recognition primer `AAGCAGTGGTATCAACGCAGAGT` and an example aliquot-index primer `AATGATACGGCGACCACCGAGATCTACACAACGTGATACACTCTTTCCCTACACGACGCTCTTCCGATCT` (the eight-base example i5 index is `AACGTGAT`). Subsequent library PCR uses a partial P5 primer `AATGATACGGCGACCACCGAGA` and an example i7 primer `CAAGCAGAAGACGGCATACGAGATCGCATGTTACGTGACTGGAGTTCAGACGTGT` (example i7 index `CGCATGTT`). These are published RNA-route reagents, so a blanket claim that the aliquot primers are unavailable is incorrect.

🔴 The complete 10x bead oligo and the kit-provided Multiome pre-amplification primers and ATAC sample-index oligos are not printed in the paper. The example RNA oligos do not establish every member of the index families or the proprietary ATAC branch.
