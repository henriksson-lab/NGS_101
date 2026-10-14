# SaferSeqS — research notes

**Defining paper:** Cohen *et al.*, “Detection of low-frequency DNA variants by targeted sequencing of the Watson and Crick strands,” *Nature Biotechnology* (2021), DOI [10.1038/s41587-021-00900-z](https://doi.org/10.1038/s41587-021-00900-z).

## Sources read

- 🟢 Defining open author manuscript, Methods and Fig. 1 (PMCID PMC8627329).
- 🟢 Publisher workbook `41587_2021_900_MOESM3_ESM.xlsx`, especially Supplementary Tables 5 and 7. Table 7 is the authority for adapter, nested-PCR and custom sequencing-primer oligos; Table 5 supplies target-specific GSPs.

## Construct-changing steps

1. 🟢 Genomic DNA is dephosphorylated and blunted during end repair.
2. 🟢 A 3′ adapter containing a 14-nt UID is ligated to the two 3′ fragment ends.
3. 🟢 Extension converts the UID into a double-stranded barcode, and ligation of the 5′ adapter completes the library.
4. 🟢 Library PCR creates redundant copies of both original strands.
5. 🟢 Separate strand-specific hemi-nested PCRs enrich Watson and Crick products; a second nested PCR adds sample barcodes and specificity.
6. 🟢 Paired-end sequencing reads the Watson UID from the first 14 bases of Read 1 and the Crick UID from the first 14 bases of Read 2. Custom sequencing primers are mixed with the standard Illumina primers.

## Recovered oligos

🟢 Supplementary Table 7 gives both strands of the 3′ adapter, including the 14-nt
random UID, `/5Phos/`, `/3SpC3/`, deoxyuridines and `/33dT/`; the 5′ adapter; both
library-amplification primers; the Watson and Crick nested-PCR shells; and the custom
Read 1, index and Read 2 sequencing primers. These are now transcribed verbatim on the
schematic page. The only variable bases are the stated eight-base indices and the
target-specific GSP portions listed separately in Supplementary Table 5.

## Interpretation

🟡 “Identical” barcode means that the UID installed at one fragment end is copied into both complementary strands before amplification. Watson and Crick families are separated by oriented target PCR and rejoined computationally through that shared molecular identity.
