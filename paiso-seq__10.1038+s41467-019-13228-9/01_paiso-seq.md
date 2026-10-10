# PAIso-seq — research notes

## Sources read

- 🟢 Liu et al. 2019, DOI `10.1038/s41467-019-13228-9`, PMC full text and Supplementary Table 1.
- 🟢 Liu et al. 2022, DOI `10.1038/s41596-022-00704-8`, detailed protocol metadata.

## Construct path

🟢 A barcoded oligo containing strategically placed deoxyuridines supports end extension over the RNA poly(A) tail. USER cleavage removes that oligo before reverse transcription. Template switching adds the opposite PCR handle, and full-length cDNA is converted to a PacBio SMRTbell. This preserves tail length and internal non-A residues in the same CCS read as the transcript isoform.

🟢 The TSO, RT primer and PCR primer on the schematic are transcribed from Supplementary Table 1. The 2dU barcode family remains in the research record rather than presenting one barcode as universal.
