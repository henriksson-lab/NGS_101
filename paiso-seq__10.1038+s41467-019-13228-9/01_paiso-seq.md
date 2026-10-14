# PAIso-seq — research notes

## Sources read

- 🟢 Liu et al. 2019, DOI `10.1038/s41467-019-13228-9`, PMC full text and Supplementary Table 1.
- 🟢 Liu et al. 2022, DOI `10.1038/s41596-022-00704-8`, detailed protocol metadata.

## Construct path

🟢 A barcoded oligo containing strategically placed deoxyuridines supports end extension over the RNA poly(A) tail. USER cleavage removes that oligo before reverse transcription. Template switching adds the opposite PCR handle, and full-length cDNA is converted to a PacBio SMRTbell. This preserves tail length and internal non-A residues in the same CCS read as the transcript isoform.

🟢 The TSO, RT primer, PCR primer, and all fifteen 2dU-BC extension oligos on the schematic are transcribed from Supplementary Table 1. BC3 is the member explicitly named for the bulk library in the defining paper's Methods; BC1–BC15 barcode individual oocytes.

🟢 Amplification uses the single published PCR primer. Template switching and reverse transcription place its common handle at opposite ends of the cDNA, so the same physical primer species binds both inward-facing sites. The sample barcode is carried by the extension oligo; it is not a second PCR primer.
