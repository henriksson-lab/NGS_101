# PAIso-seq — research notes

## Sources read

- 🟢 Liu et al. 2019, DOI `10.1038/s41467-019-13228-9`, PMC full text and Supplementary Table 1.
- 🟢 Liu et al. 2022, DOI `10.1038/s41596-022-00704-8`, detailed protocol metadata.

## Construct path

🟢 A barcoded oligo containing strategically placed deoxyuridines supports end extension over the RNA poly(A) tail. USER cleavage removes that oligo before reverse transcription. Template switching adds the opposite PCR handle, and full-length cDNA is converted to a PacBio SMRTbell. This preserves tail length and internal non-A residues in the same CCS read as the transcript isoform.

🟢 The TSO, RT primer, PCR primer, and all fifteen 2dU-BC extension oligos on the schematic are transcribed from Supplementary Table 1. BC3 is the member explicitly named for the bulk library in the defining paper's Methods; BC1–BC15 barcode individual oocytes.

🟢 **The TSO linker is `ACAT`, published by PAIso-seq itself.** Supplementary Table 1 row
`TSO` reads `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G` verbatim, i.e. the 23-nt SMART handle, the
four-base linker `ACAT` and the `rGrG+G` tail. This matters because SMART-seq2 has two
published linkers — `AC` in Picelli et al. 2013 Supplementary Table 2 for the chosen
`rGrG+G` oligo, `ACAT` in the *Nature Protocols* form — and a protocol that only cites
SMART-seq2 cannot tell them apart. PAIso-seq prints its own, so no base is inherited by
assumption. `tools/chem.py` therefore takes the bases from `smart-seq2__10.1038+nmeth.2639/tools/smartseq2` (handle,
`DT_LINKER` for the RT primer's `AC`, `TSO_LINKER`, `TSO`) rather than retyping them, and
`tools/selftest.py` pins all three ordering strings against Supplementary Table 1 so a
change in the shared constant would fail here rather than silently redraw this page.

🟢 Amplification uses the single published PCR primer. Template switching and reverse transcription place its common handle at opposite ends of the cDNA, so the same physical primer species binds both inward-facing sites. The sample barcode is carried by the extension oligo; it is not a second PCR primer.
