# GoT — research notes

Defining source: Nam *et al.*, *Nature* (2019), [10.1038/s41586-019-1367-0](https://doi.org/10.1038/s41586-019-1367-0) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6782071/)).

## 1. What it is

Genotyping of Transcriptomes (GoT) adds targeted genotyping to a 10x single-cell transcriptome by making an amplicon that reaches from a chosen transcribed variant to the already attached 10x cell barcode and UMI.

## Sources read

- 🟢 The defining paper's “Genotyping of Transcriptomes (GoT)” and “Circularization GoT” Methods, run configurations and supplementary-primer descriptions were read.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1038/s41586-019-1367-0`; the open manuscript is the readable primary source.

## Protocol variants kept explicit

- 🟢 Standard 10x Chromium 3′ v2, 3′ v3 or 5′ chemistry runs through emulsion break and first-strand-cDNA recovery.
- 🟢 For weakly detected 3′ targets, a gene-specific primer may be spiked into cDNA amplification. For 5′ libraries, a gene-specific RT primer and later additive primer improve recovery.
- 🟢 A small aliquot of amplified cDNA enters GoT; the remainder continues to the standard expression library.

## Construct-changing steps

1. 🟢 Target PCR starts about 50 bases upstream of the variant and amplifies toward the 10x end, thereby preserving cell barcode and UMI.
2. 🟢 Further amplification/indexing makes the genotype library. It may be sequenced alone or spiked into the expression library.
3. 🟢 “Circularization GoT” is an explicitly named alternative amplification strategy in the same paper; it is noted but not silently substituted in the main linear schematic.

## Final library and read layout

- 🟢 For 3′ v3: Read 1 is 28 cycles, Read 2 is 98 or 130 cycles, and the sample-index read is 8 cycles. The exact barcode/UMI split follows the chosen 10x chemistry.
- 🟢 The cell barcode and UMI are the molecular join key to the parent expression profile.
- 🟡 Exact 10x primer bases are not retyped from a secondary source; the final page records the barcode, UMI, target and read direction by role.
