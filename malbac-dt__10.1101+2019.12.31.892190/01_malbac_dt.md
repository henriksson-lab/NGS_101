# MALBAC-DT — research notes

Defining source: Chen *et al.*, [10.1101/2019.12.31.892190](https://doi.org/10.1101/2019.12.31.892190); published protocol details are also in [PNAS 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9907105/).

## Evidence used

- 🟢 A poly(dT) RT primer installs a well barcode and UMI before amplification.
- 🟢 Random primers drive MALBAC quasi-linear amplification, followed by additional PCR cycles.
- 🟢 Wells are pooled for TruePrep tagmentation; Read 1 reports barcode/UMI and Read 2 reports transcript sequence.
- 🟢 The [published MALBAC-DT Methods](https://doi.org/10.1073/pnas.2206938119) name `RT-An` (reverse transcription), `RT-Bn` (added during Exo I treatment), `GAT-7N` and `GAT-COM` (quasi-linear amplification), `Tru2-Gn-RT` (last five amplification cycles), and `Ix-Tru2` (substituted for a Nextera S7XX primer during library PCR). The original plate uses eight row RT primers and twelve column finishing primers; the automated plate uses 96 well-specific RT primers and one finishing primer.
- 🔴 The full sequences and modifications of those named MALBAC-DT oligos have not been verified from the preprint or published supporting information here. Their names and reaction roles must not be used to invent a complete primer list or barcode placement.

## Protocol boundary

MALBAC-DT applies MALBAC-style amplification to digital transcriptomics and is separate from the MALBAC single-cell DNA protocol.
