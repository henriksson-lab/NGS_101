# MALBAC-DT — research notes

Defining source: Chen *et al.*, [10.1101/2019.12.31.892190](https://doi.org/10.1101/2019.12.31.892190); published protocol details are also in [PNAS 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9907105/).

## Evidence used

- 🟢 A poly(dT) RT primer installs a well barcode and UMI before amplification.
- 🟢 Random primers drive MALBAC quasi-linear amplification, followed by additional PCR cycles.
- 🟢 Wells are pooled for TruePrep tagmentation; Read 1 reports barcode/UMI and Read 2 reports transcript sequence.

## Protocol boundary

MALBAC-DT applies MALBAC-style amplification to digital transcriptomics and is separate from the MALBAC single-cell DNA protocol.
