# DART-seq — research notes

Defining source: Meyer, [10.1038/s41592-019-0570-0](https://doi.org/10.1038/s41592-019-0570-0).

## Evidence used

- 🟢 DART-seq means deamination adjacent to RNA modification targets.
- 🟢 APOBEC1 fused to the m6A-binding YTH domain deaminates cytidines near m6A in cellular RNA.
- 🟢 Reverse transcription and ordinary RNA-seq report the resulting C-to-U changes as C-to-T mismatches.
- 🟢 The [defining article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884681/) states that the edited RNA is detected by standard RNA-seq; the APOBEC1–YTH editing chemistry itself installs no assay-specific library oligo. Exact RT, adapter and index oligos depend on the chosen conventional RNA-seq preparation and are not fixed by the DART-seq assay name.

## Protocol boundary

The informative chemistry occurs in the RNA before library construction; the page therefore shows both the molecular mark-writing step and its conventional sequencing endpoint.

This is Meyer's m6A-editing DART-seq, distinct from the unrelated droplet-assisted RNA-targeting method also called DART-seq.
