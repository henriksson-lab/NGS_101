# PLAC-seq — research notes

Defining source: Fang *et al.*, [10.1038/cr.2016.137](https://doi.org/10.1038/cr.2016.137).

## Evidence used

- 🟢 PLAC-seq means Proximity Ligation-Assisted ChIP-seq.
- 🟢 Crosslinked nuclei undergo restriction digestion, biotin fill-in and proximity ligation before chromatin shearing.
- 🟢 Immunoprecipitation after ligation enriches contacts anchored at the chosen protein or histone mark; biotin selects contact junctions for paired-end library construction.
- 🟢 An [original PLAC-seq GEO sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2295907) describes end repair, dA addition and Illumina adapter ligation after fragmentation, without naming a custom contact-specific oligonucleotide. The contact junction is created by restriction-end fill-in and proximity ligation, rather than an oligo bridge.
- 🔴 The record does not disclose the lot-specific adapter or index oligo sequences. The final library's standard Illumina shell is a model of the sequencing platform, not a transcription of an exact ordered adapter set from this GEO record.

## Protocol boundary

The order is defining: PLAC-seq ligates in intact nuclei before shearing and ChIP. HiChIP remains a separate named protocol.
