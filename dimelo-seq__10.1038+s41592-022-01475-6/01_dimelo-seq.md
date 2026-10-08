# DiMeLo-seq — research notes

## 1. What it is

DiMeLo-seq records antibody-targeted protein–DNA proximity as exogenous m6A directly on
native chromosomal DNA, then reads those marks on unamplified long molecules.

## 2. Sources read

- 🟢 Altemose et al. 2022, DOI
  [10.1038/s41592-022-01475-6](https://doi.org/10.1038/s41592-022-01475-6), full article
  (PMC9189060).
- 🟢 Altemose et al. 2024 detailed protocol, DOI `10.1038/s41596-023-00969-3`, full text
  (PMC11674881), used only to confirm the physical reaction order.

## 3. Molecular order

🟢 Primary antibody binds the chromatin feature; pA–Hia5 binds the antibody. After washes,
SAM activates Hia5 and nearby adenines acquire m6A. High-molecular-weight DNA is extracted
and sequenced without PCR by a modification-sensitive long-read platform.

## 4. Readout

🟢 Exogenous m6A reports target proximity while endogenous CpG methylation is retained on
the same molecule. The defining paper demonstrates Nanopore readout. 🔴 Kit-specific
motor-adapter bases are proprietary/version-dependent, so the page states pore entry
without inventing a primer or sequence.
