# Digital Telomere Measurement (DTM)

## 1. What it is

DTM ligates barcoded capture/tether duplexes directly to native human telomeres, fills the terminal gap and measures individual telomeres by Oxford Nanopore sequencing.

## 2. Sources read

- 🟢 Sanchez et al. (2024), [doi:10.1038/s41467-024-49007-4](https://doi.org/10.1038/s41467-024-49007-4), full text and methods.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), workflow comparison.

## 3. Construct and reaction order

1. 🟢 Equimolar barcoded capture oligo and sequencing tether are duplexed.
2. 🟢 The fresh duplex capture oligos are ligated to high-molecular-weight DNA overnight.
3. 🟢 Sulfolobus DNA polymerase IV fills potential gaps between capture oligo and the ds/ss telomere junction.
4. 🟢 An Oxford Nanopore sequencing adapter is ligated to the prepared DNA.

🟢 Supplementary Data 3 gives `seqTether` as `AACCTTGGAGATGCACGGAGCAAGCAAT` and the representative 5′-phosphorylated p-t3-nb01 capture oligo as `TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC`. The latter comprises a 22-nt tether-complementary arm, a 25-nt barcode and a terminal 6-nt `CTAACC` telomere arm. Eight canonical-telomere barcode variants and seven telomere-sequence-variant capture oligos are provided; the page uses p-t3-nb01 as the representative molecule.

🟢 The core Methods used high-molecular-weight DNA directly. The paper discusses later/optional restriction-digest refinements, but the page no longer invents a mechanical-fragmentation step.

## 4. Final library and sequencing

A motor-loaded nanopore adapter feeds an intact tagged strand through the pore; no sequencing primer binds.
