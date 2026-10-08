# Digital Telomere Measurement (DTM)

## 1. What it is

DTM ligates barcoded capture/tether duplexes directly to native human telomeres, fills the terminal gap, fragments the tagged DNA and measures individual telomeres by Oxford Nanopore sequencing.

## 2. Sources read

- 🟢 Sanchez et al. (2024), [doi:10.1038/s41467-024-49007-4](https://doi.org/10.1038/s41467-024-49007-4), full text and methods.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), workflow comparison.

## 3. Construct and reaction order

1. 🟢 Equimolar barcoded capture oligo and sequencing tether are duplexed.
2. 🟢 The fresh duplex capture oligos are ligated to high-molecular-weight DNA overnight.
3. 🟢 Sulfolobus DNA polymerase IV fills potential gaps between capture oligo and the ds/ss telomere junction.
4. 🟢 Tagged DNA is mechanically fragmented and converted to an Oxford Nanopore library.

The locally accessible article body establishes the architecture and reaction order but not the ordered oligo bases. Those bases remain explicit unknowns pending verification of the supplementary oligo table.

## 4. Final library and sequencing

A motor-loaded nanopore adapter feeds an intact tagged strand through the pore; no sequencing primer binds.
