# CUT&RUN

## 1. What it is

CUT&RUN maps a chromatin protein by antibody-tethering Protein A–MNase in intact permeabilized nuclei, activating cleavage with calcium, and sequencing the target-proximal fragments released into solution.

## 2. Sources and evidence

- 🟢 Skene PJ, Henikoff S. *An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites.* eLife 6:e21856 (2017). [doi:10.7554/eLife.21856](https://doi.org/10.7554/eLife.21856).
- 🟢 The defining protocol immobilizes nuclei on concanavalin-A beads, binds antibody and pA–MNase, activates with 2 mM calcium at 0 °C, stops with EDTA/EGTA, and recovers soluble two-cut fragments.
- 🟢 Extracted fragments entered KAPA Illumina library preparation without size selection and were paired-end sequenced.
- 🟡 Canonical TruSeq adapter arms stand in for the named but unprinted KAPA kit oligos.

## 3. Molecular path

CUT&RUN first makes unadapted, antibody-selected DNA fragments. Only after fragment release and extraction are ends repaired, dA-tailed, adapter-ligated and amplified. This separation is the defining contrast with CUT&Tag.
