# CUTAC — research notes

## 1. What it is

CUTAC uses antibody-tethered pA/G–Tn5 under low-salt conditions to preferentially tag
accessible DNA adjacent to transcription-associated nucleosomes.

## 2. Sources read

- 🟢 Henikoff et al. 2020, DOI
  [10.7554/eLife.63274](https://doi.org/10.7554/eLife.63274), full article
  and Methods (PMC7721439 / PMID 33191916).
- 🟢 The Methods name commercial Epicypher pAG–Tn5 and indexed P5/P7 PCR primers, with gap filling and 12 PCR cycles after tagmentation. They do not print the loaded commercial transposome's adapter oligos or the selected index sequences.

## 3. Molecular order

🟢 Antibody binds H3K4me2, H3K4me3 or another transcription-associated target in
permeabilized cells. Adapter-loaded pA/G–Tn5 binds the antibody. Low ionic strength during
magnesium activation favors tagmentation of the neighboring nucleosome-depleted region.
DNA is released and indexed PCR completes the Nextera library.

## 4. Interpretation boundary

🟢 The paper shows that low-salt products coincide with accessible sites and are shorter
than ordinary nucleosomal CUT&Tag fragments. 🟡 The page draws the antibody-to-NDR reach
as spatial geometry only; it is not a molecular bond or a precisely measured distance.
🟡 The displayed Nextera arms are a conventional Tn5-compatible library interpretation;
the 2020 paper alone does not establish the ordered oligos in the Epicypher pAG–Tn5 reagent.
