# DiMeLo-seq — research notes

## 1. What it is

DiMeLo-seq records antibody-targeted protein–DNA proximity as exogenous m6A directly on
native chromosomal DNA, then reads those marks on unamplified long molecules.

## 2. Sources read

- 🟢 Altemose et al. 2022, DOI
  [10.1038/s41592-022-01475-6](https://doi.org/10.1038/s41592-022-01475-6), full article
  (PMC9189060).
- 🔴 Maslan et al. 2024 detailed protocol, DOI
  `10.1038/s41596-024-01032-9` (PMC11674881): identified but its full text
  has not been read in this audit. The earlier DOI listed here was incorrect.
- 🟢 Author protocols.io "DiMeLo-seq: directed methylation with long-read sequencing",
  printed as Version 2 (dx.doi.org/10.17504/protocols.io.b2u8qezw; cached as
  `protocolsio_dimelo-seq-…-n2bvjxe4wlk5.pdf`): materials list and steps 10–14 read in
  the 2026-10-17 bounded pass.
- 🔴 Oxford Nanopore SQK-LSK109, SQK-LSK110 and EXP-NBD104/114 documentation of the
  2021 revision cited by the Methods: not read. End geometry and adapter bases stay open.

## 3. Molecular order

🟢 Primary antibody binds the chromatin feature; pA–Hia5 binds the antibody. After washes,
SAM activates Hia5 and nearby adenines acquire m6A. High-molecular-weight DNA is extracted
and sequenced without PCR by a modification-sensitive long-read platform.

🟢 Conditions (Methods, "DiMeLo-seq"): 1–5 M cells are permeabilized in Dig-Wash (0.02%
digitonin, 20 mM HEPES-KOH pH 7.5, 150 mM NaCl, 0.5 mM spermidine, protease inhibitor,
0.1% BSA) for 5 min on ice; primary antibody 1:50 in Tween-Wash (0.1% Tween-20) at 4 °C
for 2 h, two washes; 200 nM pA-Hia5 at room temperature for 1 h, two washes; Activation
Buffer with 800 µM SAM at 37 °C for 30 min. Protocol v2 (LMNB1/CTCF experiments marked in
Supplementary Tables 1 and 3) binds pA-Hia5 at 4 °C for 2 h, activates for 2 h, replenishes
800 µM SAM halfway and lowers spermidine to 0.05 mM. DNA is extracted with NEB Monarch
T3010S or T3050L.

🟢 Library (Methods, "Nanopore library preparation and sequencing"; protocols.io materials
list): 3 µg DNA input; end preparation with the NEBNext End Repair/dA-tailing module
(E7546S) and ligation with the NEBNext Quick Ligation Module (E6056S) inside the ONT kit
workflow. Method (2), SQK-LSK110 (CTCF, H3K9me3 and controls): end prep 1 h plus 30 min
deactivation, SRE-buffer pelleting cleanup, half-volume 1 h ligation, pelleting with two
LFB washes, 500 ng loaded with flow-cell wash and reload every 24 h. Method (1),
SQK-LSK109 with EXP-NBD104/114 (optimization and CENP-A/AlphaHOR-RES): 1 µg end-repaired
DNA per native-barcode ligation, ligations ≥20 min, 3 µg pooled into the final ligation,
LFB cleanup, 1 µg loaded. MinION v9.4 flow cells (FLO-MIN106.1), MinKNOW 21.02.1. The page draws
method (2) and names method (1) and the PacBio branch without drawing them.

## 4. Readout

🟢 Exogenous m6A reports target proximity while endogenous CpG methylation is retained on
the same molecule. The defining paper demonstrates Nanopore readout. 🔴 Kit-specific
motor-adapter bases are not established by the sources checked so far; the page states pore entry
without inventing a primer or sequence.
🟢 The defining paper names different nanopore library reagents by experiment: SQK-LSK109
with EXP-NBD104/114 native barcodes for optimization and CENP-A experiments, and
SQK-LSK110 without that barcode step for the stated CTCF/H3K9me3 experiments. It also
reports a separate PacBio SMRTbell/Sequel IIe experiment with a vendor binding primer.
Those kit adapters and the PacBio primer are distinct, version-specific roles; the
nanopore schematic does not assert an exact sequence for either platform.
