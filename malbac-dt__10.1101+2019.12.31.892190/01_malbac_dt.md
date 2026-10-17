# MALBAC-DT — research notes

Authority (published version first): Chapman, Lee, Cai, Ma, Li, Sun, Xie, "Correlated gene
modules uncovered by high-precision single-cell transcriptomics", *PNAS* 119, e2206938119
(2022), [10.1073/pnas.2206938119](https://doi.org/10.1073/pnas.2206938119),
[PMC9907105](https://pmc.ncbi.nlm.nih.gov/articles/PMC9907105/). The directory is named for
the superseded preprint (Chapman, Lee *et al.*, bioRxiv
[10.1101/2019.12.31.892190](https://doi.org/10.1101/2019.12.31.892190), v1 posted 2 Jan 2020);
an earlier version of these notes misattributed it to "Chen *et al.*".

## Sources read (2026-10-17 bounded pass)

- PNAS JATS (PMC Cloud) Materials and Methods "MALBAC-DT Protocol" (`.xml.txt` lines 195–247),
  "Automated MALBAC-DT Protocol" (lines 249–251), "Sequence Processing" (lines 257–259); Results
  "Development of MALBAC-DT" (lines 115–121); Fig. 1D legend (line 109).
- PNAS SI Appendix (`supp_pnas.2206938119.sapp.pdf.txt`): Figs S1–S10 and Tables S1–S3 only
  (performance and module tables, lines 165–276); no methods or oligo table.
- Preprint v1 Methods (lines 697–871) only to confirm supersession. The preprint says Ix-Tru2
  replaces the Nextera **S5XX** primer; the published text says **S7XX**. This is an ordinary
  revision change and the published value is used.
- `tools/scrape_primers.py` over the whole bundle: 0 sequences.

## Evidence

- 🟢 Cells are sorted into 3 µL lysis buffer (SSIV buffer, ICA-630, betaine, SUPERase In,
  0.2 µL 50 µM RT-An primer, dNTP); 72 °C 3 min, then SuperScript IV RT at 55 °C 10 min.
- 🟢 Fig. 1D: the RT primer is poly-T with a cell-specific barcode and UMI. Sequence
  Processing: UMIs follow the (HBDV)<sub>5</sub> RT-An pattern (20 nt); consensus UMIs must
  differ from the (VDBH)<sub>5</sub> RT-Bn pattern at ≥3 bases.
- 🟢 Excess RT primers are degraded with Exo I (37 °C 30 min, 80 °C 20 min); the Exo mix also
  contains 0.2 µL 50 µM RT-Bn. 🔴 RT-Bn's molecular role is not stated.
- 🟢 Amplification: GAT-7N + GAT-COM, Deep Vent (exo-) in ThermoPol; 10 MALBAC looping cycles
  (4→65 °C ramp, 95 °C, 58 °C), then PCR steps 12–16 written "Goto 13 17×"; then Tru2-Gn-RT and
  five more cycles of steps 12–15.
- 🟢 Original plate: 8 row RT3-An primers and 12 column Tru2-Gn-RT primers (inner/outer
  barcodes demultiplexed per cell). Automated plate (K562): 96 well-specific RT3-An and one
  Tru2-G-RT; 15 amplification cycles + 5.
- 🟢 1 µL per well is pooled and purified with 0.8× AMPure; the Nextera kit tagments, and Ix-Tru2
  primers replace the Nextera S7XX primer to select the barcoded 3′ ends. The automated
  plate used Vazyme TruePrep and NovaSeq 2 × 150 bp with a custom Read 2 primer.
- 🔴 No sequence or modification is printed for RT-An, RT-Bn, GAT-7N, GAT-COM, Tru2-Gn-RT or
  Ix-Tru2 in the article or SI Appendix. Do not invent them from the names.

## Model (derived)

- 🟡 Final library: Nextera P5–i5–S5–ME, transcript cDNA, then the barcoded 3′ end on the
  Ix-Tru2 side. Barcode/UMI therefore lie on the Read 2 side (custom Read 2 primer); Read 1
  reads transcript. This corrects the previous model, which put the barcode/UMI after the
  Read 1 mosaic end and drew a symmetric Nextera product.
- 🟡 The order poly(T)/UMI/inner barcode/RT handle inside RT-An, Tru2–Gn–RT order (from the
  primer name) and the P7 end are inferred and drawn as `<inf>` role tokens.

## Open questions

- Oligo sequences and barcode lengths; the cited automated protocol (reference 33 in the
  Methods resolves to an unrelated p53 paper in the reference list).
- Whether "Goto 13 17×" means 17 or 18 PCR cycles; read lengths of the original-plate run.

## Protocol boundary

MALBAC-DT applies MALBAC-style amplification to digital transcriptomics and is separate from
the MALBAC single-cell DNA protocol.
