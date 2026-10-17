# NEBNext Ultra II FS DNA PCR-free

## 1. What it is

A PCR-free Illumina DNA workflow combining enzymatic fragmentation/end preparation with full-length, dual-index UMI adaptors.

## 2. Sources and evidence

- 🟢 NEB E7430/E7435 instruction manual **Version 2.0_6/25** (official PDF `manuale7430_e7435.pdf`, read 2026-10-17), and [product page](https://www.neb.com/en-us/products/e7435-nebnext-ultra-ii-fs-dna-pcr-free-library-prep-with-sample-purification-beads). Its revision history lists only 1.0 (9/21) and 2.0 (6/25). An earlier version of this note cited "v3.0, September 2025"; that revision was not located and is no longer cited.
- 🟢 The separate [NEB UDI-UMI Adaptor Sets 1–4 manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7395_e7874_e7876_e7878.pdf), **Version 4.0_6/26**: 96 8-base unique dual index UMI adaptors per set, a 12-base UMI, adaptor plate at 20 µM, compatible with TA single-base-overhang ligation and PCR-free workflows; per-well expected i5/i7 index reads (Section 3); run segments Index 1 (i7) 8 cycles without UMI or 20 with UMI, Index 2 (i5) 8 cycles (p. 23).
- 🔴 Its revision 3.0 (9/25) entry: the adaptor sequence was deleted "because it was inaccurate. The correct adaptor sequence is proprietary." No historical sequence is used here.
- 🟢 Both manuals give the TruSeq-like trimming sequences Read 1 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and Read 2 `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`.

## 3. Molecular path (E7430 Version 2.0_6/25)

High input, 100–500 ng (Section 1):

1. 🟢 Fragmentation/End Prep: 25 µl DNA (in TE), 7 µl FS Reaction Buffer, 1 µl DTT, 2 µl FS Enzyme Mix (35 µl), heated lid 75 °C. 350-bp inserts: 10 min 37 °C; 450-bp inserts: 8 min 37 °C; both then 30 min 65 °C, hold 4 °C.
2. 🟢 Adaptor ligation directly in the FS mix: 2.5 µl NEBNext UMI adaptor, 30 µl Ultra II Ligation Master Mix, 1 µl Ligation Enhancer (68.5 µl); 20 °C, 15 min.
3. 🟢 Cleanup: +40 µl 0.1X TE, 50 µl SPRIselect/Sample Purification Beads, supernatant discarded, two 80% ethanol washes, elution in 102 µl 0.1X TE (100 µl recovered).
4. 🟢 Size selection: 65 µl (350 bp) or 58 µl (450 bp) beads added to the eluate, supernatant discarded, washes, elution in 22 µl (20 µl recovered); qPCR quantitation and Bioanalyzer/TapeStation check. No PCR.

Low input, 50–99 ng (Section 2), 🟢: 26 µl DNA, 1 µl DTT; 20 min 37 °C then 30 min 65 °C for ~300-bp inserts; adaptor diluted 1:10 in UMI Dilution Buffer (2 µM working); one 50 µl bead cleanup after +50 µl 0.1X TE, eluted in 17 µl. Not drawn on the page.

🟡 The split of the single FS reaction into fragmentation, end repair and dA tailing is a model decomposition; the manual names the step "Fragmentation/End Prep" and does not state what the 65 °C incubation does. The A/T ligation junction follows from the adaptor manual's "TA single base overhang ligation".

## 4. Final library and read layout

The full-length adaptor already supplies every sequence needed for cluster generation and sequencing; unlike the USER-hairpin Ultra II page, there is no loop opening and no PCR completion step.

- 🟢 The PCR-free product keeps single-stranded adaptor regions: E7430 Figures 1–2 captions (slower Bioanalyzer migration "due to the single strand regions of the adaptors") and the adaptor manual Figure 1 (forked PCR-free product). The page's final duplex is a fully paired simplification for primer landing and says so.
- 🟢 Index 1 is 20 cycles with UMI (8-base i7 + 12-base UMI); Index 2 is 8 cycles.
- 🟡 The order inside Index 1 (i7 first, then UMI) is drawn but not stated in the inspected text; the low-resolution Figure 1 was not used to assign it.

## 5. Open questions

- Complete adaptor sequence and its fork geometry (proprietary per adaptor manual revision 3.0).
- Order of i7 and UMI within the 20-cycle Index 1 read.
- High-input adaptor working concentration (the E7430 high-input section gives 2.5 µl without a dilution; the adaptor plate is 20 µM).
- Selected adaptor set/well is user-variable.
