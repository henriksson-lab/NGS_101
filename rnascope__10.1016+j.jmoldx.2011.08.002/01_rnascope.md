# RNAscope

## 1. What it is

An in situ RNA-hybridization assay that uses adjacent double-Z target probes and branched-DNA signal amplification to visualize individual RNA molecules while retaining tissue morphology.

## 2. Sources and evidence

- 🟢 Wang F, Flanagan J, Su N, et al. *RNAscope: A Novel In Situ RNA Analysis Platform for Formalin-Fixed, Paraffin-Embedded Tissues.* Journal of Molecular Diagnostics 14:22–29 (2012). [doi:10.1016/j.jmoldx.2011.08.002](https://doi.org/10.1016/j.jmoldx.2011.08.002).
- 🟢 [Bio-Techne's current technology description](https://www.bio-techne.com/spatial/reagents/rnascope-ish-technology/get-started) confirms 18–25 target-complementary bases and a 14-base tail per Z probe, a 28-base paired-tail site, and 20 probe pairs in a standard set.
- 🟢 The defining paper states 20 amplifier sites per preamplifier and 20 label-probe sites per amplifier: {{= rnascope.MODEL.labels_per_tree }} labels per double-Z pair. It typically targets ~1 kb with 20 pairs (theoretical 8,000 labels) but states "We chose 10 to 20 pairs"; the model therefore uses 10–20 pairs, {{= rnascope.MODEL.labels_per_target[0] }}–{{= rnascope.MODEL.labels_per_target[1] }} theoretical label sites per target RNA (corrected 2026-10-17 from a fixed 20 pairs).
- 🟢 [US7709198B2](https://patents.google.com/patent/US7709198B2/en), the core double-Z patent assigned to Advanced Cell Diagnostics, gives exact human GAPD probe examples in Tables 1–2 and an example preamplifier as SEQ ID NO:1. In its double-Z pair, each probe contains a target-binding region, a five-thymidine spacer and one 14-nt tail; the two tails together pair to the terminal 28 bases of the preamplifier.
- 🟡 The patent is direct evidence for a working double-Z sequence implementation, but not evidence that current commercial RNAscope products retain those exact universal sequences.
- 🔴 The exact oligo pools and amplification-tree sequences in current commercial kits have not been established by this review. The article calls its own target-probe and amplification-tree sequences proprietary; its five supplements and linked vendor material still need independent checks. The schematic uses patent bases only in a separately labelled historical example and does not assign them to current products.

## 3. Assay conditions (Wang 2012 Methods, PMC3338343)

Sections read in the 2026-10-17 bounded pass: RNAscope Design Strategy, Figure 1 legend, RNAscope Assay Procedure for RNA Detection (extracted text lines 1199–1222).

- 🟢 Cultured cells: 4% formaldehyde 60 min; protease 2.5 µg/mL at 23–25 °C.
- 🟢 FFPE: 5-µm sections, xylene deparaffinization, ethanol series, 10 mM citrate pH 6 at 100–103 °C for 15 min, water rinse, 10 µg/mL protease (Sigma-Aldrich) at 40 °C for 30 min in a HybEZ oven. (The text prints "10 nmol/L" citrate; read here as the conventional 10 mM retrieval buffer — 🟡 unit interpretation.)
- 🟢 All hybridizations at 40 °C, each followed by three room-temperature washes in 0.1× SSC, 0.03% lithium dodecyl sulfate (LDS): target probes in buffer A (6× SSC, 25% formamide, 0.2% LDS, blocking reagents) 3 h; preamplifier 2 nM in buffer B (20% formamide, 5× SSC, 0.3% LDS, 10% dextran sulfate, blocking reagents) 30 min; amplifier 2 nM in buffer B 15 min; label probe 2 nM in buffer C (5× SSC, 0.3% LDS, blocking reagents) 15 min.
- 🟢 Labels: Alexa Fluor 488/546/647/750 conjugates (epifluorescence) or HRP conjugate with DAB and hematoxylin counterstain (bright field); the design section also names alkaline phosphatase/Fast Red. Multiplex uses equimolar reagents per amplification system.
- 🟢 Controls: UBC probe set (RNA integrity/assay) and bacterial dapB (background).
- 🔴/open: reagent compositions labelled "blocking reagents", protease identity for cells, and the current-kit pretreatment/amplification steps (e.g., later AMP1–6 / HRP-channel workflows in vendor manuals) are not resolved by the article; vendor manuals Cat. 322310/323100 are not on disk in this bundle (curl error 92) and remain deferred.

## 4. Molecular path

Two target probes bind adjacent RNA regions. Their separate 14-base tails become one stable preamplifier docking site only when juxtaposed. The preamplifier recruits amplifiers, each amplifier recruits labeled probes, and microscopy records the resulting fluorescent or chromogenic dot. This is signal amplification rather than RNA or library amplification.

## 5. Public sequence example

The patent's double-Z GAPD probes use the tails `GAAGTTACCGTTTT` and `CTGAGTCAAAGCAT`. Their concatenation is complementary to the example preamplifier's terminal 28-nt docking site. This makes the cooperative junction reconstructable base by base even though the current commercial reagent sequences are not established.

The patent Table 2 rows were re-read on 2026-10-17 (extracted en.html.txt lines 5497–5545): GAPD217 `ccagtggactccacgacgtacTTT`+`TTgaagttaccgtttt` and GAPD218 `ttctccatggtggtgaagacgTTT`+`TTctgagtcaaagcat`, i.e. 21-nt target arm, T5 spacer, 14-nt tail; SEQ ID NO:1 preamplifier at line 4667 ends in the 28-nt complement of both tails.
