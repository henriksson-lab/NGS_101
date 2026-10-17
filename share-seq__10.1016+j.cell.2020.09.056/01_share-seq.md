# SHARE-seq — research notes

## 1. What it is

SHARE-seq jointly profiles accessible chromatin and RNA in fixed single cells using the same three-round split-pool barcode combination for both molecular classes.

## 2. Sources read

- 🟢 Ma et al., *Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin*, Cell 183, 1103–1116.e20 (2020), DOI [10.1016/j.cell.2020.09.056](https://doi.org/10.1016/j.cell.2020.09.056), open manuscript [PMC7669735](https://pmc.ncbi.nlm.nih.gov/articles/PMC7669735/).
- 🟢 The paper's Methods and Figure S1 establish the reaction order, three 96-well hybridization/ligation rounds, biotin separation and two final libraries.
- Later author protocol (not the defining revision; not read in the 2026-10-17 pass): [protocols.io SHARE-seq v2.2](https://www.protocols.io/view/share-seq-protocol-v2-2-g2y6byfzf.pdf).
- 🟢 Supplementary Table S1 gives the RT primer, template-switching oligo, RNA PCR primer, all three barcode-oligo constant arms, all three linkers and the blocking oligos. The model's sequences were re-checked against the publisher mmc1 workbook (bulk oligos, Ad1, Plate_R1–R3 sheets) on 2026-10-17; agreement with the scg_lib_structs drawing is not used as evidence.

## 3. Construct-changing steps

1. Tn5-tagment accessible chromatin in fixed permeabilized cells. The ligatable Tn5 product and the RT primer expose the same 15-nt `GTCTCGTGGGCTCGG` arm.
2. Reverse-transcribe RNA with the 5′-phosphorylated primer carrying a Nextera-derived handle, 10-nt UMI, internal biotin-dT and oligo-dT/VN tract.
3. Assemble the barcode chain through three 96-well hybridization rounds. After each round, a complementary blocking oligo sequesters free adapters before cells are pooled and redistributed. Each linker is exactly the reverse complement of the incoming molecule arm followed by the reverse complement of that round's barcode-oligo 3′ arm; the model constructs and validates this geometry rather than storing it as a caption.
4. After the third hybridization/blocking round, add T4 DNA ligase once to seal all three backbone nicks.
5. Reverse crosslinks and capture biotinylated cDNA on streptavidin beads.
6. Construct separate RNA and ATAC sequencing libraries whose shared barcode triplet identifies the same cell.

## 3a. Library amplification and read layout (2026-10-17 bounded pass)

🟢 Published JATS (PMC7669735) STAR Methods, lines 241–249:

- scATAC-seq: streptavidin supernatant → MinElute → PCR with NEBNext, 0.5 µM
  "library-specific Ad1 primer" and 0.5 µM P7 primer (72 °C 5 min gap fill, 5 cycles,
  then qPCR-determined extra cycles).
- cDNA: on-bead template switching (2.5 µM TSO, Maxima H Minus, 42 °C 90 min), then
  KAPA HiFi PCR with 400 nM P7 primer and 400 nM RNA PCR primer (5 cycles + qPCR),
  0.8× (cell lines) / 0.6× (tissue) AMPure.
- scRNA-seq: Tn5 loaded with the Read1 oligo only fragments 50 ng cDNA (55 °C 5 min);
  Zymo cleanup; NEBNext PCR with P7 and "Ad1 primer with sample barcodes", 7 cycles;
  0.7× AMPure.
- Sequencing: NextSeq 150-cycle High-Output (Read 1 30, Index 1 99, Index 2 8, Read 2 30)
  or NovaSeq S1 200-cycle (50/99/8/50). "For each read, there are four sets of barcodes
  (eight bases each) in the indexing reads"; RNA reads must carry TTTTTT at Read 2
  bases 11–16.

🟢 Table S1: `P7` = `CAAGCAGAAGACGGCATACGAGAT`, which is also the 5′ arm of every round 3
barcode oligo; sheet `Ad1` prints 96 primers `P5–[8-nt]–TCGTCGGCAGCGTCAGATGTGTAT`.

🟡 Derived layout (computed by the model): Index 1 primed from the mosaic-end/S7 anchor
reads round 1, 2 and 3 barcodes with their 15-nt arms in exactly
{{= sum(len(r.barcode_3_arm) + r.barcode_length + len(r.barcode_5_arm) for r in share_seq.ROUNDS) - len(share_seq.ROUNDS[-1].barcode_5_arm) }}
cycles, matching the published 99-cycle Index 1; Index 2 reads the Ad1 barcode; Read 2 of
the RNA library begins with the 10-nt UMI then the dT tract. The paper does not name the
sequencing primers; standard Nextera primers are assumed.

**Model correction:** earlier versions drew generic TruSeq (RNA) and Nextera dual-index
(ATAC) libraries with the barcodes inside the insert and omitted the cDNA PCR and
Read1-only cDNA tagmentation. The off-target Tn5 products are now unresolved fates (not
amplified by Ad1/P7) rather than physical discards.

## 4. Evidence boundary

🟢 The defining paper and protocol specify the constant arms, linkers, modifications and 96 barcode values per round. The public page draws the exact constant geometry and modifications. 🟡 It uses structured 8-nt barcode placeholders rather than embedding all 288 well-specific values in every reaction state.

## Source-version discrepancy

🟢 The June 2020 [author preprint](https://doi.org/10.1101/2020.06.17.156943), Methods section (1), PDF text lines 718–725, describes an incoming **3′ overhang** and names the **round 2 linker** for the round 1 plate. The published [Cell paper](https://doi.org/10.1016/j.cell.2020.09.056), STAR Methods “Preparing oligonucleotides for ligations” (JATS line 213), instead describes a **5′ overhang** and the **round 1 linker** for that plate. 🟢 Published Supplementary Table S1, “bulk oligos” rows 3–5, prints separate round 1, 2 and 3 linkers. The model uses those published oligos and the final paper's round assignment. 🟡 Whether the preprint meant a different physical intermediate or contained wording errors cannot be determined from that text alone.
