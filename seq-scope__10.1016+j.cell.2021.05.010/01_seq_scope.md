# Seq-Scope — research notes

Defining source: Cho *et al.*, [10.1016/j.cell.2021.05.010](https://doi.org/10.1016/j.cell.2021.05.010).

## Evidence used

- 🟢 Seq-Scope uses solid-phase amplification of randomly barcoded HDMI oligos on an Illumina flow cell.
- 🟢 A first sequencing/imaging run maps every HDMI cluster to a physical coordinate.
- 🟢 Tissue RNA is captured through oligo(dT), reverse-transcribed on the decoded surface, and recovered with spatial barcode and UMI for transcript sequencing.
- 🟢 The open defining paper prints the complete HDMI-DraI and HDMI32-DraI seed oligos, including P5/P7, TruSeq Read 1, constrained-random HDMI, HR1, oligo(dT), DraI site and DraI-adapter regions.
- 🟢 It prints the full `Read1-DraI` custom sequencing primer used for the coordinate-decoding run and the TruSeq Read2-conjugated 9-nt random primer used for second-strand synthesis; the latter random tract serves as the UMI.
- 🟢 Both amplification rounds are source-exact, including their 3′ phosphorothioate bonds and the 8-mer sample-index window in the final reverse primer. The completed library has the standard TruSeq sequencing-primer sites.
- 🟢 The paper also prints the experimental HDMI-TruEcoRI seed, Read1-EcoRI custom sequencing primer, and the 30-dA/8-nt-UMI oligo used to test a UMI-encoded array. These are retained as a distinct alternative rather than blended into the main DraI workflow.

## Bounded source pass 2026-10-17 (round breadth-first-18)

Read: published author manuscript PMC8238917 (JATS XML twin), STAR Methods lines 211–285
(seed library, 1st-Seq, HDMI-array processing, disassembly, tissue, RT, Exo I, tissue
digestion, second strand, library PCR), lines 289–323 (UMI-encoded EcoRI array), analysis
lines 341–343 and Discussion line 173. The published supplements (NIHMS1704034 supplements
1–12) remain behind the PMC download gate; the bioRxiv preprint is not used for active claims.

- 🟢 1st-Seq: HDMI-DraI loaded at 100 pM (HDMI32-DraI 60–80 pM) on MiSeq v3, manual mode,
  25-cycle (HDMI-DraI) or 37-cycle (HDMI32-DraI) single read with Read1-DraI; no denaturation,
  indexing or re-synthesis; flow cell retrieved. The FASTQ holds the HDMI followed by 5 bases of TR1.
- 🟢 Array processing: DraI (1 U, CutSmart, 37 °C overnight) "to completely cut out the P5
  sequence and expose oligo-dT"; Exo I 45 min (the Read1-DraI duplex protects the P7-bound
  strand); water, 0.1 N NaOH (3 × 5 min, removes Read1-DraI), Tris and water washes.
- 🟡 DraI cuts TTT|AAA bluntly, so the exposed capture end is the printed 27 T plus TTT (30 T).
- 🟢 Tissue: 10 µm section, 4 % formaldehyde, H&E imaging, collagenase I then pepsin/HCl release;
  RT with Maxima H− (Ficoll PM-400, actinomycin D, 42 °C overnight); Exo I; proteinase K digestion;
  NaOH washes remove mRNA.
- 🟢 Second strand: TruSeq Read2-conjugated N9 random primer with exo− Klenow, 37 °C 2 h; water
  washes; 2 × 30 µL 0.1 N NaOH elution, KOAc neutralization; AMPure XP 1.8×.
- 🟢 PCR1: KAPA HiFi HotStart, 13–15 cycles, AMPure 1.2×. PCR2: 8–9 cycles from 10 µL of 2 nM,
  gel excision 400–850 bp (Zymoclean), AMPure 0.6–0.7×. Paired-end 100–150 bp on Illumina and BGI.
- 🟢 Read layout: Read 1 positions 1–20 (HDMI-DraI) or 1–30 (HDMI32-DraI) are the HDMI, matched to
  the reverse-complemented 1st-Seq whitelist; the first 9 bases of Read 2 (random primer) are the UMI.
- 🟡 The final library is computed from the printed oligos:
  P5–TruSeq Read 1–HDMI–HR1–30 T–transcript–N9–TruSeq Read 2–i7–P7. It is single-indexed (the
  printed P5 forward primer has no i5 window). The previous model used a generic dual-index shell
  with the UMI next to the HDMI; both errors are corrected.
- **Published label vs bases (shown on the page):** the backbone (Methods line 217) labels the seed
  ends "P5 sequence" (5′) and "P7 sequence" (3′), but the printed bases (lines 221–225) are the
  Illumina P7 adapter and the reverse complement of P5. The model names segments by their bases.

Open questions: the published supplements are still unread; i7 index-read length and the exact
read lengths per dataset are not stated; HDMI32-DraI analysis uses 30 of its 32 HDMI bases (line 341),
reason not given; on-instrument bridge amplification and linearization are not drawn.

## Protocol boundary

The defining two-pass use of the Illumina flow cell is retained rather than reducing Seq-Scope to a generic spatial oligo(dT) array.

All paper-published oligos are now listed on the schematic page. The main final-library view follows the DraI implementation used for the liver and colon data; the EcoRI set is explicitly an experimental alternative.
