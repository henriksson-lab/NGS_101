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

## Native graph migration (2026-10-18, batch 12)

Read one original: published author-manuscript `j.cell.2021.05.010_PMC8238917.xml.txt`, STAR Methods lines 211–285 (seed/first sequencing/array processing/disassembly/tissue/RT/digestion/second-strand elution and both PCR rounds), and analysis lines 337–345 (coordinate join, reverse-complemented bottom-tile whitelist, N9 processing). Full bundle oligo scraper run; only the published original was used for the positive seed, Read1-DraI, N9 primer and phosphorothioate-protected PCR oligos, verified at lines 221, 233, 271, 279 and 281. Preprint scraper hits do not supply active claims.

The selected native route remains HDMI-DraI20 for liver/colon. Removed the legacy protocol-spec wrapper and shared source-selecting molecular recipe. Native ssDNA seed, on-surface capture DNA, independent SBS daughter, tissue RNA, extended original capture strand, N9-primed second DNA strand and PCR daughters are separate polymers/populations. First-run strand attachment and on-instrument linearization remain explicitly unresolved. The reference second-strand copy does not assume complete extension into the P7 seed arm near the surface; its further 3′ extent remains unassigned.

DraI/Exo I chemical products remain in the mixture until actual wash fractions remove them. Protein digestion likewise does not itself discard tissue; actual washes remove debris and RNA/hydrolysate populations. Disassembly explicitly retains bottom tiles and destroys top tiles. Array-anchored original first-strand DNA remains on the unused surface when NaOH elutes the independent second-strand copy. Cleanups expose physical fractions without proving complete class purity or single-copy recovery.

Corrected PCR1 input: the actual eluted single-stranded second-strand mixture supplies PCRGeometry, rather than a reconstructed surface duplex. PCR2 uses a physical 10 µl aliquot of 2 nM cleaned PCR1 material; the remaining preparation has unresolved fate. Exact primer sequences are retained, while native first-cycle/mixed-class PCR geometry is unresolved. Molecules directly incorporating original protected primers retain the phosphorothioate bonds; copying those bases does not copy the chemical bonds. The source-derived final read-site reference remains separate from incoming templates. Removed the arbitrary 28-base transcript placeholder length.

The first-run 25-cycle HDMI20 read and bottom-tile reverse-complement whitelist join stay distinct from the second-run HDMI20/N9/i7 acquisition. Source paired-end range is 100–150 bases; exact dataset cycles, i7 read length, pooled companion libraries/ratios and library recovery by molecular class remain unassigned. Copying N9 into Read1, polyA trimming, mapping and UMI collapse are computational operations, not physical separation.

Deferred: on-instrument bridge/linearization and complete anchored-end chemistry; full native SBS/cleavage geometry; source-specific PCR first-cycle products and PS bond inheritance by individual daughter; detailed capture/second-strand recovery; dataset acquisition and actual pooled-library identities/ratios; unread published supplements and complete linked-source review. HDMI32/EcoRI alternatives are not expanded. This bounded migration leaves existing absence dispositions unchanged.
