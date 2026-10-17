# Nextera Mate Pair

## 1. Protocol boundary

This page follows Illumina's Nextera Mate Pair Library Prep kit, document 15035209 v02.
The defining topology is a long tagmented fragment that is circularized, re-sheared and selected
for the marked old circle junction before a conventional single-index TruSeq library is made.

## 2. Sources read

🟢 Illumina, *Nextera Mate Pair Library Prep Reference Guide*, document 15035209 v02 (2016).

🟢 Illumina, *Data Processing of Nextera Mate Pair Reads on Illumina Sequencing Platforms*.
This technical note prints the junction and external adapter sequences and explains RF mapping.

## 3. Molecular state changes

Locations are in the reference guide `.txt` twin (v02, January 2016) unless the technote is named.

🟢 Mate Pair Tagment Enzyme fragments and tags long genomic DNA (55 °C, 30 min; Zymo column
purification). The technote (Figure 1) states that it attaches a **biotinylated junction adapter**
to both ends; the guide itself does not print the word biotin.

🟢 Strand Displacement (20 °C, 30 min) repairs "a short single-stranded gap" left by tagmentation.
A ~0.4x AMPure XP step then "remove[s] short fragments (< 1500 bp)" (gel-free: 50 µl reaction +
50 µl water + 40 µl beads). Earlier model text called the gap "nine-base" and discarded short
products at tagmentation; both were corrected (the gap length is not stated in the inspected guide).

🟢 Gel-free: 1 µg input, 2–15 kb, median 2.5–4 kb. Gel-plus: 4 µg, 100 µl beads on 200 µl, then
Pippin Prep or agarose size selection and Zymo recovery. Only the gel-free path is drawn.

🟢 Circularization is "a blunt ended intramolecular ligation", up to 600 ng in 300 µl, 30 °C
overnight; more DNA raises yield but also chimeric pairs. Exonuclease (9 µl) is added directly
to the circularization reaction: 37 °C 30 min, 70 °C 30 min, then Stop Ligation Buffer.

🟢 Covaris S2/S220 shearing (40 s, 6 °C) gives ~300–1,000 bp fragments with 3′ or 5′ overhangs.
Dynabeads M-280 Streptavidin (20 °C, 15 min, four Bead Wash Buffer and two RSB washes) retain the
adapter-containing fragments; everything after this is on bead until PCR.

🟢 The common duplicated junction printed by Illumina (technote Table 1) is
`CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG`, derived in the model as `ME_RC + ME`.
Single junctions (`CTGTCTCTTATACACATCT` / `AGATGTGTATAAGAGACAG`) occur less frequently. The
technote also states that incomplete enrichment leaves non-junction FR pairs as a < 1,000 bp
insert peak, and that shearing inside the junction adapter yields a fragment with only one
original end. These are shown as a page highlight, not as complete removal.

🟢 On bead: End Repair Mix (30 °C, 30 min), A-Tailing Mix (37 °C, 30 min; adapters carry a
complementary 3′ T), Ligation Mix with one of twelve TruSeq DNA LT Set A adapters (30 °C, 10 min,
Ligation Stop Buffer, washes). Each adapter has a 6-base index (seventh cycle base in parentheses
is not part of the index), sequenced in a single-index workflow (guide Index Adapter Sequences).
PCR with Enhanced PCR Mix and PCR Primer Cocktail is run on the beads, 10 or 15 cycles
(98 °C 30 s; 98 °C 10 s, 60 °C 30 s, 72 °C 30 s; 72 °C 5 min). The 45 µl supernatant is cleaned
with AMPure XP to remove library fragments < 300 bp.

🔴 The PCR Primer Cocktail and full TruSeq LT adapter sequences are not printed in the inspected
guide; the final TruSeq sites come from the technote external adapter sequences and `lib/illumina`.

## 4. Read layout

Read 1 and Read 2 enter genomic sequence from opposite ends of the short captured fragment. When
the fragment spans the junction, the genomic pieces on either side were the distant ends of the
original long molecule and align in outward-facing RF orientation. A sufficiently long read may
cross the internal junction; that sequence is trimmed during analysis.

## 5. Open questions

- Biotin position within the junction adapter and the full transposome adapter structure.
- PCR Primer Cocktail sequences (historical TruSeq DNA LT guide not inspected).
- Acquired read lengths are user-defined; none is assigned.
