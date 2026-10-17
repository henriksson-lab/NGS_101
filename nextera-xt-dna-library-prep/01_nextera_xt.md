# Nextera XT DNA Library Prep

## Sources read

🟢 Illumina, *Nextera XT DNA Library Prep Product Documentation*, document 15031942 v07, pp. 1–15. It defines one-step enzymatic fragmentation/partial-adapter transfer, neutralization, limited-cycle index PCR, cleanup and normalization.

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, supplies the canonical Nextera mosaic end, S5/S7 entry points and sequencing primers through `lib/nextera.py`. Its section "Sequences for Illumina Prep and Nextera Kits" (extracted text lines 266–283) prints the Nextera transposase adapters "used for Nextera tagmentation" and the Index 1/Index 2 PCR primers; the N7xx/S5xx tables (lines 13318–13345, 13505–13521) give 8-base Nextera XT Index Kit v2 indexes.

🟢 Illumina, *Nextera XT Library Prep: Tips and Troubleshooting*, publication 770-2015-015, pp. 1 and 5–6, confirms that transposomes carry partial adapters, limited-cycle PCR completes the sequencing arms, and bead normalization yields single-stranded library. The schematic endpoint is the double-stranded library after PCR cleanup, before optional bead normalization.

## Molecular workflow

1. 🟢 Tagment: 10 µl TD + 1 ng DNA (5 µl) + 5 µl ATM (Amplicon Tagment Mix, the transposome); TAG program 55 °C for 5 min, hold 10 °C (guide pp. 11–12). The guide's "15 minutes" headline (p. 1) is the overall tagmentation stage, not the incubation.
2. 🟢 Add 5 µl NT (Neutralize Tagment Buffer) immediately at 10 °C "because the transposome is still active"; 5 min at room temperature.
3. 🟢 Amplify: 5 µl i7 + 5 µl i5 index adapter (or 10 µl from a plate) + 15 µl NPM; NXT PCR 72 °C 3 min, 95 °C 30 s, 12 × (95 °C 10 s, 55 °C 30 s, 72 °C 30 s), 72 °C 5 min (pp. 12–13). 🟡 The 72 °C opening step is consistent with filling the transposition gaps, but the guide does not say so; the model draws the tagmented intermediate as a closed duplex.
4. 🟢 Clean up: 50 µl PCR + 30 µl IPB (3:2) single-sided bead purification, two 80% EtOH washes, elute 52.5 µl RSB, keep 50 µl (pp. 13–14; 300–500 bp amplicons use 1.8×, > 500 bp 0.6×). The technical note (p. 5) states fragments < 200 bp, residual adapter and primer dimers are discarded. Libraries are then bead-normalized or quantified and pooled.

## Scope and endpoint

🟢 The guide recommends Nextera XT for genomes below 5 Mb, plasmids, double-stranded cDNA and PCR amplicons longer than 300 bp; Illumina DNA Prep is the alternative for large or complex genomes.

🟡 The product guide and technical note describe partial-adapter transfer but do not list the ATM loaded strands. The model takes the transferred Read 1/Read 2 transposase adapters (S5+ME, S7+ME) from the generic Nextera section of the adapter guide; their applicability to the XT ATM lot, and the non-transferred (ME-complement) strand, are not stated in the inspected documents and remain open. The final library has the public Nextera S5/ME–insert–ME/S7 architecture. The model displays the heterologous S5/S7 product because it can acquire both P5 and P7 and amplify.

🟢 The guide offers Nextera XT v2 and other compatible dual-index kits. The index-reference document prints the N7xx and S5xx families, while the pooling guide gives valid combinations. A particular pair is selected for each sample; it is not a fixed sequence of the prep chemistry. The model's PCR uses the printed Index 1 (P7–[i7]–S7) and Index 2 (P5–[i5]–S5) primers with 8-N index windows, annealing at S7/S5.
