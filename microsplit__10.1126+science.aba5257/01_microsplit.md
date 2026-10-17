# microSPLiT — source note

The complete source audit and oligo transcription are retained in the companion
[SPLiT-seq / microSPLiT research note](../split-seq__10.1126+science.aam8999/01_split-seq.md).

**Defining source:** Kuchina A *et al.* “Microbial single-cell RNA sequencing by
split-pool barcoding.” *Science* 371, eaba5257 (2021).
doi:[10.1126/science.aba5257](https://doi.org/10.1126/science.aba5257).

The schematic uses a mirrored copy of the authors' 2021 *Science* Supplementary
Table S3 (`scg_aba5257_table_s3.xlsx`, sheet `General Oligonucleotides`,
lines 1–19). The workbook identifies itself as Table S3; its embedded metadata
names Anna Kuchina as the last editor. The original publisher attachment is still
gated, so byte identity of this mirror has not been verified. The table supplies
the split-pool oligos, BC_0243/BC_0244 adapter duplex, and BC_0027 plus indexed
BC_0076–BC_0083 terminal PCR pairs. The page instantiates BC_0076 (TSBC07)
as a representative final-library index. BC_0062 and BC_0108 cDNA-amplification
primers are printed in the later 2024 detailed protocol, not in the mirrored 2021
Table S3.

**Adapter source-version difference.** The mirrored 2021 Table S3, General
Oligonucleotides row BC_0243 (line 15), prints
`ACACTCTTTCCCTACACGACGCTCTTCCGATC*T` with the adjacent note
`*phosphorothioated base` 🟢. Its BC_0244 opposite strand (line 16) is
`GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` 🟢; the first 32 bases pair, leaving
BC_0243's terminal T as the 3′ overhang 🟡. The 2019 bioRxiv v2 Auxiliary
Table (`869248_media-2.xlsx.txt`, lines 14–15) and 2024 *Nature Protocols*
Table 3 (`...MOESM2_ESM.xlsx.txt`, lines 121–122) print the same BC_0243 bases
without the asterisk 🟢. The molecular model uses the 2021 Table S3 terminal
C*T phosphorothioate annotation. The reason for the notation difference is
unknown; the official 2021 attachment remains inaccessible for direct comparison.

**Later primary-source cross-check (2024).** The authors' detailed Nature Protocols
procedure, doi:[10.1038/s41596-024-01007-w](https://doi.org/10.1038/s41596-024-01007-w),
explicitly describes fragmentation followed by ligation of the BC_0243/BC_0244
adapter duplex and PCR with BC_0027 plus an indexed reverse primer (Figure 2 and
Part 2, steps 134–157). Its Supplementary Table 3 (`...MOESM2_ESM.xlsx`, sheet
`Table 3 General Oligonucleotide`) prints the same representative BC_0027,
BC_0062, BC_0108, BC_0076, BC_0127 (`rGrG+G`), BC_0243 and BC_0244 sequences
used in the model 🟢. It additionally gives **24** indexed final-PCR primers:
BC_0076–BC_0083 and BC_0246–BC_0261 🟢. Thus the 24 sublibraries mentioned
in the earlier source can be assigned distinct published primer sequences;
the model's BC_0076 remains a representative rather than a complete index list.
No conflicting shorter BC_0076 exists in the 2024 workbook: its `Table 1 Barcodes`
cost sheet lists BC_0076–BC_0083 only with concentration, `Custom IDT`, price (19.2 USD)
and cost per run, not sequences, and every BC_0076 sequence in the bundle (2021 Table S3,
2024 Table 3, 2018 SPLiT-seq Table S12) is the same 64-nt primer 🟢 (checked 2026-10-17).
The 2024 procedure confirms this later implementation. The mirrored 2021 table
allows a targeted sequence comparison, while direct verification against the
publisher-hosted 2021 attachment remains blocked 🔴.

**Bounded pass 2026-10-17 (round breadth-first-12): published Methods order.**
The published *Science* author manuscript (PMC8269303 JATS, `science.aba5257_PMC8269303.xml.txt`)
was re-read at Methods lines 149 (formaldehyde fixation, Tween-20 and lysozyme
permeabilization), 153 (in-situ polyadenylation with NEB *E. coli* Poly(A) Polymerase I,
ATP, 37 °C 30 min), 165 (round-1 in-cell RT), 169 (T4 DNA ligase rounds 2–3) and 175
(sublibrary aliquots, −80 °C freeze, SDS lysis buffer + proteinase K at 55 °C 2 h) 🟢.
Line 165 states that every RT well holds both random-hexamer and anchored poly(dT)15
barcoded primers at 1:2 (2.5 µM + 5 µM), with Maxima H Minus RT, 10 min 23 °C then
50 min 50 °C 🟢. The graph previously started RT directly on the permeabilized cells
with only the dT primer; it now declares the PAP tailing step, both RT primer families
and the lysis/proteinase K step before streptavidin capture (model correction).

The mirrored 2021 Table S3, sheet `Round 1 RT Barcodes`, rows E1–H12 (lines 70–117),
prints the hexamer primers, e.g. Round1_49 `\5Phos\ACTGTGGCTGCTTTGTNNNNNN` 🟢: the
7-nt handle, an 8-nt barcode, one printed T, then N6. The graph uses Round1_49 as the
representative hexamer primer. The 2024 *Nature Protocols* Table 3 instead gives a
separate `R1 Random Hex Barcodes` sheet whose Round1_01 is `/5Phos/ACTGTGGACTCGTAANNNNNN`
(line 226, no T, same barcode as dT Round1_01) and reuses Round1_49's barcode for a dT
primer (line 176). This is a later-protocol design change, not a published discrepancy;
the page follows the 2021 implementation. Only the dT-primed first strand is drawn (the
2024 Figure 2 legend, line 1661, makes the same simplification).

**Read lengths.** The model's Read 1 74 nt, Read 2 86 nt and 6-nt index come from the
bioRxiv v2 supplement (`869248_media-1.pdf.txt` lines 135–136) and agree with the 2024
procedure step 172 (`s41596...xml.txt` line 1206) 🟢. The published 2021 Supplementary
DOCX that would state the published run is still gated, so the 2021 acquired run
lengths remain unverified against the publication 🔴.
