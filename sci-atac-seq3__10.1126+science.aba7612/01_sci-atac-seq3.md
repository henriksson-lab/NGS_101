# sci-ATAC-seq3 — research note

Domcke S, Hill AJ, Daza RM *et al.* A human cell atlas of fetal chromatin accessibility.
*Science* 370:eaba7612 (2020). doi:10.1126/science.aba7612 (PMC7785298, author manuscript).

The earlier combined evidence trail is kept in the
[sci-ATAC-seq research note](../sci-atac-seq__10.1126+science.aab1601/01_sci-atac-seq.md).
This note records the sci-ATAC-seq3-specific reading behind the schematic.

## Sources read (bounded pass, 2026-10-17)

All three originals are in the shared family bundle
`$CHEM_DATA/sources/sci-atac-seq-family__10.1126+science.aab1601/`.

- Published article, PMC JATS `science.aba7612_PMC7785298.xml.txt`: lines 89 (design),
  197–209 (Online Methods: library preparation and sequencing), 293 (data and materials
  availability), 1037 (Fig. 1A legend). 🟢
- Table S7 (author workbook, mirrored original) `scg_aba7612_domcke_table-s7.xlsx.txt`:
  legend line 2, sheet `3lv2_splintandseqprimers` lines 3–10, first rows of the N5, N7,
  P5 and P7 sheets. 🟢
- Cited protocols.io protocol, doi:10.17504/protocols.io.be8mjhu6 (PDF export dated
  Nov 02 2020, "Last Modified: November 02, 2020", integer ID 35821; fetched as
  `protocolsio_sci-atac-seq3-ewov18xn7gr2.pdf.txt`): steps 5–13, lines 610–1160. 🟢
- Not read: the gated PMC supplementary PDF (NIHMS1652075 manuscript_supplementary.pdf;
  manual URL in the bundle MANIFEST).

## Oligos (Table S7, verbatim)

- `3LV2_N5_Splint` 5′-GCCGACGACTGATTA/3ddC/-3′ 🟢
- `3LV2_N7_Splint` 5′-CACGAGACGACAAGT/3ddC/-3′ 🟢
- N5 ligation oligo, e.g. `3LV2_N5_Oligo_A01_Plate1` 5′-CACCGCACGAGAGGTaagtagactaGTAATCAG-3′
  (384 members, 10-nt lowercase barcode) 🟢
- N7 ligation oligo, e.g. `3LV2_N7_Oligo_A01_Plate1` 5′-CAGCACGGCGAGACTggacgcctaaGACTTGTC-3′
  (384 members) 🟢
- P5 PCR, e.g. `3LV2_P5_PCR_A01` 5′-AATGATACGGCGACCACCGAGATCTACACGACCTCCTTGCACCGCACGAGAGGT-3′
  (96 members, 10-nt index) 🟢
- P7 PCR, e.g. `3LV2_P7_PCR_A01` 5′-CAAGCAGAAGACGGCATACGAGATcaaggcattcCAGCACGGCGAGACT-3′
  (96 members) 🟢
- Sequencing primers: `3lv2_seq_R1` TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,
  `3lv2_seq_IDX1` CTCCGAGCCCACGAGACGACAAGTC,
  `3lv2_seq_IDX2` ACACATCTGACGCTGCCGACGACTGATTAC,
  `3lv2_seq_R2` GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG. 🟢 The protocols.io step 13 list
  prints the same four sequences.

No 5′ phosphate is printed on the N5/N7 ligation oligos; ligation joins their 3′ end to
the Tn5-transferred 5′ end after T4 PNK phosphorylation. 🟡

## How they interlock (computed in `tools/sci3_graph.py`)

- Each splint is the reverse complement of the 8-nt ligation-oligo tail plus the first
  8 nt of the Nextera s5 or s7 end; the selftest derives both splints and compares them
  with Table S7. The 3′-ddC blocks splint extension. 🟡
- The P5/P7 PCR primers anneal through the 15-nt N5/N7 heads; only N5…N7 fragments
  carry both landing sites. 🟡
- `3lv2_seq_IDX1` and `3lv2_seq_IDX2` end on the N7 and N5 tails, so each index read
  starts on the 10-nt ligation barcode; the selftest reconstructs both from the model. 🟡

## Steps (article lines 209; protocols.io steps 5–13)

1. Tagment fixed nuclei with Nextera v2 enzyme (2.5 µl/well, 55 °C 30 min); stop with
   EDTA/spermidine, 37 °C 15 min; pool and wash. 🟢
2. T4 PNK with 1 mM rATP, 37 °C 30 min, nuclei aliquoted into four 96-well plates. 🟢
3. N5 ligation: T7 ligase buffer, N5 splint and T7 DNA ligase added to the PNK reaction,
   then one of 384 N5 oligos per well; 25 °C 1 h; EDTA/spermidine stop; pool and wash. 🟢
4. N7 ligation: master mix with N7 splint and T7 DNA ligase, 384 N7 oligos, 25 °C 1 h;
   stop; pool. 🟢
5. Resuspend in EB, dilute to 1,000–3,000 nuclei per well (four 96-well plates); reverse
   cross-links with EB/Proteinase K/1% SDS, 65 °C 16 h. 🟢
6. Test PCR with SYBR to set the cycle number; then NEBNext High-Fidelity 2X master mix,
   BSA, indexed P5 (column) and P7 (row) primers: 72 °C 5 min, 98 °C 30 s, 10–20 × (98 °C
   10 s, 63 °C 30 s, 72 °C 1 min). 🟢 The initial 72 °C step is modelled as gap fill-in. 🟡
7. Pool, Zymo Clean & Concentrator-5, then 1X AMPure; TapeStation; 2 nM pool. 🟢

## Read layout

Article line 209 (NovaSeq 6000, custom recipe and Table S7 primers): Read 1 51, Read 2 51,
Index 1 10 + 15 dark + 10, Index 2 10 + 15 dark + 10 cycles. 🟢 The protocols.io NextSeq
QC recipe (`3LV2_sciATAC_high`, step 13) lists 50/50 genomic bases and the same
10 + 15 dark + 10 index structure; this is a separate QC instrument recipe, not a
conflict with the NovaSeq production runs. The page uses the article values.

## Open questions

- The gated published supplementary PDF has not been read; it may add reaction detail.
- Whether a separate Tn5-gap fill-in is intended at the PCR's initial 72 °C step is not
  stated; the model labels that role inferred.
- The 15 dark cycles cover the N5/N7 heads; which recipe step makes them dark is in the
  custom recipe, not in the inspected sources.
