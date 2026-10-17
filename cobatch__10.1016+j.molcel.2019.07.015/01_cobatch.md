# CoBATCH — research notes

Defining source: Wang et al., 2019, [doi:10.1016/j.molcel.2019.07.015](https://doi.org/10.1016/j.molcel.2019.07.015).

## Sources read

- 🟢 Superseded author preprint 10.1101/590661 Methods (historical provenance; the published STAR Methods now govern): antibody targeting, PAT binding, barcoded tagmentation, pooling, redistribution and indexed PCR.
- 🟢 Preprint Supplementary Table 2 workbook (`590661_media-3.xlsx`): Mosaic-A/B, phosphorylated common ME oligo, PAT-T5 and PAT-T7 barcode series, Nextera i5/i7 PCR primers.
- 🟢 Recovered original *Molecular Cell* Table S1 (`mmc2.xlsx`), sheets Mosaic-A&B / PA-Tn5_barcodes / Nextera_index: confirms the adapter and PCR families previously read in preprint Table 2.
- 🟢 Original supplemental PDF (`mmc1.pdf`) contains Figures S1–S4 and captions. Original Table S2 (`mmc3.xlsx`) records antibodies, library-preparation labels and sequencing statistics; Table S3 (`mmc4.xlsx`) compares CoBATCH/scChIC-seq performance, and Table S4 (`mmc5.xlsx`) gives tissue-digestion recipes. These are not additional sequencing-primer tables.
- 🟢 Author-deposited GEO records (round 2026-10-17-deferred-followup-02): EC_CoBATCH sub-series GSE129327 (`geo-GSE129327-self-quick.txt`, series summary/design/relations, SubSeries of GSE129335, PubMed 31471188) and fixed EC sample GSM3705255 (`geo-GSM3705255-self-quick.txt`, lines 14–27, 37).
- 🟢 Published *Molecular Cell* article with supplemental information (`1-s2.0-S1097276519305453-mmc7.pdf`, hand-fetched 2026-10-17 from the Elsevier CDN; printed title, authors, Mol Cell 76:206–216 and DOI checked), round deferred-followup-03: STAR Methods lines 1418–1528 of the `.txt` (Quality control for the assembled PAT transposome, Low-input in situ ChIP, Nextera library preparation and sequencing strategy, Fixed and Native single-cell CoBATCH), Key Resources "Oligonucleotides" (line 1346: refers to Table S1), Additional Resources (Methods S1).
- 🟢 Published Methods S1 bench protocol (`1-s2.0-S1097276519305453-mmc6.pdf`, steps 1–27 and troubleshooting).
- 🔴 Predecessor sources (Chen et al. 2016, Picelli et al. 2014 Tn5 assembly) were not re-read; no missing-information conclusion is accepted.

## Molecular path

🟢 Antibody-bound Protein-A–Tn5 deposits barcoded adapters at target chromatin in a first plate. The cells are pooled and redistributed; i5/i7 PCR indexes encode the second plate. A cell is identified by the combined T5, T7, i5 and i7 values.

🟢 The first-plate oligos are complete transposon top strands, not bare barcode windows. T5 has S5 + `TCCACGC` + 8-nt T5 barcode + `GCGATCGAGGACGGC` + ME. T7 has S7 + `CTGTCCCTGTCC` + 8-nt T7 barcode + `CACCGTCTCCGCCTC` + ME. Both anneal only over their terminal ME to the 5′-phosphorylated common bottom `CTGTCTCTTATACACATCT` before pA–Tn5 loading.

🟢 The source workbook contains 24 T5 and 25 T7 barcode oligos. The schematic draws the exact T5-1 and T7-2 representatives and attaches whitelist feature metadata; it does not render repeated placeholder bases as if they were a published barcode.

🟢 After targeted transfer and gap fill, the first-plate codes sit between their outer and inner connectors at both fragment ends. Indexed PCR adds P5–i5–S5 and S7–i7–P7, making the four-part cellular identity explicit in the completed duplex.

## Reaction conditions (published STAR Methods; round deferred-followup-03)

🟢 The published STAR Methods (mmc7 pp. e4–e5, "Fixed single-cell CoBATCH" and
"Nextera library preparation and sequencing strategy") and Methods S1 (mmc6
steps 11–27) give the same fixed single-cell conditions as the preprint read
earlier: 3 µg/ml PAT-T5 and PAT-T7 per well with 2,000 cells, 10 µl TAPS-NaOH
pH 8.3 / 5 mM MgCl2 / 10% DMF at 37 °C for 1 h, 40 mM EDTA stop, pooling,
DAPI, 30 µm strainer, 20–25 cells into 4 µl SDS/proteinase K lysis buffer,
65 °C 5–6 h, 85 °C 15 min, 1 µl 1.8% Triton X-100 at 37 °C 1 h, then in-well
2X KAPA PCR: 72 °C 5 min, 98 °C 45 s, 10–18 cycles of 98 °C 15 s / 63 °C 30 s
/ 72 °C 1 min, 72 °C 5 min. New from Methods S1 step 24: all wells are combined
after PCR and purified on a TIANGEN column (50 µl elution) before the 1X AMPure
cleanup and 0.5X/0.5X double-sided selection (steps 25–26); the graph now draws
this pooled column step. STAR Methods (p. e4) transposome assembly: each T5 or
T7 oligo mixed with MErev, 95 °C 5 min, 0.1 °C/s to 75 °C, cooled to RT; 37.5 µM
adaptor + 37.5 µM PAT at 25 °C for 60 min. Native single-cell CoBATCH differs
by no fixation, tagmentation at 25 °C 1 h and lysis at 65 °C 10 min (not drawn).

Earlier preprint-scoped record (round 2026-10-17-breadth-first-08): the
author preprint (10.1101/590661v1, `590661_v1.full.pdf.txt` lines 470–495 and
423–452) describes, for fixed single-cell CoBATCH, PAT-T5/PAT-T7 binding in
96-well plates (2000 cells per well), activation in TAPS/MgCl2/DMF at 37 °C
for 1 h, EDTA stop, pooling, DAPI staining and straining, sorting of 20–25
cells into each second-plate well containing SDS/proteinase K lysis buffer
(65 °C, then 85 °C inactivation, Triton X-100 quench), then PCR directly in
the well with 2X KAPA master mix: 72 °C 5 min (gap fill), 98 °C 45 s,
10–18 cycles, final 72 °C 5 min; 1X AMPure XP cleanup and 0.5X/0.5X
double-sided selection for 200–1000 bp. The graph now places gap fill inside
the PCR well after redistribution (it was previously drawn before pooling)
and represents lysis and bead selection with discarded outputs.

## Published Table S1 annotation offset

🟢 In the original publisher workbook (`mmc2.xlsx`, sheet Mosaic-A&B) the
sequences sit in cells B3–B5 (Mosaic-A, Mosaic-B, common annealing primer
`5Phos/CTGTCTCTTATACACATCT`), but the equivalence notes sit in C4–C6:
"equally to Tn5ME-A_FC-121-1030" beside Mosaic-B, "equally to
Tn5ME-B_FC-121-1031" beside the common primer, and "equally to
TnExtMErev_primer" on an otherwise empty row 6. By sequence Mosaic-A is the
S5/Tn5ME-A arm, Mosaic-B the S7/Tn5ME-B arm, and the common primer the ME
reverse strand (the T5/T7 rows say they anneal with TnExtMErev_primer). The
model uses the printed sequences; the page shows this discrepancy.

## Sequencing evidence correction

🟢 Original published Table S1 prints the T5/T7 transfer oligos and i5/i7 PCR
primers. Its original sheet rows agree with the representative complete oligos
assembled on this page. The author preprint Methods, lines423–452, specify
paired-end 150-base reads for the Nextera library strategy.

🟡 The former page declared inner-connector-plus-ME sequences as physically used
custom genomic read primers, and showed 50-base reads. Those were unsupported
repository assignments. The displayed adapter subsequences alone cannot
establish the loaded sequencing primers. They have been removed; read length
is now taken from the published Methods S1 step 27 (paired-end 150 bp,
HiSeq X Ten or NovaSeq 6000). The published STAR Methods and Methods S1 print
only the i5/i7 PCR primers; they do not print loaded sequencing primers or
index-read lengths, so exact sequencing-entry/index-read geometry stays
unresolved (bounded-pass observation, not an accepted absence claim).
This correction is not attributed to a publication error.

## GEO sample record (deferred follow-up 02, 2026-10-17)

🟢 GSM3705255 (H3K27ac, left forebrain endothelial cells; "Fixed by formaldehyde")
states: extraction "Following CoBATCH protocol"; "After PCR enrichment using 2X KAPA
master mix, products were purified and selected by AMPure XP beads for 200-1000 bp"
(consistent with the preprint-scoped conditions above); instrument Illumina NovaSeq
6000; and "For spliting single cell, we cut the first 21 bp of the read1 and the first
27 bp of the read2" (lines 18–21, 37).

🟡 Derived with the model: if Read 1 began immediately 3′ of S5, its first 21 cycles
would cover the {{= len(cobatch.T5_OUTER) }}-nt T5 outer connector, the 8-nt T5 barcode
and {{= cobatch.split_overrun("Read 1") }} nt of the inner connector; if Read 2 began
3′ of S7, its first 27 cycles would cover the {{= len(cobatch.T7_OUTER) }}-nt T7 outer
connector, the 8-nt T7 barcode and {{= cobatch.split_overrun("Read 2") }} nt of the
inner connector. This is consistent with S5/S7-ending read primers but does not
establish the loaded primers, so none are declared. The record is one tissue (EC)
sample; whether other sub-series used the same split is not checked.

Retrieval history (deferred follow-up 02): Unpaywall/OpenAlex list it as
bronze open access at `http://www.cell.com/article/S1097276519305453/pdf` (HTTP 403
Cloudflare challenge, also for the cell.com full-text and ScienceDirect article/pdfft
pages); Elsevier article API without a key returns metadata only; Europe PMC lists it
as subscription-only (no PMCID). The bioRxiv 590661 version is the already-read
preprint. The maintainer later hand-fetched the published article PDF (mmc7)
and Methods S1 (mmc6) from the Elsevier CDN; both were read in round
deferred-followup-03.


## Critical re-review, 2026-10-17

Non-library tagmentation end-class candidates remain unresolved. PCR selectivity is not a represented physical separation.

Inspected scope: Local tagmentation-to-PCR declarations and terminal-fate audit; no new reaction-source assertion. This bounded correction does not complete the linked-source or two-reader absence audit.
