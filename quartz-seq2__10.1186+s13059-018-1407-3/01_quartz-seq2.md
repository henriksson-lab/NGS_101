# Quartz-Seq2 — research-note pointer

Quartz-Seq2 is documented from Sasagawa *et al.*, 2018,
doi:[10.1186/s13059-018-1407-3](https://doi.org/10.1186/s13059-018-1407-3).

The detailed source record was originally written jointly for the two related protocols
and is retained as provenance in the
[combined Quartz-Seq / Quartz-Seq2 research note](../quartz-seq__10.1186+gb-2013-14-4-r31/01_quartz-seq.md#b-quartz-seq2-2018).
The public schematic itself uses only the Quartz-Seq2 paper and its supplementary oligo
table for its molecular construct.

## Published Table S4 discrepancy

🟢 Additional file 5, Table S4 (MOESM5 XLSX text lines 389–393), heads the v3.2
sheet as “SeqLv15mer” but its next line says “14 mer.” The first printed v3.2
barcode, eMDRT0001 `ATCAATCCTATCTGC`, is 15 bases; the defining article's
Methods (“Design of cell barcodes”) likewise specifies 15 nt for v3.2.
🟡 The schematic models a 15-base v3.2 barcode because the actual oligos and
Methods agree. The 14-mer line appears to be a table heading error; a formal
correction notice was not found in the inspected sources.

## Bounded curation round 2 (2026-10-16)

Inspected: Original PMC5845169 JATS60,72,236–240,256–268, stored in Quartz-seq source bundle.

Added 11-cycle suppression PCR. Source confirms selected v3.2 read mask23/6/63; v3.1 is a different mask. Public notice qualifies RNase H activity order in the tailing mixture.

Deferred: Original primer tables and supplements; combined RNase H/TdT graph activity; dedicated source bundle absent, original read from existing external Quartz-seq collection. This capped pass does not establish a complete source review or publication-wide absence.
