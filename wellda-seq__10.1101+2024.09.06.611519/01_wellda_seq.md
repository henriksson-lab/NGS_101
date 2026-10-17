# wellDA-seq

Wang K, Yan Y, Elgamal H *et al.*, Navin N. "Single cell genome and epigenome co-profiling reveals hardwiring and plasticity in breast cancer." *bioRxiv* 2024, [doi:10.1101/2024.09.06.611519](https://doi.org/10.1101/2024.09.06.611519), v1 posted 2024-09-10 (PMC mirror PMC11418942). On 2026-10-17 the source fetcher linked no published version, so this preprint and its supplement (media-1) are the authority.

## Sources read

- 🟢 Preprint v1 full text (`2024.09.06.611519_v1.full.pdf.txt`): Results "Overview of the wellDA-seq method" (lines 113–131); Methods "Tn5 transposome assembly" and "wellDA-seq procedure" Steps 1–3 (lines 816–911); "Preprocessing" (lines 926–929).
- 🟢 Supplement media-1: Extended Data Fig. 1 (text lines 20–82, with the rendered page 2 checked), Supplementary Table 5 (ATAC barcode primers, from line 3384) and Table 6 (DNA barcode primers, from line 3868).
- Bundle oligo scraper run (49 TSV rows); every oligo below was checked against the passages above.

## Chemistry (step by step)

1. 🟢 **Tn5-1 (ATAC).** DEpi_ME_t1 `GCCTCCCTCGCGCCATAGATGTGTATAAGAGACAG` and DEpi_ME_t2 `CTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG` are each annealed to ME reverse `5′-/phos/CTGTCTCTTATACACATCT-3′` (95 °C 3 min, 0.1 °C/s to 25 °C), diluted to 25 µM and loaded on Diagenode Tagmentase (C0107001020) with glycerol.
2. 🟢 Cells are permeabilized (NP40/Tween-20/0.01% digitonin, 4–5 min) and about 30–35k nuclei are bulk-tagmented with Tn5-1 in 2× TD buffer, 37 °C 30 min.
3. 🟢 Tagmented nuclei (DAPI) are dispensed into ICELL8 350v chips (5,184 nanowells). Imaging and CellSelect keep wells holding one nucleus.
4. 🟢 In each selected well, protease lysis removes chromatin (54.5 °C 30 min, then 75.3 °C 15 min to inactivate).
5. 🟢 **Tn5-2 (DNA).** A second tagmentation uses Illumina TDE1 (20034198; Nextera S5/S7 adapters), 54.5 °C 10 min, tagging the previously closed genome.
6. 🟢 A neutralization mix (EDTA, dNTP) carries one well-specific pair of S5XX primers (ATAC_CB1, DNA_CB1). A second mix carries the N7XX pair (ATAC_CB2, DNA_IN CB2). KAPA HiFi master mix is then added. On-chip program: 72.1 °C 8 min (gap fill), then 12 PCR cycles. Wells are collected by centrifugation and cleaned with 1.5× AMPure XP.
7. 🟢 30 ng aliquots are enriched separately. ATAC uses Bioo-PCR-F (P5) + Bioo-PCR-R (P7) for 12 cycles. DNA uses Bioo-PCR-F + wellDA_DNA_OUT_N7XX (P7–DNA_CB3–`CTGAGTCGGAGACACGCA`) for 7 cycles. Each is followed by 1.2× AMPure.
8. 🟢 Sequencing: NextSeq 2000. The ATAC library uses four spiked-in custom primers: Read 1 DEpi_ME_t1, Read 2 DEpi_ME_t2, Index 1 Idx7 `CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAG`, Index 2 Idx5 `CTGTCTCTTATACACATCTATGGCGCGAGGGAGGC`. 🟡 Idx7 and Idx5 are the reverse complements of the two DEpi adapters (enforced in the module). bcl2fastq splits ATAC and DNA by index, allowing no mismatches.

## Final libraries

- 🟢 ATAC (Extended Data Fig. 1a): P5–TGA–ATAC_CB1–DEpi handle 1–ME–insert–ME–DEpi handle 2–ATAC_CB2–P7. With the custom primers, Index 1 reads ATAC_CB2 and Index 2 reads ATAC_CB1 (derived on the page).
- 🟢 DNA (Extended Data Fig. 1b): P5–DNA_CB1–Nextera S5–ME–insert–ME–S7–DNA_CB2–DNA_IN link–DNA_CB3–P7.
- 🟢 Supplementary Tables 5–6 list 72 ATAC-CB1, 72 ATAC-CB2, 72 DNA-CB1, 72 DNA-CB2 and 11 DNA enrichment (CB3) primers. The selftest counts these families against the module's flanking sequences.

## Published discrepancies (shown on the page)

- Methods Step 2 (4)–(5) cites the barcode primers as "Supplementary Table 4-5". In the supplement, Table 4 holds tumour QC metrics and the barcodes are in Tables 5 (ATAC) and 6 (DNA); Extended Data Fig. 1 also cites S5/S6. The model cites Tables 5 and 6.
- Extended Data Fig. 1a draws the ATAC products as `…ATAC_CB2-ATCACGATCTCGTATGCCGTCTTCTGCTTG-3′`, with six extra bases (ATCACG) before the P7 complement. The Methods N7XX template, the figure's own primer list and all 72 Table 5 ATAC-CB2 primers lack them. The model follows the primers. The source does not say whether the figure's bases are an artefact.

## Corrections in round 2026-10-17-breadth-first-21

The previous role-level model drew one indexed Tn5, a "separation" of tagged products from residual genome, and a strand-displacing whole-genome amplification. The source describes none of these steps. The model now follows the source: two Tn5 transposomes (bulk Tn5-1, then in-well Tn5-2 after protease lysis), a single four-primer on-chip barcoding PCR, a pooled cleanup and two modality-specific enrichment PCRs. It also has exact ATAC and DNA final libraries and the ATAC custom sequencing primers.

## Open questions

- 🔴 Read lengths are not stated. Neither are the DNA library's sequencing primers or its index-read cycles. Reading both DNA_CB2 and DNA_CB3 from one index read would need more than 8 cycles; this is not stated.
- 🟡 The role of DNA_CB3 relative to the CB1/CB2 well pair is not stated. The supplement only calls it a DNA cell barcode, and 11 primers exist.
- 🟡 Cross-modality products are possible in the shared on-chip PCR, for example ATAC_S5 on one end and DNA_IN on the other. Such a hybrid would carry the DNA enrichment primer site. The source does not discuss them, so they are not drawn.
- 🟡 Some fragments carry the same DEpi adapter at both ends (t1/t1, t2/t2). In the on-chip PCR they would take the same barcode primer at both ends and lack the P5/P7 pair needed for sequencing. The source does not discuss their fate, so they are not drawn.
