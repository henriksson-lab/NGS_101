# ICICLE-seq — source notes

Defining paper: Swanson et al., *eLife* 2021;10:e63632 (published version),
[doi:10.7554/eLife.63632](https://doi.org/10.7554/eLife.63632). The published bundle
(PMC8034981 JATS and eLife v2 supplements) is shared with TEA-seq; it was fetched to the
ICICLE-seq source slug on 2026-10-17 (supplement 5 text identical to the TEA-seq copy).

## Sources inspected (2026-10-17 bounded pass)

- PMC8034981 JATS `.xml.txt`: Results line 246 (ICICLE-seq description); Figure 3—figure
  supplement 1 legend (a–e), line 267; Methods "ICICLE-seq" lines 702–722
  (oligonucleotide sources, Tn5 complexing, antibody staining, library preparation,
  sequencing) and data preprocessing lines 724–728.
- Supplementary file 5 (`supp5-v2.xlsx.txt`, sheet "Table 5 - ICICLE Oligos"); the
  bundle oligo scraper was run and its Table 5 hits followed to the original rows.

## Oligos (Supplementary file 5, verbatim)

- 🟢 MOSAIC_Bot `5'-/5phos/CTGTCTCTTATACACATCT-3'`
- 🟢 Tn5ME-s7_Top `5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-3'` ("Compatible with Nextera i7 primers")
- 🟢 Poly-A Top-L `5'-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVNAGATGTGTATAAGAGACAG-3'`
- 🟢 SI-P5-22 `5'-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3'`
- 🟢 F BC Primer `5'-CTACACGACGCTCTTCCGATCT-3'`
- 🟢 ADT-Rev-AMP `5'-CCTTGGCACCCGAGAATTCC-3'`
- 🟢 TruSeq R1 Seq Primer `5'-ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3'`
- 🟢 ADT-i7 Primer `5'-CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCC*A-3'`
  (PAGE; Methods: phosphorothioate bond linking the final two 3′ bases)

The selftest compares every complete base sequence with the extracted table.

## Construction (Methods and Fig. 3—S1)

- 🟢 Poly-T Top-L/MOSAIC_Bot and Tn5ME-s7_Top/MOSAIC_Bot are annealed separately
  (16 µM, 2× Dialysis Buffer, 95°C 5 min then 5°C per 2 min to 20°C), mixed 1:1 and
  complexed with Tn5 (Beta Lifescience TN5-BL01 + NaCl) at 25°C for 60 min.
- 🟢 PBMCs are blocked with TruStain FcX, stained with a 46-antibody TotalSeq-A panel
  (Supplementary file 6), washed three times and permeabilized with 0.01% digitonin.
- 🟢 20,000 cells in TD1 are tagmented with 6 µL complexes, 37°C 120 min. Fig. 3—S1a:
  about half of the fragments carry one Poly-T and one R2N end and can be barcoded and
  indexed; the model keeps same-end fragments in the mixture, not as a separated product.
- 🟢 GEMs: Chip G, Chromium Single Cell 3′ v3.1 Gel Beads, barcoding mix of NEBNext
  Ultra II Q5, Reducing Agent B, F BC Primer and ADT-Rev-AMP (0.2 µM each); 72°C 5 min
  gap fill, then 98°C 30 s and 12 cycles 98/42/65°C. ATAC is copied linearly with F BC
  Primer; ADTs are amplified exponentially with ADT-Rev-AMP (Fig. 3—S1c).
- 🟢 Recovery Agent, SILANE cleanup, 0.6×/2.0× SPRIselect; split 3:1 ATAC:ADT.
- 🟢 ATAC PCR: SI-P5-22 + Chromium i7 Multiplex Kit N Set A, Amp Mix, 7 cycles;
  0.4×/1.2× SPRI. ADT PCR: SI-P5-22 + ADT-i7, KAPA HiFi, 15 cycles; 1.2× SPRI.
- 🟡 The 10x 3′ v3.1 bead oligo is not printed: partial Read 1 (implied by F BC
  Primer), 16-base cell barcode (v3 whitelist) and 12-base UMI (28-cycle Read 1 − 16),
  then oligo(dT)30VN. It is drawn as inferred. TotalSeq-A barcode length (15) and the
  Chromium i7 Kit N primer sequence are vendor values outside this pass.

## Final libraries and reads

- 🟢 NovaSeq 28 bp Read 1 (cell barcode and UMI), 8 bp i7, 100 bp Read 2 (ATAC or ADT);
  the Table 5 TruSeq R1 Seq Primer is spiked in (0.3 µM) as Custom Read 1 primer.
- 🟢 Fig. 3—S1e: ATAC i7 and Read 2 use the Nextera Read 2 priming site (single-end
  gDNA read); ADT uses a "TruSeq small RNA/i7 Read 2 Priming site". The ADT loaded
  primer sequence is not printed, so its I1/R2 binding sites are not asserted.
- 🟡 ATAC library (derived): P5–TruSeq Read 1–CB16–UMI12–T30–VN–ME–gDNA–ME′–s7′–i7–P7.
  ADT library (derived): P5–TruSeq Read 1–CB16–UMI12–T30–VN–antibody barcode′–ADT-i7 arm′–i7–P7.

## Model corrections (2026-10-17)

The earlier model reused shared 10x scATAC/targeted-tagmentation and generic
antibody-tag endpoints: it placed the cell barcode beside an S5 handle with an invented
8-base assay barcode, used a 10-base UMI and Nextera Read 1. These were repository
modelling errors; the libraries are now declared locally from Table 5 and Fig. 3—S1.

## Open questions

- Tn5 complexing follows Mulqueen et al. 2019 (cited, not re-read).
- 10x 3′ v3.1 user guide (bead oligo, Chromium i7 Kit N sequences) not read in this pass.
- Exact ADT Read 2/Index 1 sequencing primer; TotalSeq-A member sequences.
