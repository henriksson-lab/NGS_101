# SNARE-seq — research notes

## Sources read

- 🟢 Chen, Lake & Zhang, *Nature Biotechnology* 37:1452–1457 (2019), DOI
  `10.1038/s41587-019-0290-0`, read as the PMC6893138 author-manuscript JATS (Methods:
  tagmentation line 60, nuclei barcoding 64, library preparation 68, NGS sequencing 72,
  preprocessing 76), with Supplementary Table 4 (`sup_tab4.xlsx`, identical to Springer
  MOESM6) and Table S3 (MOESM5).
- The other publisher supplements are a figure/legend file (MOESM1, NIHMS docx) and the
  Reporting Summary (MOESM2, image PDF, inspected rendered). There is no separate
  Supplementary Protocol; an earlier version of this note wrongly listed one.
- Drop-seq (Macosko et al. 2015, PMC4481139) is the cited predecessor for droplets, beads
  and the Nextera XT cDNA library; only its main text is on disk (Extended Experimental
  Procedures and Table S6 are `(manual)` in the manifest).
- The Teichlab SNARE-seq construct page was used for orientation only; it is not source evidence.

🟢 **Published table cross-reference error:** Chen et al. Methods (original PMC
JATS line 64) directs the Nextera-R1-rc-polyA splint to “Table S3”. The actual
Table S3 (`sup_tab3.xlsx.txt`, lines 1–4) contains differential-accessibility
results; it does not list oligos. Original Table S4 (`sup_tab4.xlsx.txt`, sheet
“Oligonucleotide sequences”, lines 1–27) prints
`/5Phos/CTGTCTCTTATACACATCTGACGCTGCCGACG` followed by 31 A (63 nt) for that splint,
plus blockers, TSO, PCR, index and custom Read 1 primers. The model uses the Table S4
strand. There is no disagreement about the printed bases; the error is the Methods
table number.

🟡 Model correction (2026-10-17): the first of those 31 A completes the S5 complement,
so the splint is a 33-nt Nextera arm plus a 30-nt poly(A). The earlier model added an
extra A (64 nt).

## Construct path

🟢 About 200,000 nuclei are tagmented with TDE1 (37 °C, 30 min) before droplets. The
splint is added to the Drop-seq lysis buffer; nuclei and barcoded beads (ChemGenes
Macosko201110) are co-encapsulated and heated to 72 °C for 5 min to lyse nuclei and
release Tn5. 🟡 The splint's 33-nt arm anneals across an S5+ME end and its poly(A)
binds the bead dT30, so a 12-base cell barcode is shared by both branches. 🟢 The
preprocessing assigns cell barcode and UMI from Read 1 for both RNA and chromatin data.

🟢 After droplet breakage, beads are washed in 6X SSC with the blocking oligos
Nextera-R1-bk and Nextera-Ad2-bk (both `/5Phos/ … /3InvdT/`); their role is not
explained. One RT/ligation mix (T7 ligation buffer, TSO, Hemo KlenTaq, T7 ligase,
reverse transcriptase) is incubated at room temperature for 30 min and 42 °C for
90 min, followed by Exonuclease I (37 °C, 45 min). 🟡 The captured chromatin is
gap-filled and the phosphorylated splint-side nick can be sealed; the other nick is not
shown as covalently sealed, because the ordered S5 strand has no stated 5′ phosphate.

🟢 10,000 beads are PCR-amplified with Nextera-R2 / Tso-PCR for 16 cycles in total and
column-purified; the product is split into two halves. cDNA half: 0.6X bead selection,
Tso-PCR-alone enrichment, 0.6X selection, then Nextera XT "as described previously".
Chromatin half: P5XX-Tso / Ad2.X index PCR, then PAGE excision of 225–1000 bp and
column purification.

## Final libraries and reads (Methods line 72)

- 🟢 Chromatin: P5–i5 (P50x-Tso)–`GCCTGTCCGCGG`–bead SMART handle–cell barcode–UMI–dT30–
  S5–ME–gDNA–ME–S7–i7 (Ad2.x)–P7 (derived from the Table S4 primers). HiSeq 2500:
  Read1CustSeqB Read 1 30 bp, HP10 Index 1 75 bp, HP11 Read 2 75 bp, plus an 8-bp Index 2.
  🟡 The pipeline's paired chromatin reads are "read2 and read3", so the 75-bp Index 1
  read is genomic; HP10/HP11 contents are not printed, and their Nextera Read 1/Read 2
  sites are inferred. The Index 2 priming site is not printed; the read table on the page
  is source-declared for this library.
- 🟢 cDNA: Read1CustSeqB Read 1 30 bp (barcode 1–12, UMI 13–20) and HP11 Read 2 80 bp; no
  index read is stated. 🟡 The Nextera XT P5-side primer and i7 adapter are not printed
  in the SNARE-seq sources; the page draws P5, the `GCCTGTCCGCGG` linker required by
  Read1CustSeqB, the ME/S7 end and i7 as inferred.

🟢 Table S4 prints six `P50x-Tso` i5 primers (with terminal phosphorothioates `T*A*C`),
six `Ad2.x` i7 primers and `Read1CustomSeqB`. Which index pair each library used is not
stated; P501-Tso and Ad2.1 are shown as representatives.

## Open questions

- Drop-seq Extended Experimental Procedures / Table S6 (bead sequence, Nextera XT custom
  P5 primer) are not on disk; the bead constant sequence and 8-T spacer are not re-verified.
- Purpose of the blocking oligos; the polymerase used for the PCRs; per-library index
  members; the HiSeq Index 2 priming route.
