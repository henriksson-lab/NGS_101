# SNuBar-ATAC — research notes

Defining source: Wang *et al.*, [10.1016/j.molcel.2021.09.026](https://doi.org/10.1016/j.molcel.2021.09.026).

## Evidence used

- 🟢 Each sample is first tagmented with Illumina TDE1 and then receives an SBO oligo before nuclei are pooled for 10x single-cell ATAC (original STAR Methods, JATS lines 185–191).
- 🟢 The oligo contains a sequence complementary to the transposome universal tail, a sample barcode, UMI and PCR handle.
- 🟢 The sample tag is molecularly distinct from the later 10x droplet cell barcode.
- 🟢 The [open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC8611914/) identifies the SNuBar-ATAC barcode oligo (`SBO`) as a CST arm complementary to the Tn5 transposome universal tail, a sample barcode plus UMI, and a PCR handle. Its separate ATAC and RNA co-assay oligo (`SBOC`) adds a poly(A) capture tail; that tail is not evidence for the ATAC-only SBO.
- 🟢 The original article's STAR Methods (JATS lines 187–193) print the ATAC-only SBO family as `CCTTGGCACCCGAGAATTCCANNNNNNNN-(N)18-CTACTCTGCGTTGATACCAGACGCTGCCGACGA`: an eight-base UMI, eighteen-base sample-barcode window and source-known CST arm. They also print the TruSeq RPIX barcode-library primer `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`, P5 primer `AATGATACGGCGACCACCGAGATCTACAC`, custom Read 1 and i5 primers. The oligo family and two PCR primers are displayed on the public page.
- 🟢 After GEM cleanup the ATAC and SBO barcode templates are separated, PCR-amplified as distinct libraries, and mixed 8:2 for sequencing (JATS lines 191–193). The public final-library panel shows the ATAC branch; the reaction graph also shows the separate SBO library.
- 🟢 Table S1 and Methods S1 are gated on the PMC download route but were obtained from the publisher (Document S1, see the 2026-10-17 follow-up below). Table S1 lists 96 individual eighteen-base sample barcodes; the generic model keeps the N18 window. The author code repository is for demultiplexing and does not itself establish synthesis modifications.

## Bounded pass 2026-10-17 (round breadth-first-19)

Inspected: PMC8611914 JATS Results line 85 (SBO parts, delivery), STAR Methods lines 185–193 (tagmentation/SBO, 10x loading, AMPure separation, ATAC and barcode PCRs, 8:2 mix, spiked primers), line 205 (FASTQ layout), key resources lines 769–777.

- 🟢 Barcode-library reads (line 205): the 16-bp cell barcode comes from the Index 2 FASTQ; the first 18 bp of Read 1 are the sample barcode, followed by the 8-bp UMI. The ATAC library uses a single "Chromium i7 Sample Index N, Set A" with SI-PCR Primer B (line 193).
- 🟢 Templates are separated by 1.2X AMPure after Dynabead cleanup: ATAC template on beads, barcode supernatant raised to 1.8X (line 191).
- 🟡 Model correction: both final libraries were previously generic shells — the ATAC branch had an invented 8-nt i5 sample index with the cell barcode inside the insert, and the SBO branch a TruSeq Read 1/dA/Read 2 dual-index library. The SBO library is now assembled from the printed oligos: P5 – 16-bp cell barcode – CST complement (`TCGTCGGCAGCGTCTGGTATCAACGCAGAGTAG`, printed as the spiked custom Read 1 primer, line 193) – N18 sample barcode – N8 UMI – handle complement – RPIX arm – 6-nt i7 – P7. The spiked I5 index primer equals the SBO CST arm, so it reads the cell barcode. The ATAC branch places the cell barcode in the i5 position.
- 🟡 Derived: only the 3′-terminal 14 nt of the CST arm ({{= nextera.S5_RC }}) are complementary to the Nextera S5 tail carried by TDE1/Tn5-ME-Ap; the remaining CST bases are not complementary to the printed transposome adapters. The article's "complementary to the transposome universal tail" therefore refers to that 14-nt overlap.
- 🟡 Inferred, not inspected: the cell-barcode junctions with P5 and S5/CST (10x gel-bead oligo, guides CG000168 Rev A / CG000209 Rev D) and the ATAC i7 index length. Not found in the inspected Methods: the SBO library's Index 1 primer and acquired cycle counts. *(Superseded by the follow-up below: junctions, i7 length and run lengths are now source-stated; the RPIX arm carried a transcription error.)*

## Bounded follow-up 2026-10-17 (round deferred-followup-02)

Obtained from the Cell Press CDN (the article is free to read at the publisher):
Document S1 `1-s2.0-S1097276521007954-mmc1.pdf` (Figures S1–S5, Tables S1–S2,
Methods S1), Document S2 `…-mmc3.pdf` (version of record plus supplements) and Table S3
`…-mmc2.csv` (chromVAR motif results; no chemistry). Inspected: Document S1 Figure S1A
(text lines 28–88), Table S1 (lines 1896–2295), Methods S1 SNuBar-ATAC (lines
2561–2777); Document S2 STAR Methods lines 1672–1698 and 1750–1751.

- 🟢 Figure S1A prints the 10x cell-barcode primer as
  `AATGATACGGCGACCACCGAGATCTACAC-N16-TCGTCGGCAGCGTC`, the ATAC i7 primer as
  `CAAGCAGAAGACGGCATACGAGAT-I8-GTCTCGTGGGCTCGG` ("-I8-: 8bp index") and the RPIX primer
  with `-I6-`. The cell-barcode junctions and the 8-bp ATAC i7 index are therefore no
  longer drawn as inferred; the 10x guides themselves remain unread.
- 🟢 Methods S1 step 10: libraries mixed 4:1–20:1 (STAR Methods: 8:2) for NextSeq; custom
  Read 1 and I5 index primers spiked; "Read 1: 50bp, Read 2: 50bp, Index 1: 8bp, Index 2:
  16bp" — now the page's declared run lengths. Step 7: ATAC SI-PCR 9–12 cycles then
  0.4X/1.2X; barcode PCR with KAPA HiFi HotStart ReadyMix, RPIX and P5, 8–10 cycles, then
  1.5X SPRIselect. GEM incubation: 72 °C 5 min, 98 °C 30 s, 12 cycles of 98/59/72 °C.
- 🟢 Table S1 caption: each SBO is a 21-bp P7-PCR handle, 8-bp UMI, 18-bp sample barcode,
  19-bp universal sequence and 14-bp CST. All 96 printed rows have that 80-nt layout. The
  14-bp CST is the reverse complement of the Nextera S5 tail ({{= nextera.S5 }}), which
  settles the earlier derived 14-nt overlap; the model now splits the tail accordingly.
- 🟡 Model correction: the RPIX inner arm was transcribed with an extra C
  (`GTGACTGGAGTTC` before the handle); the printed primer has `GTGACTGGAGTT`
  (JATS line 771, Methods S1, Document S2 line 1694). The arm is now derived from the
  printed primer by construction.
- 🟢/🟡 Published discrepancy, shown on the page: Figure S1A's final barcode library
  prints `…CCAAGGTTGAGGTCAGTG-I6-…`, the 12 bases after the handle written in reverse
  relative to the RPIX primer printed in the same figure, STAR Methods and Methods S1
  (which give `…CCAAGGAACTCCAGTCAC-I6-…` on that strand). The model uses the primer.
- 🔴 Not named in Methods S1: an Index 1 or Read 2 primer for the SBO library (its RPIX
  index is 6 nt while Index 1 runs 8 cycles); those reads are not drawn on it.
- Deferred: 10x CG000168 Rev A / CG000209 Rev D, SNuBar-ARC branch, author code, a second
  independent pass on Documents S1/S2.

## Protocol boundary

SNuBar-ATAC adds pre-pooling sample identity to the standard 10x ATAC library; the parent 10x chemistry remains separately documented.
