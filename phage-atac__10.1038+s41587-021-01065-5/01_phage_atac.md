# PHAGE-ATAC

**Defining source:** doi:10.1038/s41587-021-01065-5.

## Molecular boundary

🟢 Nanobody-displaying M13 phage carry a PAC/RD1 capture sequence; 10x ATAC droplet barcoding labels phagemids and tagmented chromatin together.

## Evidence boundary

🟢 Separate protein-derived-tag and ATAC libraries retain the same cell barcode; mitochondrial ATAC fragments can provide native clonal markers.

## Published oligos and remaining boundary

🟢 The linked bioRxiv preprint (doi:10.1101/2020.10.01.322420) supplies four
workbooks. The `media-3.xlsx` workbook lists the EF primer sequences
for PAC-tag construction and protein-derived-tag (PDT) library preparation;
`media-4.xlsx` lists the nanobody/hashtag readout
sequences. The methods explicitly use EF147 and EF91 for PDT-specific PCR, then
10x SI-PCR Primer B and an i7-index primer for final indexing (preprint extracted
text, lines 2283–2320). The PCR oligos as written in Table 1 are:

| PDT PCR oligo | 5′→3′ sequence |
|---|---|
| EF147 | AATGATACGGCGACCACCGAGA |
| EF91 | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGgataccgcggtgtattattgc |

🟢 The lower-case tract of EF91 is lower-case in the source workbook and is the
nanobody-facing annealing portion; it is not an inferred modification. The paper
reports paired-end 34 + 34 sequencing for the pooled PDT and ATAC libraries.

🟡 The present page represents those library branches by role and omits the
source-published EF primer bases and individual nanobody tags. The exact 10x
gel-bead oligo, SI-PCR Primer B and chosen i7 primer require an authoritative
version-matched 10x specification; the paper gives their names but no exact
sequences for its run. The broad page limitation therefore mixes a **modeling
omission** (EF147/EF91 and readout tags) with a **vendor-source gap** (the
specific 10x kit oligos).
