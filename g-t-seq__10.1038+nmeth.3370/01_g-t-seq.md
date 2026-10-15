# G&T-seq — research notes

## Sources read

- 🟢 Macaulay et al. 2015, DOI [10.1038/nmeth.3370](https://doi.org/10.1038/nmeth.3370), public abstract and supplementary figures/results files. The full original Methods are not yet in the local source bundle.

## Molecular order

🟡 The accessible abstract and later G&T-seq implementations support this outline, pending the original Methods: a lysed single cell is mixed with a biotinylated oligo(dT) primer. Poly(A) RNA is
captured on streptavidin beads and magnetically separated from the genomic-DNA-containing
supernatant. The bead fraction undergoes Smart-seq2-like reverse transcription,
template-switching and WTA. The DNA fraction undergoes either MDA or PicoPLEX WGA. The
two amplified fractions are converted into separate indexed sequencing libraries.

## Protocol boundary

MDA and PicoPLEX are alternative genome-amplification choices tested by the paper. The
transcriptome and genome products never become one chimeric construct. Their common cell
identity is inherited from physical separation of a single-cell lysate.

🟡 The accessible Springer supplementary files contain figures and
analysis data rather than an ordering table. The defining oligo is the
5′-biotin-triethyleneglycol-linked capture/RT primer, reported in a later
[direct G&T-seq implementation](https://pmc.ncbi.nlm.nih.gov/articles/PMC7581540/) as
`5′-biotin-TEG-AAGCAGTGGTATCAACGCAGAGTAC-T30VN-3′`.
That exact transcription is from a later implementation, so it should not be
treated as a verified 2015 ordering record until the original Methods or the
authors' [2016 detailed protocol](https://doi.org/10.1038/nprot.2016.138) is read in full. The downstream Smart-seq2
TSO/WTA and MDA-or-PicoPLEX WGA branches also need branch-specific oligo
verification; the current public page uses reusable canonical geometry.
