# G&T-seq — research notes

## Sources read

- 🟢 Macaulay et al. 2015, DOI [10.1038/nmeth.3370](https://doi.org/10.1038/nmeth.3370), main article and supplementary protocol/figures.

## Molecular order

🟢 A lysed single cell is mixed with a biotinylated oligo(dT) primer. Poly(A) RNA is
captured on streptavidin beads and magnetically separated from the genomic-DNA-containing
supernatant. The bead fraction undergoes Smart-seq2-like reverse transcription,
template-switching and WTA. The DNA fraction undergoes either MDA or PicoPLEX WGA. The
two amplified fractions are converted into separate indexed sequencing libraries.

## Protocol boundary

MDA and PicoPLEX are alternative genome-amplification choices tested by the paper. The
transcriptome and genome products never become one chimeric construct. Their common cell
identity is inherited from physical separation of a single-cell lysate.

🟡 The paper's detailed oligo supplement is accessible but the public page keeps the
downstream Smart-seq2 and Nextera arms at their reusable canonical geometry; it focuses on
the G&T-specific molecular operation, biotin-oligo(dT) partitioning.
