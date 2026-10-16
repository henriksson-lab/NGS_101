# G&T-seq — research notes

## Sources read

- 🟢 Macaulay et al. 2015, DOI [10.1038/nmeth.3370](https://doi.org/10.1038/nmeth.3370), [author-uploaded original PDF with Online Methods](https://www.researchgate.net/profile/Mubeen-Goolam/publication/275525502_GT-seq_Parallel_sequencing_of_single-cell_genomes_and_transcriptomes/links/5757e80f08ae5c65490724bc/G-T-seq-Parallel-sequencing-of-single-cell-genomes-and-transcriptomes.pdf), PDF lines 554–623. The supplement workbooks and figures were also checked; the large MOESM220 workbook is Fig. 1b/1c measurement data, not a laboratory protocol. The PDF is readable through the web source viewer, though direct download into the local source bundle failed.
- 🟢 Macaulay et al. 2016, [authors' detailed G&T-seq protocol](https://www.pure.ed.ac.uk/ws/files/181359881/Separation_and_parallel_sequencing_of_the_genomes_and_transcriptomes_of_single_cells_using_G_T_seq.pdf), DOI 10.1038/nprot.2016.138: oligonucleotide list at extracted lines 644–658 and capture-bead preparation at lines 718–734.

## Molecular order

🟢 The original 2015 Online Methods specify this order: a lysed single cell is mixed with a biotinylated oligo(dT) primer. Poly(A) RNA is
captured on streptavidin beads and magnetically separated from the genomic-DNA-containing
supernatant. The bead fraction undergoes Smart-seq2-like reverse transcription,
template-switching and WTA. The DNA fraction undergoes either MDA or PicoPLEX WGA. The
two amplified fractions are converted into separate indexed sequencing libraries.

## Protocol boundary

MDA and PicoPLEX are alternative genome-amplification choices tested by the paper. The
transcriptome and genome products never become one chimeric construct. Their common cell
identity is inherited from physical separation of a single-cell lysate.

🟢 The 2015 Online Methods print the capture/RT primer as
`5′-biotin-triethyleneglycol-AAGCAGTGGTATCAACGCAGAGTAC-T30-VN-3′` (PDF lines 573–574),
the TSO `5′-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3′` (lines 597–599), and the
PCR primer `5′-AAGCAGTGGTATCAACGCAGAGT-3′` (lines 604–605). The method names
MDA/PicoPLEX as alternative WGA branches (lines 612–616) and Nextera XT for
the cDNA and gDNA libraries (lines 621–623).

🟢 The authors' 2016 detailed protocol separately prints the capture/RT primer verbatim as
`5′-Biotin-TEG-AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′`, with a 30-base dT tract and
3′ `VN` anchor (PDF text lines 644–651). It also prints the TSO
`5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′` and ISPCR
`5′-AAGCAGTGGTATCAACGCAGAGT-3′` (lines 652–658). The same protocol names
Nextera XT library preparation and PicoPLEX or MDA as alternative WGA branches
(lines 171–185, 421–440, 1373–1387).

The 2015 TSO ends `…AGTACrGrG+G`; the 2016 author protocol adds `AT` before
`rGrG+G`, giving `…AGTACATrGrG+G`. This is a **source-version difference**: the
public 2015 page uses the defining paper's shorter TSO. Whether the 2016 change
was intentional or a source error is unresolved. Kit-specific Nextera XT index
choices and WGA reagent oligos remain unresolved for a selected branch; the
public page shows generic final-library geometry for those branches.
