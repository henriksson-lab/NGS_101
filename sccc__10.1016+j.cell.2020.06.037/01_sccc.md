# scCC

**Defining paper:** Moudgil et al., *Cell* (2020), doi:10.1016/j.cell.2020.06.037. The open author manuscript and methods were read from PMCID PMC7510185.

## Molecular boundary

🟢 A self-reporting piggyBac transposon lacks its downstream polyadenylation signal, so reporter transcription continues through the terminal repeat into the genomic insertion flank.

🟢 The 10x first-strand product is split. The scCC branch is specifically amplified and circularized, moving the transcript's cell barcode and UMI next to the transposon–genome junction before linear-DNA removal, shearing and junction library construction.

## Evidence boundary

🟢 The paper's Table S4 defines custom adapters and primers. The page preserves their roles and links the source rather than treating this historical Nextera Mate Pair implementation as a universal kit geometry.

## Oligo source audit

🟢 The author preprint's supplementary PDF (`538553_media-1.pdf.txt`, Supplemental Table 4) is retrievable even though the published Cell supplement `NIHMS1608211-supplement-1.pdf` is gated at PMC. It prints the scCC library oligos verbatim, including modifications and a variable sample index. The two annealed adapter components are `scCC_P5_adapter`, `AATGATACGGCGACCACCGAGATCTTCACTCATTCCACACGACTCCTTGCCAGTCTC*T`, and `scCC_P7_adapter`, `/5Phos/GAGACTGGCAAGTACACGTCGCACTCACCATGA[index]ATCTCGTATGCCGTCTTCTGCTTG`. Their paired geometry should be computed before an exact junction is shown.

🟢 The same table gives the final PCR primers `scCC_P5_primer`, `AATGATACGGCGACCACCGAGATC`, and `scCC_P7_primer`, `CAAGCAGAAGACGGCATACGAGAT`; the custom sequencing primers are `scCC_PB_CustomRead2`, `CGTGTAGGGAAAGAGTGTGCGTCAATTTTACGCAGACTATCTTTCTAG`, and `scCC_CustomIndex1`, `GAGACTGGCAAGTACACGTCGCACTCACCATGA`. Capture primers include `Bio_Illumina_Seq1_scCC_10X_3xPT`, `/5Phos/ACACTCTTTCCC/iBiodT/ACACGACGCTCTTCCGA*T*C*T`, and `Bio_Long_PB_LTR_3xPT`, `/5Phos/GCGTCAATTTTACGCAGAC/iBiodT/ATCTTTC*T*A*G`.

🟡 These recovered bases belong to the author's preprint supplement. The final Cell supplement could not be compared byte-for-byte through PMC; retain the preprint provenance on any model transcription and leave the variable `[index]` window symbolic.
