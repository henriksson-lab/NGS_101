# CITE-seq — research notes

## 1. What it is

CITE-seq converts surface-protein abundance into sequenceable antibody-derived tags while measuring the transcriptome from the same single cell.

## 2. Sources read

- 🟢 Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, Nature Methods 14, 865–868 (2017), DOI [10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380), open manuscript [PMC5669064](https://pmc.ncbi.nlm.nih.gov/articles/PMC5669064/).
- 🟢 Main Methods, especially the **Antibody-oligo sequences** subsection, and [author supplementary protocol, pp. 6–7](https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4380/MediaObjects/41592_2017_BFnmeth4380_MOESM3_ESM.pdf), establish the poly(A)-tailed antibody tag, DNA-dependent extension during RT, size separation and independent ADT library PCR. The protocol has an additional v2_BC7 oligo absent from the paper's listed CBMC subset.

## 3. Published antibody oligos

🟢 For the 10x Single Cell 3′ **v2** CBMC experiment, the author protocol prints 14 ordered oligos with `/5AmMC12/` at the 5′ end, the common `CCTTGGCACCCGAGAATTCC` TruSeq small-RNA Read-2 handle, a shared `A` followed by a **six-base** antibody identifier, one degenerate `B`, and a poly(A) tail. The seven-base `A` + identifier regions in BC1–BC14 order are `AATCACG`, `ACGATGT`, `ATTAGGC`, `ATGACCA`, `AACAGTG`, `AGCCAAT`, `ACAGATC`, `AACTTGA`, `AGATCAG`, `ATAGCTT`, `AGGCTAC`, `ACTTGTA`, `AAGTCAA`, and `AAGTTCC`. All 14 lengths were checked against the printed sequence rows; the defining paper prints 13 of these.

🟢 A complete example (v2_BC1) is `/5AmMC12/CCTTGGCACCCGAGAATTCCAATCACGBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA`. The antibody identifier is the six bases `ATCACG` in Supplementary Table 2 (PMC supplement-1, lines 45–55); the preceding `A` is shared by the v2_BC1–14 family (Springer author protocol, lines 403–434). The terminal `B` and poly(A) tract are retained as ordered; they are not antibody identifiers. The paper separately prints Drop-seq and 10x 3′ v1 species-mixing designs with Nextera Read-2 and Read-1 handles, respectively. Those are different branches of CITE-seq, not replacements for the v2 CBMC construct.

🟢 The author protocol also prints both ADT library PCR oligos for the 10x v2 branch: PE 1.0 P5 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T` and Small RNA RPI1 `CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`. The two `*` positions in P5 are terminal phosphorothioate bonds. It separately gives a different P5-SMART-PCR hybrid oligo for the Drop-seq branch.

🟢 The paper's Online Methods (lines 94–132 of the PMC JATS text) name the two species-mixing oligos as BC6 and BC12 for both Drop-seq and 10x 3′ v1, with different platform-specific PCR handles. Its Supplementary Table 2 (PMC supplement-1, lines 178–196) assigns `GCCAAT`/BC6 to human CD29 and `CTTGTA`/BC12 to mouse CD29. These assignments apply to the species-mixing experiments, not to the 10x v2 CBMC panel.

🟢 The paper calls its 10x v2 CBMC assay a **13-antibody panel** (JATS line 70), prints the v2 family under “CBMC profiling” (lines 106–132), and maps the antibodies to identifiers in Supplementary Table 2 (lines 45–176). Crosswalking the table to the author's full v2 sequence list (protocol lines 403–434) resolves the panel as CD3e→v2_BC1, CD19→BC2, CD16→BC3, CD4→BC4, CD45RA→BC5, CD8a→BC6, CCR5→BC8, CCR7→BC9, CD11c→BC10, CD34→BC11, CD56→BC12, CD14→BC13, and CD10→BC14. The author's v2_BC7 (`ACAGATC`) is an additional printed oligo; Table 2 assigns no panel antibody to this identifier. The paper's 13 listed v2 oligos likewise omit BC7. This crosswalk identifies the CBMC panel without assuming that every oligo in the author protocol was used.

🟡 The four CD8a FACS-sort comparison pools used an oligo-labelled CD8a antibody and Drop-seq (JATS lines 68 and 156), and Table 2 gives the CD8a identifier `GCCAAT`. The Methods do not name the complete handle-bearing oligo separately for each pool. Assigning the CBMC v2_BC6 construct to those Drop-seq validation runs would therefore be an inference.

## 4. Construct-changing steps

1. Conjugate antibody identity oligos containing PCR handle, antibody barcode and poly(A) tail.
2. Label cells and wash away free tags.
3. Co-encapsulate cells and barcoded beads; oligo-dT captures both mRNA and ADT poly(A).
4. Reverse transcription/extension gives both molecule classes the cell barcode and UMI.
5. Size-separate large cDNA from short ADTs and make separate libraries.

## 5. Evidence boundary

🟢 The defining paper and author protocol print the antibody oligos and ADT PCR primers. The v2 handle is transcribed into the molecular model with its shared `A` followed by a six-base antibody identifier from Supplementary Table 2. 🟡 The final page uses role placeholders for panel-specific antibody barcodes and normalizes platform arms for computed primer placement. The complete 10x bead oligo is a separate vendor dependency.
