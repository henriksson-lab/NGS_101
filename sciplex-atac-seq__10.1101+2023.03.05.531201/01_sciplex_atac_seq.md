# sciPlex-ATAC-seq

**Defining preprint:** Booth et al. (2023), doi:10.1101/2023.03.05.531201; published in *BMC Genomics* (2023), doi:10.1186/s12864-023-09832-1. The bioRxiv v1 PDF and Supplementary Tables 3–4, published JATS and author [two-level protocol](https://www.protocols.io/view/sciplex-atac-2-level-cz56x89e.pdf) were read.

## Molecular boundary

🟢 Permeabilized nuclei absorb unmodified, sample-specific single-stranded DNA hash oligos. An indexed primer-extension step makes polyadenylated hashes recoverable in the ATAC workflow.

🟢 Two-level combinatorial indexing gives each nucleus an ordered cell identity independent of its sample hash. Hash reads and ATAC reads are sequenced separately in the two-level implementation; the later three-level implementation handles them differently.

## Evidence boundary

🟢 The hash is modelled as a sample index, never as a cell barcode. The two combinatorial well indices jointly identify a nucleus.

🟢 The defining preprint Methods (local `2023.03.05.531201_v1.full.pdf.txt`, lines 649–673) prints a hash oligo with a 10-base sample ID and A30 tail, and a two-level capture primer with N8 UMI, 10-base well ID and T30VN. Its Supplementary Table 3 (`media-2.xlsx.txt`, rows 1–4) supplies sample hash IDs; Table 4 (`media-3.xlsx.txt`, rows 1–22) prints the eight indexed P5 and twelve indexed P7 Tn5 oligos. The page displays one P5 example and the original hash/capture family.

🟢 **Source-version difference:** the preprint's two-level capture oligo (Methods lines 667–673) and author two-level protocols.io PDF (pp. 3–4, lines 119–125) print an unmodified poly(T)30 tract. The published *BMC Genomics* JATS (local `s12864-023-09832-1_PMC10696879.xml.txt`, line 179) prints `TTT + TTT + …` and defines `+` as a locked nucleic acid modification. The page follows the defining preprint's plain-T form, while its graph represents only the capture and indexing roles. The available sources do not establish which version was used for each published sample; neither version is silently substituted for the other.

🔴 The defining preprint and published article also describe a distinct three-level sciPlex-ATAC3 implementation. Its linked detailed protocol and all version-specific barcode choices require separate review before an exact three-level molecule is drawn.
