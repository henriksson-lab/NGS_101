# scCAT-seq — source intake

Defining paper: Liu et al., *Nature Communications* (2019), DOI
[10.1038/s41467-018-08205-7](https://doi.org/10.1038/s41467-018-08205-7).

🟢 The protocol begins by physically separating cytoplasmic RNA and the intact nucleus
from one cell, followed by Smart-seq2 and scATAC-seq library preparation.

🟢 The defining paper's Methods give the RNA-branch oligos verbatim: modified
oligo-dT `AAGCAGTGGTATCAACGCAGAGTACT30VN`, template-switching oligo
`AAGCAGTGGTATCAACGCAGAGTACATrGrG+G`, and PCR primer
`AAGCAGTGGTATCAACGCAGAGT`. These are the canonical SMART-seq2 handle system.

🟢 **The TSO linker is `ACAT` on this paper's own authority.** Liu et al. 2019 Methods
(Nature Communications, PMC6349937; the same sentence is in the bioRxiv preprint 316208)
print `5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′` — 23-nt SMART handle, `ACAT`, `rGrG+G`.
That distinction matters because SMART-seq2 published two linkers: Picelli et al. 2013
Supplementary Table 2 gives the chosen `rGrG+G` oligo as handle + `AC` + `rGrG+G`, while
the *Nature Protocols* form carries `ACAT`. scCAT-seq states which one it used, so
`lib/smartseq2.TSO_LINKER` is reproducing a published base here, not an assumption;
`tools/selftest.py` pins the shared constant to the sentence above.

🟢 The same Methods give both ATAC preamplification primers:
`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` and
`GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`.

🟢 The final chromatin PCR instead names only a “universal primer” and a
“barcode primer” without their bases. The RNA branch refers onward to SMART-seq2 for
library preparation, and both products undergo a further BGISEQ-500 ssDNA-circle/DNB
conversion whose oligos are likewise not listed. Those endpoint sequences are a real
remaining source gap; they should not obscure the exact upstream oligos above.
