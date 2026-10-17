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
`smart-seq2__10.1038+nmeth.2639/tools/smartseq2.TSO_LINKER` is reproducing a published base here, not an assumption;
`tools/selftest.py` pins the shared constant to the sentence above.

🟢 The same Methods give both ATAC preamplification primers:
`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` and
`GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`.

🟢 The final chromatin PCR instead names only a “universal primer” and a
“barcode primer” without their bases. The RNA branch refers onward to SMART-seq2 for
library preparation, and both products undergo a further BGISEQ-500 ssDNA-circle/DNB
conversion whose oligos are likewise not listed. Those endpoint sequences are a real
remaining source gap; they should not obscure the exact upstream oligos above.

## Bounded pass 2026-10-17 (round breadth-first-16)

Sources read: published PMC6349937 JATS `.xml.txt` Methods lines 201 (bulk ATAC), 219–229
(scCAT-seq separation, RT/TSO/ISPCR, nucleus tagmentation, carrier, ATAC preamp and final
PCR, RNA tagmentation, BGISEQ-500 conversion and sequencing) and reference list (ref 7
Picelli 2013 Smart-seq2; ref 18 ATAC predecessor; ref 36 BGISEQ-500 reference dataset);
Supplementary Information (MOESM1) Supplementary Figure 1 legend lines 290–293; MOESM3
supplementary-data descriptions (count matrices and statistics, not chemistry). Scraper:
21 hits in 3 files, all the Methods oligos above plus the two qPCR accessible-region primers.

🟢 Steps now modelled from Methods line 225–227: nucleus tagmentation 15 min at 37 °C;
EDTA stop at 50 °C for 15 min; 30 ng carrier plasmid DNA with RLT Plus, 15 min at 37 °C;
1:1.8 XP purification; NEBNext HF 2x mix with 20 µM adapter 1/2 primers, 72 °C 5 min
then 8 cycles; 1x XP; another 10–16 cycles (15 for one cell) with a “universal primer”
and a “barcode primer”; XP size selection for 150–700 bp; RNA tagmentation of 2 ng cDNA
at 55 °C for 10 min with 0.1% SDS release, then PCR and library preparation “according to
the Smart-seq2 method”.

🟢 Line 229: all libraries were converted to BGISEQ-500 ssDNA circles (2 pmol pool),
rolled into DNBs and sequenced paired-end 50 bp. 🟡 The previous model drew both
endpoints as Illumina P5/i5/S5–ME…ME–S7/i7/P7 libraries read by Nextera primers; that
was a model error, not a source claim. The page now shows only the Tn5-adapter-flanked
core (RNA branch adapters inferred, since Methods does not name the RNA-branch Tn5
adapters), a warned BGISEQ circle/DNB boundary, and a PE50 read table without primer
landing.

🔴/open: universal and barcode primer bases and any platform tails; RNA-branch PCR primers
(Smart-seq2 Nat Methods 2013 and its protocol not read in this pass); BGISEQ-500
circle/DNB oligos and sequencing-primer entry (ref 36 not read). These remain named
deferrals, not confirmed absences.
