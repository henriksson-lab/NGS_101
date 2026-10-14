# scM&T-seq — research notes

## Sources read

- 🟢 Angermueller et al. 2016, DOI `10.1038/nmeth.3728`, PMC full text and supplementary protocol.
- 🟢 SMART-seq2 and scBS-seq defining papers and supplements for the two downstream branches.

## Construct path

🟢 Biotinylated oligo(dT) beads capture cellular mRNA. The bead-bound fraction is processed by Smart-seq2 and Nextera XT; the supernatant and washes retain genomic DNA for scBS-seq. 🟢 The common cell identity is plate provenance, not a covalent shared barcode, so the schematic does not invent one.

🟢 The schematic reuses the oligos from the branch-defining sources cited by
scM&T-seq: Smart-seq2 oligo-dT30VN, TSO and ISPCR; canonical Nextera N5xx/N7xx landing
regions; and scBS-seq oligo1/oligo2 handles with PE1/indexed iPCRTag. The scM&T paper's
Methods explicitly change both bisulfite strand-synthesis primers to **3′ random
hexamers**. The original scBS-seq oligo1/oligo2 use N9; those N9 regions must not be
inherited into scM&T. The DNA branch now models N6 on each primer, while retaining the
source-printed scBS-seq handles.
