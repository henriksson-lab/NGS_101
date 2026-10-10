# TAB-seq — TET-assisted bisulfite sequencing

## Sources read

- 🟢 Yu *et al.*, “Base-resolution analysis of 5-hydroxymethylcytosine in the mammalian genome,” *Cell* (2012), DOI `10.1016/j.cell.2012.04.027` (defining paper).
- 🟢 Yu *et al.*, “Tet-assisted bisulfite sequencing of 5-hydroxymethylcytosine,” *Nature Protocols* (2012), DOI `10.1038/nprot.2012.137` (open author manuscript, PMCID: PMC3641661).

## Molecular interpretation

- 🟢 Sonicated genomic DNA receives C, 5mC and 5hmC conversion controls.
- 🟢 T4 β-glucosyltransferase protects 5hmC as 5gmC. Excess recombinant mTet1 then oxidizes 5mC to 5caC while 5gmC remains protected.
- 🟢 Bisulfite converts unmodified C and 5caC to U. After amplification, original C and 5mC read as T, whereas original 5hmC reads as C.
- 🟢 Whole-genome libraries follow Illumina’s whole-genome bisulfite/TruSeq construction. The protocol requires measuring protection and oxidation efficiencies from spike-ins.

## Final readout

TAB-seq directly reports 5hmC as retained C. It has the opposite biological target from oxBS-seq, which retains 5mC as C.
