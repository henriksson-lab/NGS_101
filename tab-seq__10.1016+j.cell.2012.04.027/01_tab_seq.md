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

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC3589129 JATS70–74, TAB-Seq of model DNA and specific loci.

Corrected pre-PCR bisulfite state: C is U, oxidized5mC is5caU and protected5hmC remains5gmC; C/T belongs to the later amplified sequencing readout. Replaced unsupported source-wide adapter/primer absence claims with unverified model geometry.

Deferred: Original Extended Experimental Procedures and WGBS precursor, full enzyme/adapter order, spike-in controls and read acquisition; later published protocol scraper hits not promoted as defining2012 oligos. This capped pass does not establish a complete source review or publication-wide absence.
