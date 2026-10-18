# TAB-seq — TET-assisted bisulfite sequencing

## Sources read

- 🟢 Yu *et al.*, “Base-resolution analysis of 5-hydroxymethylcytosine in the mammalian genome,” *Cell* (2012), DOI `10.1016/j.cell.2012.04.027` (defining paper).
- 🟢 Yu *et al.*, “Tet-assisted bisulfite sequencing of 5-hydroxymethylcytosine,” *Nature Protocols* (2012), DOI `10.1038/nprot.2012.137` (open author manuscript, PMCID: PMC3641661).

## Molecular interpretation

- 🟢 Sonicated genomic DNA receives C, 5mC and 5hmC conversion controls.
- 🟢 T4 β-glucosyltransferase protects 5hmC as 5gmC. Excess recombinant mTet1 then oxidizes 5mC to 5caC while 5gmC remains protected.
- 🟢 Bisulfite converts unmodified C and 5caC to U. After amplification, original C and 5mC read as T, whereas original 5hmC reads as C.
- 🟢 The companion protocol delegates whole-genome library construction to historical Illumina guide 15021861_A; its native adapter architecture is unassigned in this model. The protocol requires measuring protection and oxidation efficiencies from spike-ins.

## Final readout

TAB-seq directly reports 5hmC as retained C. It has the opposite biological target from oxBS-seq, which retains 5mC as C.

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC3589129 JATS70–74, TAB-Seq of model DNA and specific loci.

Historical model correction (superseded below): this pass labelled oxidized 5mC as 5caU. The current physical bisulfite state is U, U and protected 5gmC; C/T belongs to the later amplified sequencing readout. Replaced unsupported source-wide adapter/primer absence claims with unverified model geometry.

Deferred: Original Extended Experimental Procedures and WGBS precursor, full enzyme/adapter order, spike-in controls and read acquisition; later published protocol scraper hits not promoted as defining2012 oligos. This capped pass does not establish a complete source review or publication-wide absence.


## Bounded native migration (2026-10-18, batch07)

Inspected 2 originals: Defining PMC3589129 ExperimentalProcedures144–146; published Nature Protocols2012 PMC3641661 overview37–44, steps15–18/text484–520, whole-genome25B/text675–688 and publication metadata. Mandatory bundle scraper reviewed.

Replaced hybrid/default library recipes with local source-owned graph. Primary glucosylation60min and oxidation90min, PK/BioSpin30/QIAquick PCRPurification are explicit; companion80min and NucleotideRemoval are separately scoped ordinary protocol differences. Correct physical bisulfite substrate U/U/5gmC feeds PCR; normal daughter T/T/C and retained modified originals are distinct. No canonical TruSeq/index/primer shell or locus-specific PfuTurboCx enzyme is assigned to whole-genome libraries. Companion publication is2012;2013 is later PMC deposition.

Deferred: Historical Illumina WGBS guide15021861_A explicitly linked by25B, not inspected; native endprep/adapter order, exact amplification reagents and sequencing acquisition remain bounded; Defining original extended procedures, complete control repertoire and linked assay details beyond inspected sections. This capped pass does not establish complete review or source-wide absence.
