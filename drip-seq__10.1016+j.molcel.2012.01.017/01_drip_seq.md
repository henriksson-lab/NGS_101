# DRIP-seq — DNA:RNA immunoprecipitation sequencing

## Current source scope

- 🟢 Ginno *et al.*, “R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters,” *Molecular Cell* (2012), DOI `10.1016/j.molcel.2012.01.017`. The defining published article was read as PMC3319272 JATS, especially Experimental Procedures, “Genome-wide R-loop profiling using DRIP-seq and DRIVE-seq,” line 113.
- 🟢 The publisher's supplementary PDF **mmc1** is now available in the source bundle. This pass inspected **Figure S5 caption, lines 73–87**. S5A describes an **in-vitro validation**: transcribed R-loop substrate mixed with dsDNA, S9.6 capture on agarose A/G beads, washes and proteinase K/phenol–chloroform recovery. Those validation reagents are not assigned wholesale to genomic DRIP-seq.
- 🟢 The publisher's **mmc2 ZIP** is also available. Its original `molcel_4184_mmc2.xlsx` workbook was inspected **in memory**, without copying source material into the checkout. The sheet identities and first header rows describe all GC-skewed regions, promoters with positive GC skew, and gene ontology classifications, including Hg18 coordinate annotations. This is a limited workbook identity/header check, not a complete cell-by-cell inspection or an accepted source-absence finding.
- The older PMC NIHMS supplement endpoints were gated. That historical route failure does not mean the published supplements are unavailable: the publisher originals above are on disk. Complete mmc1 inspection, including its OCR-bearing page 8, remains deferred. Equivalence of the old PMC package names to individual publisher files is not assumed.

## Defining 2012 chemistry

- 🟢 Total nucleic acids were extracted from Ntera2 cells using SDS/proteinase K at 37 °C, phenol–chloroform and ethanol precipitation.
- 🟢 HindIII, EcoRI, BsrGI, XbaI and SspI fragment the DNA. RNase H pretreatment is a **parallel control**, not an RNA-removal operation applied to the experimental captured population.
- 🟢 Genomic DRIP uses S9.6 immunoprecipitation following the cited MeDIP procedure **without denaturation**. The MeDIP predecessor, Weber *et al.* (2005), has not been inspected here. The graph leaves its detailed genomic capture support/wash implementation unassigned.
- The model retains the **RNA:DNA hybrid and displaced DNA through immunoprecipitation**. Capturing a hybrid does not physically remove its RNA component. After sonication, the model records hybrid integrity and RNA fate as unresolved.
- 🟢 The pulled-down material, RNase H control and input DNA were **sonicated, size-selected and ligated to Illumina barcoded adaptors**, then sequenced on **GAIIx and HiSeq**. No numeric size window, adapter geometry, primer assignment or acquired read/index mask is supplied by this model.
- The transition from hybrid-associated material to the DNA library substrate is an explicit unresolved boundary. It does not select an RNA-destruction enzyme, an end-repair/dA-tailing route or a library PCR program. The adaptor-ligation state has unavailable physical end geometry.

## Protocol boundaries and model corrections

The authors' Sanz and Chédin *Nature Protocols* (2019) method, DOI `10.1038/s41596-019-0159-1`, remains a **separate later protocol**. Earlier notes had inspected that method, but the preceding model incorrectly included its end-repair, dA-tailing and PCR steps inside the 2012 graph while describing them as separate in prose. Those imported steps are now removed. The 2019 primer sequences do not establish the 2012 library ends.

The earlier Affymetrix Prep2Seq and single-end assignments were unsupported by the inspected 2012 passage and remain unassigned. The replacement graph also removes the generic 18-base adapters and 26-base insert shell; unspecified molecular regions use role tokens with unknown lengths. A final PCR is not invented merely to reach a familiar Illumina construct.

## Remaining work

- Inspect the complete publisher mmc1 PDF, including OCR-bearing page 8, and determine the exact chemistry-relevant supplementary scope.
- Complete the mmc2 workbook content review if required for a source-disposition decision; the current check covers identity and headers only.
- Read the cited Weber 2005 MeDIP original for genomic IP/recovery conditions.
- Establish the RNA fate, DNA-end processing, selected adaptor/index identities, any amplification, physical library geometry and sequencing run mask from applicable originals.

This bounded migration changes source ownership and corrects the represented chemistry. It is neither a complete protocol review nor an accepted missing-information audit.
