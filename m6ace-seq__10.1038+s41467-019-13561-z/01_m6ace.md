# m6ACE-seq

Defining original: https://doi.org/10.1038/s41467-019-13561-z

Koh, Goh and Goh (2019), *Atlas of quantitative single-base-resolution N6-methyl-adenine methylomes*.

🟢 Published original `s41467-019-13561-z_PMC6904561.xml.txt`, Results39–44, discussion130 and Methods168–170; Supplementary Data1 original workbook `MOESM4_ESM.xlsx`, named methylated-spike and library-oligo rows; original `MOESM3_ESM.docx` opening source-file index. The JATS identifies a non-manuscript published article. Other supplements/result workbooks were not chemically reviewed. Full-bundle oligo scraper run after retrieval; relevant workbook hits checked at their source rows.

The distinguishing operation is antibody/UV protection from 5′→3′ XRN1. The source single-m6A spike places m6A at position21 and reports retained RNA starting exactly at m6A: upstream nucleotides are digested, the protected nucleotide remains. This is a selected protected-substrate event; the model does not predict crosslink efficiency, off-target protection or digestion completeness.

## Oligos as printed

🟢 Supplementary Data1: 3′ DNA adapter `/5rApp/TGGAATTCTCGGGTGCCAAGG/3ddC/`; RT primer `GCCTTGGCACCCGAGAATTCCA`; 5′ RNA adapter `GUUCAGAGUUCUACAGUCCGACGAUCNNNNNNNB`. The final B restricts the eight-base UMI; it is not N8. Its identifier role is structured locally.

The terminal ddC is an additional cytidine nucleotide with no 3′OH. Its C participates in the template: the RT primer is computed as the reverse complement of the adapter body plus terminal C, including the G that pairs with ddC. The original adapter retains its blocked end when used as a template. The local adapter object preserves nucleotide identity, sugar and blocked-end chemistry; it is not reduced to an annotation discarded by sequence arithmetic.

🟢 Methylated spike: `UCUCCUUAGUAGCUCCUAAGm6AUCGUCGAGUUACACGACGACAUUGUUCCGACUGACG`. It is added already 3′-ligated after the bulk RNA adapter cleanup. The inspected workbook does not specify its starting 5′end; the model does not silently apply the earlier bulk PNK workup to this separately added spike. Synthetic event tests explicitly supply the required 5′phosphate condition.

## Source reaction order and local boundary

🟢 Methods170: poly(A) selection; fragment120–150nt at70°C7.5min; T4PNK37°C30min then ATP/30min; 3′adapter ligation by truncatedT4RNA ligase2 in ATP-free buffer25°C2h; AMPure; add preligated spike; anti-m6A overnight4°C; UV254nm six0.15J/cm² exposures; set aside1% input after UV; capture remaining sample with proteinA beads; wash; XRN137°C1h; wash; proteinaseK50°C1.5h; ethanol; 5′RNA adapter ligation16°C16h; SSIII50°C1h; PhusionHF14 cycles with named TruSeq PCR primers; NextSeqPE75.

The graph covers the bead-bound, crosslinked bulk-RNA prefix through XRN1, physical washing, proteinase elution and ethanol precipitation. Released nucleotides remain molecular material until the washes. Proteinase treatment does not assert perfect removal of crosslink remnants. The input split is recorded above but not assembled as a graph branch yet. RT RNA degradation is not borrowed from the separate RT-qPCR validation assay.

## Readout, open work and ancestry

🟢 The protected RNA 5′end provides a read-start signal, distinct from miCLIP mutation/truncation mapping. Quantitative recovery requires calibration/input normalization, not a universal site-detection probability. The source also maps m6Am; intact-cap/XRN1 geometry is outside this ordinary5′p event model.

🟡 Local constructor paths are conditional mechanistic relations, not sequenced library reconstruction. Both adapter ligations, RNA/cDNA state through RT, exact source-version PCR/index chemistry and acquired read spans remain unassembled. Notes status; no finished schematic or absence claim.

Ancestry suggestion: source Results39 explicitly cites the earlier DNA6mA antibody/crosslink/exonuclease strategy (ref31, DOI10.1093/nar/gky1104). Scope is the protection strategy; transferring DNA chemistry or adapters is not justified. The predecessor has not been independently read in this capped pass. Similar small-RNA adapter motifs do not establish derivation from any catalogue protocol.
