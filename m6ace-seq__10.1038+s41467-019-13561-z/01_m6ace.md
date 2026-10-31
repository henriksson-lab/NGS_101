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

The graph now covers actual shared RNA preparation, the UV-aliquot array and both input and enriched DNA-library producer paths. Released nucleotides remain molecular material until the washes. Proteinase treatment does not assert perfect removal of crosslink remnants. The actual post-UV 1% input split is represented as its own retained library branch. RT RNA degradation is not borrowed from the separate RT-qPCR validation assay.

## Readout, open work and ancestry

🟢 The protected RNA 5′end provides a read-start signal, distinct from miCLIP mutation/truncation mapping. Quantitative recovery requires calibration/input normalization, not a universal site-detection probability. The source also maps m6Am; intact-cap/XRN1 geometry is outside this ordinary5′p event model.

🟡 Local constructor paths are conditional mechanistic relations, not sequenced library reconstruction. Both physical adapter ligations, RT and source PE75 are represented; exact occupied RNA/cDNA molecules, historical PCR/index primer bodies and full read landing remain unassembled. The public page is a source-qualified partial published instance; no absence claim is accepted.

Ancestry suggestion: source Results39 explicitly cites the earlier DNA6mA antibody/crosslink/exonuclease strategy (ref31, DOI10.1093/nar/gky1104). Scope is the protection strategy; transferring DNA chemistry or adapters is not justified. The predecessor has not been independently read in this capped pass. Similar small-RNA adapter motifs do not establish derivation from any catalogue protocol.

## Bounded complete-preparation model

The same closed published Methods 170 and Supplementary Data 1 support the
actual upstream and downstream material stages. No DNA predecessor recipe or
modern TruSeq constructor was imported. Bulk RNA end treatment precedes its
activated DNA-adapter ligation and cleanup. The methylated spike is added
already 3′-ligated afterwards; its starting 5′ end remains unassigned, so the
earlier bulk PNK treatment is not retroactively applied to that spike.

Actual 50 µl pre-UV aliquots are represented as a first aliquot and an additional
aliquot collection, with exact count and final-partial-portion geometry
unassigned. They receive the source six UV exposures as distinct physical
classes and are recombined after treatment. The subsequent physical split
retains 1% input and 99% material for capture/XRN1. The input branch does not
inherit enriched-branch Protein A capture or exonuclease treatment.

Actual BSA-blocked Protein A preparation, molecular BSA/SDS/PEG contributions,
wash fractions and soluble digestion products retain physical provenance. Both
branches continue through proteinase recovery, ethanol precipitation, RNA5
ligation, SSIII and source RT stop. RT supply components and Phusion master-mix
components come from their own actual supplies without duplicate doses or
invented stock packaging. Historical whole TruSeq primer bodies remain
unassembled. Fourteen-cycle PCR, AMPure purification, PE75 loading and stock
remainders complete the source-qualified DNA producer paths. No extra library
QC aliquot or pool between controls is invented.

A narrow source-format exception preserves the physical terminal ddC nucleotide
that the generic source parser omits. Whole modified source-order checks and
the RT primer computed across that additional C verify its identity. The
original template's blocked end is not propagated into PCR copies. N7B is the
source-constrained eight-base UMI, not N8.

Normal publication checks pass. Whole occupied transcript/adduct molecules,
crosslink remnants, protection efficiency and read-primer landing remain
unassembled. The selected protection event requires an explicitly eligible
5′-phosphate/crosslinked substrate; it does not make all unknown spike molecules
eligible or predict global recovery. Source m6Am/cap geometry remains outside
that selected ordinary 5′-phosphate event.
