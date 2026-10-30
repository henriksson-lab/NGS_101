# Constellation-Seq: selected published Drop-seq branch

Defining paper: Vallejo et al., *Resolving cellular systems by ultra-sensitive and economical single-cell transcriptome filtering* (2021), DOI10.1016/j.isci.2021.102147.

🟢 The distinct preparative change is targeted single-primer linear amplification after cDNA synthesis. This is a method entry, not an additional target-panel count. The selected Drop-seq branch omits TSO during reverse transcription; the10X branch begins from separately prepared amplified cDNA and is not merged into this graph.

🟢 Transparent Methods761–813 adds the last14bases of the SMART primer, `TATCAACGCAGAGT`, to each target-specific primer. The whole published TableS2 rows are declared by `chem.SOURCE_PRIMERS`; representative LDHB ordering is `TATCAACGCAGAGTGGTGTGAATGTGGCAGGTGTT`. Actual table orderings are retained, including repeated rows; general design-length guidance is not used to truncate them.

🟢 Selected operations: prepared no-TSO RT bead stock → aliquot2000beads(~100STAMPs) →50ul KAPA reaction with hybrid-primer pool →20rounds60°C linear extension → standard Drop-seq PCR/library continuation →two1:0.6AMPure purifications →analytical Bioanalyser →tagmentation/library preparation →NextSeq500 paired20×50bp. The source's10uM pool wording does not establish a different stock/final dosing convention. One actual ReadyMix addition supplies co-present enzyme, structural nucleotide and buffer components.

🟡 Prior bead/cDNA construction, target binding sites, full products, historical PCR/tagmentation stocks and read landings are unassembled. A single known common handle does not authorize inheriting the whole Drop-seq construct. The aliquot remainder and analytical material have separate fates. No source absence, complete chemistry validation or public schematic readiness is claimed.

## Selected sources and review scope

- j.isci.2021.102147_PMC7900351.xml: Published main Results91/97/111 and supplemental-role list283–289; selected targeted-linear-amplification boundary.
- j.isci.2021.102147_supp_mmc1.pdf: Published transparentMethods761–813: primer design; no-TSO Drop-seq branch and standard comparator separated; 10X branch824–831 read only as alternate scope.
- j.isci.2021.102147_supp_mmc3.xlsx: Published TableS2 CDS52 actual whole source rows4–55; original workbook cells and extracted orderings checked.

This capped pass is not a complete source, linked-information, chemistry or absence audit.
