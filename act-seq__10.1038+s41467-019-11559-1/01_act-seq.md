# ACT-seq: antibody complex precedes chromatin targeting

Carter et al. (2019), [Mapping histone modifications in low cell number and single cells using antibody-guided chromatin tagmentation (ACT-seq)](https://doi.org/10.1038/s41467-019-11559-1).

## Sources and corrected identity

One reader,28October2026; three originals: published PMC6702168 JATS (Results and all preparation/library/sequencing Methods); corrected publisher Supplementary Data2 workbook, file41467_2019_11559_MOESM3_ESM.xlsx (Sheet1, all named bulk and single-cell oligo rows); [2020 Author Correction](https://doi.org/10.1038/s41467-020-18309-8), full notice. Mandatory oligo scraper run; hits/modifications verified in corrected workbook. No preprint chemistry adopted. Capped pass, not a completed review or missing-information audit.

🟢 **Published correction matters:** the original Supplementary Data2 duplicated Data1. The current publisher workbook (md5 faf90e8f9fa86272817fa1c3ccda94d7) contains the oligo list. The locally fetched PMC-named workbook (md5 7bc556f080952a6bdd587304ace63fef) instead contains mapping-statistics sheets. This is a correction and stale mirror, not conflicting oligo values. Use the corrected publisher workbook; retain both document identities. No adapter architecture may be inferred from the stale workbook.

## Scope

🟢 ProteinA is fused to the N terminus of Tn5. Adapter-loaded PA-Tnp is incubated with the selected antibody *before* exposure to permeabilized cells. Binding and washing precede Mg activation. Unlike the source-owned CUT&Tag sequence of separate chromatin antibody recruitment, ACT specifies a preassembled antibody/transposase complex. Similar final adapters do not authorize importing its recipe.

The schematic represents the published bulk method and its printed example index pair. iACT is retained as a source-described branch in these notes, not silently equated with the bulk endpoint.

## Corrected source oligos, verbatim 5′→3′

| Workbook name | Sequence |
|---|---|
| pMENTS, annealed to all complex barcodes | 5′Phos-CTGTCTCTTATACACATCT |
| Bulk-cell5′ Complex Barcode | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG |
| Bulk-cell3′ Complex Barcode | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG |
| Bulk-cell5′ Library Barcode | AATGATACGGCGACCACCGAGATCTACAC TAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG |
| Bulk-cell3′ Library Barcode | CAAGCAGAAGACGGCATACGAGAT CATGCCTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG |
| Single-cell5′ Complex Barcode#1 | ACACTCTTTCCCTACACGACGCTCTTCCGATCTCATTGGTAGATGTGTATAAGAGACAG |
| Single-cell3′ Complex Barcode#1 | GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGTAATAGGAGAGATGTGTATAAGAGACAG |
| Single-cell Universal5′ Library Barcode | AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT |
| Single-cell Unique Library Barcode#1 | CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T |

Spaces within bulk PCR rows preserve source index grouping. Workbook calls `*` phosphothiolate bonds; the molecular modification is represented as the printed terminal phosphorothioate bond, without changing the source spelling. The single-cell workbook has16 first-arm variants,12 second-arm variants and96 second-round index primers; those are distinct component families, not192 or96 additional protocols.

## Bulk reaction order

1. 🟢 PA-Tnp loads the5′/3′ complex adapters (each annealed to phosphorylated pMENTS); 10min25°C in complex formation buffer, then antibody60min25°C.
2. 🟢 Permeabilize HEK293T cells10min on ice in complex buffer (Tris,150mMNaCl,0.05%Triton,12.5%glycerol). Prepare1000-cell aliquots, bind antibody-bound PA-Tnp60min at room temperature, three washes and rinse.
3. 🟢 Activate with10mMMgCl2,30min37°C. Stop with EDTA/SDS/ProteinaseK60min55°C. Phenol/chloroform and ethanol precipitation recover DNA.
4. 🟢 Phusion PCR begins72°C5min, then15cycles98°C10s/65°C30s/72°C15s, final72°C5min. Initial extension is represented separately from PCR amplification. Gel-excise desired fragments and QIAquick purify; source assigns no universal fragment-size range in this section.
5. 🟢 HiSeq2500 paired-end sequencing. Exact cycle lengths and loaded sequencing-primer oligos are unassigned in this capped model. Printed PCR and loaded adapter sequences establish their own binding relationships; run primers are not automatically the PCR primers.

🟡 A/B, A/A and B/B end classes follow from the two loaded adapter identities. Their frequencies, actual Tn5 dimer composition and same-end physical fates are not established by the source. Draw the inferred productive A/B branch without calling all molecules dual-index libraries. Derive endpoint and index orientation from PCR primers, not colour or side labels.

## iACT branch and caveats

🟢 Results describe5000cells/well, whereas Methods “Single-cell chromatin binding” describes10000cells/well. Both use a first indexed96-well complex-binding plate, pooled cells then FACS18cells/well into a second96-well plate. Mg transposition occurs *after* redistribution. Second-round universal5′ and distinct3′ index PCR uses18cycles, followed by pooled MinElute and gel cleanup. The difference in first-well density is retained as a published discrepancy and not resolved by borrowing bulk settings. Exact barcode-subset/well pairing and single-cell read-cycle assignment are deferred.

## Ancestry and deferred sources

🟢 Article Results explicitly calls iACT barcode multiplexing similar to references14/15 (Rosenberg2018 SPLiT-seq; Cusanovich2015 combinatorial indexing). This supports conceptual combinatorial-index ancestry, not inheritance of their adapter sequences or target chemistry. Article Methods points to the Illumina Adapter Sequences document only for design background. Exact oligos are independently supplied by corrected Data2. Linked ProtocolExchange DOI10.1038/protex.2019.026 and source-linked Addgene121137 are deferred within the three-original cap; no chemistry-completeness claim depends on unread versions. Other single-cell-family members, detailed run masks and second independent reader remain open; no absence claim accepted.
