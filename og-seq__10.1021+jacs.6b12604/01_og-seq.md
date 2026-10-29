# OG-seq

Defining DOI: https://doi.org/10.1021/jacs.6b12604

Ding et al. (2017), Sequencing the Mouse Genome for the Oxidatively Modified Base 8-Oxo-7,8-dihydroguanine by OG-Seq.

## Source scope

🟢 Associated accepted author manuscript PMC5440228, rationale and chemical enrichment (text24–32), genomic library/readout32–36. Manifest identifies author manuscript=True; full publisher comparison and chemistry-relevant Supporting Information remain unreviewed. The Supporting Information DOCX request returned the PMC download gate: https://pmc.ncbi.nlm.nih.gov/articles/instance/5440228/bin/NIHMS860150-supplement-Supporting_Information.docx . No oligo absence is inferred from that failure or the zero scraper hits in the main text.

## Selected genomic enrichment route

1. 🟢 Genomic DNA from mouse embryonic fibroblasts is sonicated near150bp before chemical labeling. WT and Ogg1-knockout are biological input branches, not an in-vitro glycosylase pretreatment invented here.
2. 🟢 Mild oxidation with K2IrBr6 in the presence of an amine-terminated biotin-PEG reagent (BTN) traps an oxidized-OG intermediate as a stable amine conjugate. The source's test conditions are5mMK2IrBr6/20mMBTN,100mMNaPi pH8,75°C30min; source says the reaction was applied to fragmented genomic DNA. Exact genomic reaction/workup parameters require Supporting Information. OG is a reactive base class, not an assigned genomic coordinate; no conversion yield or complete selectivity is predicted. The adduct is not drawn as unchanged OG and no detailed constitutional structure is invented.
3. 🟢 Recover biotinylated strands using streptavidin magnetic beads. The selected DNA and its complementary strand are distinct molecules associated by pairing, not one covalent polymer.
4. 🟢150mMNaOH20°C10min releases the complementary strands while the labeled adduct strand remains bead-associated. These complementary strands, not the modified biotin-bearing strand, provide the NGS library substrate. Amine adducts are source-described as stable to release conditions. Capture/alkali fractions remain explicit; complete purity is not assumed.
5. 🟢 A single-stranded DNA adapter kit prepares approximately10ng released strands for HiSeq125-cycle paired-end v4 sequencing. Kit identity/revision, adapters, amplification and full run-primer/index geometry are unassembled. The untreated fragmented-WT input control is source-positive but not merged into this selected enrichment branch.
6. 🟢 Enrichment is interpreted at approximately150bp region scale against input, not at exact OG base resolution. Read letters come from the complementary DNA and do not chemically encode a retained OG adduct. Alignment, repeat processing and peak calling occur computationally.

## Boundary and lineage

The private partial model preserves oxidized-OG/amine-biotin tagging, association and complementary-strand release. It cannot validate complete adduct chemistry, capture yield or whole-library sequences. It has an uncovered Supporting Information source item and no finished marker.

The main text explicitly builds on the Greenberg chemical-tagging approach (ref17, Bajacan et al.2014 DOI10.1021/tx500120p), with commercial BTN and a selected K2IrBr6 oxidant here. This is component-scoped chemical precedent, not an imported NGS preparation; that original remains deferred. Antibody control enrichment is a separate comparison, not a parent recipe. Full Supporting Information/VOR, historical ssDNA-kit documents, remaining controls and independent review remain open. No confirmed missing-information claim is accepted.
