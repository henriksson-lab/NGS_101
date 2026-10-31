# rG4-seq

Kwok et al.(2016), *rG4-seq reveals widespread formation of G-quadruplex structures in the human transcriptome*. DOI https://doi.org/10.1038/nmeth.3965 .

## Originals inspected

🟢 Published SI MOESM355 Fig1/2 workflow/copying direction and Notes2/3 limitations inspected; original published Table1 MOESM356 workbook adapter/RT/PCR rows and modification notation inspected. The full defining main Methods remain inaccessible in the standard source-tool pass (manual defining DOI row), not source-read. Mandatory scraper run; selected hits verified against original table cells. Other data/analysis supplements were fetched automatically but unreviewed. One capped pass, no absence audit.

## Molecular route and actual oligos

🟢 SI Fig1 explicitly orders RNA3′adapter ligation, RNA folding under separate Li+, K+ or K++PDS conditions, RT stalling/cDNA length differences, cDNA5′adapter ligation, PCR and NGS. Li+ is a nonstabilizing comparison; it does not physically forbid every RT stall. G-rich sequence or potassium alone does not establish a specific folded molecule, unique motif or in-vivo structure. The local graph follows every known operation while leaving unreviewed enzyme identities/workups/end geometry and run layout explicit.

Published Table1 ordering strings, with modifications:

- DNA3′adapter: `/5rApp/AGATCGGAAGAGCACACGTCTG/3SpC3/`.
- Library RTprimer: `CAGACGTGTGCTCTTCCGATCT`.
- DNA5′adapter: `/5Phos/AGATCGGAAGAGCGTCGTGTAGCTCTTCCGATCTNNNNNN/3SpC3/`.

The model computes RTprimer/3′adapter complementarity. Both terminal C3 moieties remain non-nucleotide blockers, not extra bases. The source randomN6 window is retained without assigning it UMI/sample-index semantics. Indexed PCR reverse-primer variants exist in Table1; the graph does not choose a study-specific member or invent its complete landing. Source adapter nomenclature alone does not prove the exact ligation-junction orientation. In particular, no later rG4-seq2.0 dU adapter recipe is imported.

🟡 `FoldedRNA` holds caller-specified RNA, condition and optional selected stall position; it does not predict folding/activity. The selected RT-stop example copies the RNA3′suffix into complementary firstDNA5′→3′, preserving the full RNA template. It is a coordinate/copying relation, not a fully assembled library. Unligated substrates/RNA/templates remain in reaction populations until an actual source workup can be assembled; no imaginary purification is inserted. Source SI establishes NGS but not the acquired read-role/cycle layout in this selected scope, so acquisition remains unspecified.

## Limits and lineage

🟢 Published SI Note2 states that the assay maps in-vitro rG4 structures; these may differ from cellular folding. Note3 identifies coverage-based false-negative/interpretation limitations. Validation NAI/SHALiPE gel assays are not silently appended to the NGS route. G4RP-seq later cites this primary assay as a prior strategy but uses a different capture frontend; citation/shared G4 target does not establish recipe inheritance. Full main/library workup/conditions, study index assignment and exact primer/read geometry remain deferred, not absent.

## Bounded unlinked-page adjudication (2026-10-30)

🟢 Published Table1 supplies the whole PCR forward primer and indexed reverse members I2/I4/I5/I6/I7/I12. The model now retains their known polymer bodies using canonical sequences only where they match the actual source orderings. The default reverse-primer state has an unassigned six-base source sample-index window; its placeholder is not a claim of degenerate synthesis or an assignment to an experimental library. All six printed members are checked against the original source text. Complete primer landing and whole library geometry remain unassembled. The distinct randomN6 in the DNA5adapter still has no assigned UMI or sample meaning.

The K++PDS folding branch now includes pyridostatin as a molecular contributor. Its concentration, occupancy and detailed binding state remain unreviewed; it is not a covalent RNA adduct or guaranteed G4-folding reaction. Li+, K+ and K++PDS stay separate conditions.

The original blocked_main SourceItem remains uncovered. The inspected SI Fig3 specifies four independent biological replicates for each condition; the current representative graph does not establish those producer arrays, input/fragment/end preparation, workups, QA, enzyme recipes or full run layout. Known oligo states and a library endpoint cannot substitute for the actual required Methods scope. No source coverage is generated from graph nodes, no newer rG4-seq2.0 recipe is imported, and no finished page is claimed. The source pass remains closed on the published SI and Table1 only.
