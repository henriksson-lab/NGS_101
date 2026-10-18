# ChIPmentation — research notes

## 1. What it is

ChIPmentation replaces conventional post-ChIP end repair and adapter ligation with one
Tn5 tagmentation reaction performed on bead-bound immunoprecipitated chromatin.

## 2. Sources read

- 🟢 Schmidl et al. 2015, DOI
  [10.1038/nmeth.3542](https://doi.org/10.1038/nmeth.3542), defining paper and Methods.
- 🟢 The Online Methods identify Tagment DNA Enzyme from the Nextera DNA Sample Prep Kit and state that library amplification used custom Nextera primers **as described for ATAC-seq** (Buenrostro et al. 2013, Supplementary Table 1). This is an explicit source link for Ad1_noMX and the Ad2.x index-primer family, rather than an unpublished ChIPmentation-specific oligo set. The ChIPmentation supplements describe three ChIP variants and provide figure and sequencing summaries, but no separate oligo ordering table.

🟢 The original NIHMS64423 Supplementary Table 2, *Sequencing Summary*, lists
ChIPmentation runs in both 50-cycle single-read and 100-cycle paired-end
configurations (44 and 9 rows respectively). The native endpoint now uses unassigned primer-end roles; these source run variants are retained without assigning one universal mask.

## 3. Molecular order

🟢 Crosslinked chromatin is sonicated and immunoprecipitated. While the complexes remain
on beads, Nextera Tn5 fragments and adapter-tags the selected DNA. Proteinase treatment
releases DNA and reverses crosslinks; indexed PCR completes the flow-cell arms.

## 4. Relationship to neighboring protocols

ChIPmentation selects target-associated chromatin before free Tn5 acts. CUT&Tag instead
tethers Tn5 to the antibody inside permeabilized cells. That difference in targeting
geometry is why these remain separate protocol pages despite related adapter families; no shared native endpoint assembly is assumed.

## Bounded curation round 2 (2026-10-16)

Inspected: Original PMC4589892 JATS58–76.

Added bead-bound37°C/10min tagmentation and qPCR-Cq-dependent enrichment labels. Distinguished purified-DNA55°C/5min variant; candidate Tn5 products have unresolved recovery rather than presumed discard.

Deferred: Original Supplementary Note and sequencing Table2; exact custom primer/index assignment; post-PCR size-selection graph states. This capped pass does not establish a complete source review or publication-wide absence.


## Bounded native migration — 2026-10-18 batch09

Inspected 3 originals: Published PMC4589892 JATS Methods57–76 ChIP/beadbound and purifiedDNA tagmentation, qPCR/enrichment/recovery; publisher MOESM631 SupplementaryTable2 SequencingSummary heading and ChIPmentationrows25–26 SR50/PE100; publisher MOESM629 SupplementaryFigure1 opening/workflowlegend. Mandatory bundle scraper completed0hits, not treated as absence.

Removed chromatin_epigenetics endpoint/scene defaults. Native source positive beadbound37°C10min route, washes, chosenChIP finalwash/elution/reversecrosslink/purification, analytical1µlqPCRaliquot, KAPA preactivation/actual72°Cnicktranslation and qPCR-selected enrichment are explicit. Unknown variant-specific ChIP reagents, native loaded adapter/primer strands/gap geometry remain unassigned. Real taggedDNA input feedsPCR, independent residualorigins preserved. Explicit0.7×retain-supernatant long-fragment removal followed by2×recovery; source differingSR/PE examples retained without canonical indexwidth/acquisition.

Deferred: Detailed chosen ChIP supplementary protocol variant and linked ATAC customprimer originals; Native Nextera loaded strand/gap and first-cycle nicktranslation geometry; Exact study primer/indexmember/width and full run-specific acquisition masks; remainingTable2rows not re-curated. This capped pass does not accept source-wide absence or establish complete experimental validation. Earlier dated/provisional interpretations above are historical; this entry describes the current source scope.
