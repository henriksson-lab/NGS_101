# TOP-seq

Defining DOI: https://doi.org/10.1016/j.molcel.2016.12.012

## Sources read and scope

🟢 `j.molcel.2016.12.012_PMC5291818.xml.txt`: Associated accepted manuscript NIHMS838518; mainDesign64, Methods153–177 and KeyResources333–361. Defining original inspected before linked references. Genomic labeling/adapter/click/Pfu/TOP/PCR and IonProton source order; A1/A2 exact sequences; TO leading-T discrepancy and unreviewed LNA/PS placement.

Capped pass: at most10minutes/3originals, published supplements preferred over associated-manuscript aliases where available. No linked-reference recursion, two-reader absence audit or finished-protocol promotion.

## Source chemistry and model

🟢 Source genomic DNA is sonicated to a 150–200 bp peak. eM.SssI/Ado-6-azide tagging (0.5–1µM enzyme, 200µM cofactor, 30°C1h) targets selected unmodified CpG cytosines. Proteinase K 0.2mg/ml55°C30min precedes GeneJet cleanup. Protein digestion is not physical removal.

🟢 End filling uses the named historical DNA End Repair Kit followed by GeneJet cleanup. Klenow exo−/0.5mM dATP37°C45min and75°C15min inactivation precede another cleanup. The partially complementary A1/A2 adapters are annealed and ligated with T4 DNA ligase22°Covernight;65°C10min and Zymo cleanup follow. A1 is `5′-P-GATTGGAAGAGTGGTTCAGCAGGAATGCTGAG-3′`; A2 is `5′-ACACTCTTTCCCTACATGACACTCTTCCAATCT-3′`. Their source bodies and A1 phosphate are modeled independently; complete fork/cohesive-end/junction assembly is unreviewed.

🟢 CuBr:TBTA/60–65%DMSO clicks an alkyne-U-bearing tethered oligo onto azide-tagged genomic DNA at45°C2.5or6h. DMSO dilution to<1% precedes Zymo cleanup. Source-known EP binding/Pfu proximity-driven synthesis uses95°C1min/65°C1min/72°C1min. Physical genomic DNA is not replaced by the new copy. The source says the detailed nonhomologous priming mechanism is not investigated; the model's explicit-coordinate complementary-region relation does not predict its register, template switching or complete first strand.

🟢 EP body is `TGTTACATCTGGTAGTCAGTCTCCAATAAATATAT` with custom LNA modifications; exact positions are unreviewed. Source allocates5µl of50µlTOPreaction to15cycle Maxima Hot Start PCR with63nt barcode fusion and45nt trP1-A2 primers, both described as phosphorothioated. The independent remainder is retained, with no assigned later fate. Zymo cleanup and MagJet~250bp selection precede actual single-end IonProtonPI sequencing. Sample-index bodies occur in KeyResources but the exact modified primer/backbone/run geometry is unassembled. No unannotated barcode is drawn as an established final polymer.

### Source TO discrepancy

🟢 Within the same associated original, Methods165 prints `5′-T(alkyneU)TTATATATTTATTGGAGACTGACTACCAGATGTAACA-3′`; KeyResources333–336 prints `5′-(alkyneU)TTATATATTTATTGGAGACTGACTACCAGATGTAACA-3′`. The former has an extra leading T. Both alternatives are preserved; neither is silently repaired or selected. Alkyne-U remains a modified nucleotide, not a nucleotide-free spacer. The gated original chemistry supplement and final VOR applicability are unresolved. Any future affected schematic must show both source values and its choice/remaining uncertainty.

### Readout and linked sources

Computational trimming to70nt bacterial/165nt human reads and the ±3nt CpG call window are separate from physical sequencing cycles and exact priming geometry. Methylation estimates are fitted computational quantities, not chemical conversion yields. SourceDesign64 explicitly uses the earlier mTAG component (Kriukiene2013); exact original/version was not followed or imported. Ado-6-azide synthesis Masevicius2016 and historical kit/vendor originals remain direct pending targets. No recipe is inherited from RRBS, NOMe or matching endpoints.

## Supplement scopes and actual access barriers

The defining JATS numbers the PDF attachment1 and three XLSX attachments2/3/4, without captions mapping workbooks to supplemental tables. Main77/79 uses TableS1 for library characteristics; Methods195 identifies TableS2 as validation pyrosequencing oligos; Results/Fig5caption634 identifies DMR/gene/validation-result tablesS3–S6. Thus PDF/FigS1 is a relevant chemistry/priming target and TableS2 is relevant to validation oligos, but no unsupported claim is made that every workbook is an original chemistry table or that a numbered XLSX is TableS2. The other result-table contents remain uninspected.

## Actual access barriers and deferred originals

- https://pmc.ncbi.nlm.nih.gov/articles/instance/5291818/bin/NIHMS838518-supplement-1.pdf — PMC supplementary file; PMC 'Preparing to download' gate
- https://pmc.ncbi.nlm.nih.gov/articles/instance/5291818/bin/NIHMS838518-supplement-2.xlsx — PMC supplementary file; PMC download gate (not in the PMC Cloud bucket: licence TDM)
- https://pmc.ncbi.nlm.nih.gov/articles/instance/5291818/bin/NIHMS838518-supplement-3.xlsx — PMC supplementary file; PMC download gate (not in the PMC Cloud bucket: licence TDM)
- https://pmc.ncbi.nlm.nih.gov/articles/instance/5291818/bin/NIHMS838518-supplement-4.xlsx — PMC supplementary file; PMC download gate (not in the PMC Cloud bucket: licence TDM)

Manifest wording such as “probably redundant” is a fetcher hint, not an identity/equivalence conclusion. Exact inspected scopes above determine authority; ordinary deferred historical geometry is not an access failure.

## Validation

Protocol-owned `tools/chem.py` has a zero-argument graph and focused state/fate checks. It is an incomplete graph-only implementation, with no public page or completion marker. Source-known late operations/acquisition are retained; unassembled adapters, modifications and fates are not inferred from other families.


## Bounded publication review — 2026-10-30

Known source contributors are now explicit: eM.SssI catalyst, Ado6azide donor, and BSA/buffer components from one actual TNB buffer allocation. BSA is not an independently invented fresh dose. Proteinase K is an actual catalyst rather than a condition-only label. Klenow A-tailing receives separate actual dATP; TOP extension receives actual dNTP. The actual Maxima master-mix allocation exposes produced polymerase, nucleotide and buffer components consumed together by PCR, with no fresh duplicate component stocks or invented formulation. These functional decompositions are inventory interfaces, not physical separation reactions. Original `pending` remains unchanged and uncovered: actual SI/table scope, TO Methods extra-leading-T versus KeyResources alternative, LNA/PS positions and full library/read assembly remain unresolved. The actual 5µL TOP aliquot and remainder remain separate. No occupied whole library or acquired read geometry was invented.

No additional originals were read or acquired during implementation. This capped pass accepts no missing-information claim and does not publish the graph.
