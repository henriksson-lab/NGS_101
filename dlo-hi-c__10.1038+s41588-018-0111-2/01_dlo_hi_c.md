# DLO Hi-C

Lin et al.2018, *Digestion-ligation-only Hi-C is an efficient and cost-effective method for chromosome conformation capture*, DOI10.1038/s41588-018-0111-2.

## Originals inspected

🟢 Published SupplementaryInformation MOESM1, DOI/title checked: FigS1–5 captions, Table1 originalp11 and FigS2 originalp3 rendered, Table2 timeline. Mandatory scraper found15 candidates; image inspection was necessary for the FigureS2 adapters/primers. The exact publisher mainHTML was fetched, but contains abstract/figurepreview rather than fullMethods; it establishes identity without complete recipe coverage. Operational coverage below comes from the published SI. ReportingSummary MOESM2 fetched but not chemically read. No upstream drawn page used.

## Chemistry and source branch

🟢 Selected original K562/HindIII branch: EGS/formaldehyde double crosslinking, HindIII digestion/half-linker ligation with T7DNA ligase, centrifugation retaining linker-bearing chromatin pellet and excluding excess linker, in-gel T4DNA ligase proximity ligation, reverse-crosslinking/DNA purification, MmeI digestion, nativePAGE recovery of approximately80bp contact products, Illumina adapter ligation and libraryPCR. The in-situ/MseI branch has a different linker, phosphorylation operation and timings; it is not silently imported.

🟢 Table1 K562 HindIII half-linkerA: `5'-p-GATAGCTCAAGAACTCCGAC-3'` and `5'-p-AGCTGTCGGAGTTCTTGAGCTATC-3'`; B: `5'-p-CTTAGCTCAAGAACTCCGAC-3'` and `5'-p-AGCTGTCGGAGTTCTTGAGCTAAG-3'`. Table1 describes T4PNK phosphorylation then95°C5min and cooling5°C/min to25°C. Model sequences are these actualK562 rows; THP1 rows differ. Paired strands retain5′phosphate. HindIII-facing AGCT is a5′overhang. The remaining half-linker end enables subsequent linker-linker joining, not arbitrary genomic assembly.

🟢 FigureS2 presents a40bp full-linker core between two20bp genomic tags, with MmeI-derived staggered ends. Its Illumina adapters have **two-base3′NN overhangs** and5′OH ends. These must not become standardA-tailed/Y-adapter geometry merely because primer sequences resemble Illumina. FigureS2 primer2 contains a sample-index window without an assigned length/member; the exact source fixed flanks are retained, full geometry remains unassembled. The graph does not invent a barcode sequence.

## Physical fates and readout

🟡 Source selected contact class is distinct from the mixture of unligated, homotypic/heterotypic-linker, random-ligation and unreacted material. Digestion creates actual shorter DNA and remnant fragments; only the following sourcePAGE/cleanup physically separates them. Linker-class classification and computational removal of linker/mapping-tail adjustments do not alter DNA. A QC assay consumes its own aliquot, not the entire production library. Exact A/B batch allocation/mixing timing is unassembled because the full mainMethods were not obtained; the graph leaves this operation opaque rather than supplying a chosen pooling design. Cost-table SE100/SE150 comparisons do not establish a particular experimental run. Exact acquisition/primer landing and downstream complete product geometry remain unreviewed.

## Scope and linked originals

This is a private incomplete source-owned graph, not a finished schematic or accepted missing-information audit. Source-specific end/ligation/PCR/read geometry remains incomplete. No automatic inheritance from Hi-C or ChIA-PET is claimed; directly cited component predecessors require exact scope/version review. Full mainMethods, source-owned supplemental adapter/primer member assignments and directly cited original reagent details remain explicit follow-up targets; no disclosure is declared absent.
