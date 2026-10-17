# SEC-seq — source notes

Defining paper: Cheng et al., *Nature Communications* (2023), [doi:10.1038/s41467-023-39367-8](https://doi.org/10.1038/s41467-023-39367-8).

- 🟢 The PMC full text and supplementary methods were read.
- 🟢 Antibody-coated 55-µm nanovials retain cells and capture secreted IgG; PE anti-IgG followed by a barcoded anti-PE antibody reports the captured secretion in the main human B-cell experiment.
- 🟢 Cell-loaded nanovials enter a standard 10x Chromium 3′ v3.1 workflow. Transcript and antibody-derived-tag libraries are pooled 80:20 and sequenced 28-bp Read 1, 10-bp i7, 10-bp i5 and 90-bp Read 2.
- 🟢 The paper also describes a separate Raji/hybridoma compatibility experiment using TotalSeq-C streptavidin (BioLegend 405271) and anti-IgG1 (406636) with the Chromium 5′ v2 kit. Those reagents and that kit should not be conflated with the main 3′ v3.1 anti-PE readout.
- 🟢 Supplementary Table 2 (MOESM1; rendered PDF table, `SEC-seq` rows) names the main-experiment reagents: PE Mouse Anti-Human IgG (BD 555787, clone G18-145) and TotalSeq™-B0911 anti-phycoerythrin (PE) Antibody (BioLegend 408113, clone PE001). The nanovial capture coat is streptavidin plus biotinylated goat anti-human IgG Fc (Southern Biotech 2014-08) and biotinylated anti-CD27 or anti-CD45; the same anti-IgG is the 25 µg/mL blocking antibody added to cells before loading (Methods "coating", "cell loading").
- 🟢 Front-end order (Methods "cell loading", "incubation", "SEC-seq human ASCs"): load cells 1 h on ice, strain away unbound cells, incubate 30 min at 37 °C to accumulate secretion, stain PE anti-IgG then barcoded anti-PE, sort SYTOX Green-negative cell-loaded nanovials on a Nanocellect WOLF (collected independent of secretion amount, Supplementary Fig. 8a), load 2,500–10,000 per Chip G lane.
- 🟡 The 15-base TotalSeq-B0911 barcode is not printed in the inspected article or supplement; the BioLegend product documentation is a deferred original, so the feature barcode stays symbolic. The vendor guide below supplies the capture and library families; it does not establish the selected index members.
- 🟡 Raw reads are listed under GEO `GSE229042` and BioProject `PRJNA953084`. A feature-barcode consensus from those reads could help identify a used barcode after comparison with a vendor reference, but would not establish the complete antibody-conjugated oligo or its chemical modifications.

## Source-specific library correction

🟢 The human ASC Methods explicitly cite CG000317, not a generic 10x family.
The official CG000317 Rev E Appendix (pp. 71–72) prints the dual-index 3′ v3.1
GEX and cell-surface-protein feature products. The gene-expression library uses
TruSeq Read 1. The feature library uses Nextera Read 1 followed by the cell
barcode, UMI, Capture Sequence 1, an N9 flank, N15 feature barcode, N10 flank
and the TruSeq Read 2 arm. Both products have ten-base i5 and i7 indexes.
The source-known canonical arms and complete guide product order are checked
against the downloaded original in the protocol self-test.

The former shared endpoint in this repository put Nextera Read 1 on the GEX
branch, gave both branches an eight-base single index, and placed a poly(T)
tract in the secretion-tag branch. Those were repository modelling errors;
they are corrected from the defining Methods and cited vendor-guide family.
The published paper and vendor guide are not in conflict on these points.

The selected barcode/index bases and exact guide revision used in the
experiment remain unassigned. Rev E provides the displayed family geometry.
The separate Raji/hybridoma validation uses 5′ v2 and TotalSeq-C; its chemistry
has not been imported into the main human ASC model. Capture of two analytes
shares cell identity, while transcript and tag molecules have separate UMIs.

The GEX graph now includes cDNA amplification, fragmentation, end repair/A-tailing,
and Read 2 adapter ligation before dual-index PCR. The ligated junction is marked.
The exact ligation-ready PN-2000094 end geometry is a remaining model gap and
produces a validation warning; it is not claimed absent from published sources.

## Apparent cross-guide TT prefix inconsistency

🟢 The original CG000317 Rev E p.72 prints the TT forward primer and final GEX
product with `AATGATACGGCGACCACCGAGATCTACAC` (29 bases) before N10 i5.
The original Multiome CG000338 Rev A p.80 prints the GEX product with
`AATGATACGGCGACCACCGAGATCT` (25 bases) before N10. Both the top and
complementary bottom strands show the difference in the rendered originals;
it is not a text-extraction artifact. Both guides list Dual Index Kit TT Set A
PN-1000215 and plate PN-3000431 (CG000317 kit list and CG000338 p.12).

🟡 This is an apparent published cross-guide inconsistency. Different named
reagent variants cannot be assumed from these sources. The SEC model retains
the CG000317 Rev E 29-base prefix, supported by its printed forward primer
and PCR product; TEA's distinct declaration follows its cited Multiome Rev A.
The historical physical prefix used in SEC-seq remains unverified. The public
page displays this disagreement and source choice.

## Bounded pass 2026-10-17 (round breadth-first-18)

Re-read the published JATS Methods (coating, loading, incubation, staining,
mixed-species and human ASC paragraphs) and MOESM1 Supplementary Table 2 and
Fig. 8 caption. The graph now separates loading/strainer wash, 37 °C secretion
accumulation, two-antibody staining and the SYTOX Green live-nanovial sort, and
names TotalSeq-B0911 anti-PE. Still deferred: BioLegend 408113 barcode sheet,
the author analysis repository, the 5′ v2 guide for the validation branch, and
the 80:20 pool (described in prose, not drawn as a graph pool).
