# 10x Fixed RNA Profiling (Flex) — research notes

## 1. What it is

Flex measures RNA in fixed cells or nuclei through pairs of gene-specific DNA probes. Adjacent probe halves hybridize to one RNA target and are ligated; the reporter product is then cell-barcoded in GEMs.

## 2. Sources read

- 🟢 10x Genomics, *Chromium Fixed RNA Profiling Reagent Kits for Singleplexed Samples with Feature Barcode technology for Protein using Barcode Oligo Capture*, CG000674 Rev C: [support page](https://www.10xgenomics.com/support/cn/flex-gene-expression/documentation/steps/library-prep/chromium-fixed-rna-profiling-reagent-kits-for-singleplexed-samples-with-feature-barcode-technology-for-protein-using-barcode-oligo-capture).
- 🟢 Earlier CG000477 Rev E library and sequencing diagrams were checked for the same fixed-RNA probe-pair architecture.

## 3. Construct-changing steps

1. Fix and permeabilize the sample.
2. Hybridize target-specific left and right probe halves.
3. Ligate adjacent probe pairs.
4. Partition cells or nuclei into GEMs and transfer cell barcode and UMI to reporters.
5. Pre-amplify and add sample indexes.

## 4. Evidence boundary

🟢 The manuals establish adjacent-probe ligation and the read roles. 🔴 The exact commercial probe and bead-oligo sequences are not published, so the schematic does not invent them.
