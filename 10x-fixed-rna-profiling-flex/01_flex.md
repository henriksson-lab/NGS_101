# 10x Fixed RNA Profiling (Flex) — research notes

## 1. What it is

Flex measures RNA in fixed cells or nuclei through pairs of gene-specific DNA probes. Adjacent probe halves hybridize to one RNA target and are ligated; the reporter product is then cell-barcoded in GEMs.

## 2. Sources read

- 🟢 10x Genomics, *Chromium Fixed RNA Profiling Reagent Kits for Singleplexed Samples with Feature Barcode technology for Protein using Barcode Oligo Capture*, CG000674 Rev C: [support page](https://www.10xgenomics.com/support/cn/flex-gene-expression/documentation/steps/library-prep/chromium-fixed-rna-profiling-reagent-kits-for-singleplexed-samples-with-feature-barcode-technology-for-protein-using-barcode-oligo-capture).
- 🟢 Earlier CG000477 Rev E library and sequencing diagrams were checked for the same fixed-RNA probe-pair architecture.
- 🟢 [CG000674 Rev C direct PDF](https://cdn.10xgenomics.com/image/upload/v1729115031/support-documents/CG000674_ChromiumFixedRNAProfiling_SingleplexedSamples_ProteinExpression_BarcodeOligoCapture_UserGuide_Rev_C.pdf), Appendix p. 85, **Oligonucleotide Sequences**, prints the Flex gel-bead primer and both completed gene-expression and protein-expression library layouts. Its sequencing section specifies TruSeq Read 1 / Small RNA Read 2 for gene expression, and TruSeq Read 1 / Nextera Read 2 for protein expression.

## 3. Published oligo structure

🟢 The gel-bead primer is `5′-CTACACGACGCTCTTCCGATCT-N16-N12-TTGCTAGGACCG-3′`. Here `N16` is the GEM barcode and `N12` the UMI. The terminal `TTGCTAGGACCG` is printed, and should not be treated as an unknown transfer-oligo sequence.

🟢 The gene-expression library is printed with a 10-base i5 index, TruSeq Read 1 arm, GEM barcode, UMI, the bead-primer terminus, `BC8-NN-TACGTGCTAACCGCGT`, a ligated-probe insert, the Small RNA Read-2-side constant `TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC`, a 10-base i7 index and P7. `BC8` is a probe barcode window, not an extra GEM cell barcode. The full duplex appears on p. 85 of the guide; these exact constants are now surfaced on the schematic as key oligos.

🟢 The [10x Probe Set Reference CSV description](https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-probe-set-files) says downloadable Cell Ranger probe-set files contain the **combined LHS/RHS target sequence** for each probe, along with probe-barcode sequence files. 🔴 This guide does not print the complete ordered LHS/RHS oligos, including any synthesis-specific modifications. The page's compact final duplex still uses normalized platform arms and must be rebuilt against the source-exact Small RNA Read-2 end before it can be considered base-complete.

## 4. Construct-changing steps

1. Fix and permeabilize the sample.
2. Hybridize target-specific left and right probe halves.
3. Ligate adjacent probe pairs.
4. Partition cells or nuclei into GEMs and transfer cell barcode and UMI to reporters.
5. Pre-amplify and add sample indexes.

## 5. Evidence boundary

🟢 The current manual establishes adjacent-probe ligation, read roles, the gel-bead sequence and the complete library shell. The downloadable probe-set reference supplies target-region sequences, but the complete ordered probe halves have not been reconstructed here. 🟡 The current schematic preserves assay roles but still normalizes the source-specific Read-2 arm; its caveat calls that out explicitly.
