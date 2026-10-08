# Slide-seqV2 — research notes

## 1. What it is

Slide-seqV2 captures tissue RNA on a monolayer of spatially decoded 10-µm beads. A bead-specific split-pool barcode and a molecule-specific random sequence travel with each captured transcript.

## 2. Sources read

- 🟢 Stickels et al., *Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2*, Nature Biotechnology 39, 313–319 (2021), DOI [10.1038/s41587-020-0739-1](https://doi.org/10.1038/s41587-020-0739-1); open full text [PMC8606189](https://pmc.ncbi.nlm.nih.gov/articles/PMC8606189/).
- 🟢 Online Methods lines for bead sequences, puck sequencing and library preparation.

## 3. Key oligo as printed

🟢 Chemgenes bead: `5'-TTTTTTTTCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNNT30`.

The paper defines `J` as split-pool barcode bases, `N` as mixed UMI bases and `T30` as capture oligo-dT. Fifteen J positions are synthesized; fourteen are used for spatial matching.

## 4. Construct-changing steps

1. Decode bead barcodes by in-situ ligation sequencing.
2. Hybridize tissue poly(A) RNA and reverse-transcribe with a template-switch oligo.
3. Digest tissue, recover beads, remove excess primers and synthesize the second strand.
4. PCR-amplify cDNA and use Nextera library preparation with the published TruSeq-P5 hybrid primer.

## 5. Evidence boundary

🟢 The bead sequence and TSO are printed in the paper. 🟡 The schematic normalizes the final platform arms to canonical Illumina roles so sequencing-primer placement can be computed; it does not claim that normalization as another published oligo.
