# SHARE-seq — research notes

## 1. What it is

SHARE-seq jointly profiles accessible chromatin and RNA in fixed single cells using the same three-round split-pool barcode combination for both molecular classes.

## 2. Sources read

- 🟢 Ma et al., *Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin*, Cell 183, 1103–1116.e20 (2020), DOI [10.1016/j.cell.2020.09.056](https://doi.org/10.1016/j.cell.2020.09.056), open manuscript [PMC7669735](https://pmc.ncbi.nlm.nih.gov/articles/PMC7669735/).
- 🟢 The paper's Methods and Figure S1 establish the reaction order, three 96-well hybridization/ligation rounds, biotin separation and two final libraries.
- 🟢 Current author protocol: [protocols.io SHARE-seq v2.2](https://www.protocols.io/view/share-seq-protocol-v2-2-g2y6byfzf.pdf).

## 3. Construct-changing steps

1. Tn5-tagment accessible chromatin in fixed permeabilized cells.
2. Reverse-transcribe RNA with a biotinylated oligo-dT primer carrying a UMI.
3. Repeat hybridization-guided barcode ligation through three 96-well plates.
4. Reverse crosslinks and capture biotinylated cDNA on streptavidin beads.
5. Construct separate RNA and ATAC sequencing libraries whose shared barcode triplet identifies the same cell.

## 4. Evidence boundary

🟢 The defining paper and protocol specify exact plate oligos in supplementary material. 🟡 This page retains the verified oligo roles and junction order but uses barcode placeholders rather than embedding 288 well-specific sequences.
