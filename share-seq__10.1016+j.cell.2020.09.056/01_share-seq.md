# SHARE-seq — research notes

## 1. What it is

SHARE-seq jointly profiles accessible chromatin and RNA in fixed single cells using the same three-round split-pool barcode combination for both molecular classes.

## 2. Sources read

- 🟢 Ma et al., *Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin*, Cell 183, 1103–1116.e20 (2020), DOI [10.1016/j.cell.2020.09.056](https://doi.org/10.1016/j.cell.2020.09.056), open manuscript [PMC7669735](https://pmc.ncbi.nlm.nih.gov/articles/PMC7669735/).
- 🟢 The paper's Methods and Figure S1 establish the reaction order, three 96-well hybridization/ligation rounds, biotin separation and two final libraries.
- 🟢 Current author protocol: [protocols.io SHARE-seq v2.2](https://www.protocols.io/view/share-seq-protocol-v2-2-g2y6byfzf.pdf).
- 🟢 Supplementary Table S1 gives the RT primer, template-switching oligo, RNA PCR primer, all three barcode-oligo constant arms, all three linkers and the blocking oligos. The upstream scg_lib_structs reconstruction was checked base-for-base against that table before these sequences were put in the model.

## 3. Construct-changing steps

1. Tn5-tagment accessible chromatin in fixed permeabilized cells. The ligatable Tn5 product and the RT primer expose the same 15-nt `GTCTCGTGGGCTCGG` arm.
2. Reverse-transcribe RNA with the 5′-phosphorylated primer carrying a Nextera-derived handle, 10-nt UMI, internal biotin-dT and oligo-dT/VN tract.
3. Assemble the barcode chain through three 96-well hybridization rounds. After each round, a complementary blocking oligo sequesters free adapters before cells are pooled and redistributed. Each linker is exactly the reverse complement of the incoming molecule arm followed by the reverse complement of that round's barcode-oligo 3′ arm; the model constructs and validates this geometry rather than storing it as a caption.
4. After the third hybridization/blocking round, add T4 DNA ligase once to seal all three backbone nicks.
5. Reverse crosslinks and capture biotinylated cDNA on streptavidin beads.
6. Construct separate RNA and ATAC sequencing libraries whose shared barcode triplet identifies the same cell.

## 4. Evidence boundary

🟢 The defining paper and protocol specify the constant arms, linkers, modifications and 96 barcode values per round. The public page draws the exact constant geometry and modifications. 🟡 It uses structured 8-nt barcode placeholders rather than embedding all 288 well-specific values in every reaction state.
