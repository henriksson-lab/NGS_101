# Prox-seq — research notes

Defining source: Vistain *et al.*, *Nature Methods* (2022), [10.1038/s41592-022-01684-z](https://doi.org/10.1038/s41592-022-01684-z) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC11289786/)). Detailed protocol: Vistain *et al.*, *Nature Protocols* (2024), [10.1038/s41596-024-01030-x](https://doi.org/10.1038/s41596-024-01030-x).

## 1. What it is

Proximity sequencing (Prox-seq) measures extracellular proteins, pairwise protein proximity and mRNA in the same single cells. Two antibody-conjugated DNA probes produce a sequenceable product only when a universal connector can bridge and ligate them.

## Sources read

- 🟢 The defining paper's probe design, reporter architecture, ligation and single-cell library sections were read in the open manuscript.
- 🟢 The later authoritative protocol confirms probe-conjugation and branch-amplification order.
- 🟢 The paper's Supplementary Methods and Supplementary Tables 7–10 (MOESM1 PDF and MOESM3 XLSX) print the connector, branch-specific library and sequencing primers, index choices and bead oligos. The downloaded `.txt` twins were checked against the named table rows.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1038/s41592-022-01684-z`.

## Reporter architecture

- 🟢 Every protein target has probe A and probe B, each an antibody conjugated to a single-stranded oligo.
- 🟢 A universal connector allows any A probe to join any B probe. The complete PLA product is 119 bases; a 20-base region spanning the joined probe ends is hybridized to the connector during ligation. The connector is a splint, not inserted sequence.
- 🟢 The ligation connector is `5′-TTTCACGACACGACACGATTTAGGTC-3′` (Supplementary Methods, and the defining article's plate-based Methods). Its exact orientation and the 20 paired bases at the probe junction require the probe-end oligos; the connector sequence alone does not specify those arms.
- 🟢 The product contains protein barcode A, protein barcode B, a UMI, a PCR-primer site and a 3′ poly(A) tail. Both probe halves are required for an amplifiable product.
- 🟢 The poly(A) tail makes the reporter compatible with the same oligo(dT) capture used for cellular mRNA.

## Construct-changing steps

1. Bind antibody–oligo probe pairs to cell-surface proteins.
2. 🟢 Hybridize connector and ligate proximal A/B probes. This is the defining covalent junction and is explicitly marked on the schematic.
3. 🟢 Process cells by Drop-seq, Smart-seq2 or 10x. The single-cell barcode joins the PLA and transcriptome measurements computationally.
4. 🟢 PCR the PLA and cDNA branches and add platform indexes.

## Final library and uncertainty

- 🟢 Each PLA read identifies both targeted proteins and its UMI; reporter counts encode protein abundance and proximity.
- 🟢 Supplementary Table 7 gives droplet-branch custom sequencing primers: Read 1 `GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`, Read 2 `AGATGCATCGTCTCGTGGGCTCGG`, and i7 `CCGAGCCCACGAGACGATGCATCT`. The paper's Methods state 20/85/8 cycles for Read 1/Read 2/i7.
- 🟢 Supplementary Table 8 gives plate-branch custom primers: Read 1 `TATTTGGCATCGTCTCGTGGGCTCGG`, i5 `CCGAGCCCACGAGACGATGCCAAATA`, and i7 `ACTCTGCGTTGATACCACTGCTTCGGC`. The Methods state 75/8/8 cycles for Read 1/i5/i7.
- 🟢 The distinct 10x branch in Supplementary Methods (MOESM1, extracted lines 317–344) uses CG000315 Rev B for cDNA and CG000317 Rev B for PLA. Its cDNA run is Read 1 28, i7 10, i5 10, Read 2 90 cycles; its PLA run is Read 1 28, i7 8, Read 2 75 cycles, with custom Read 2 and i7 primers spiked into the Illumina primers. These are not the Drop-seq 20/8/85 recipe.
- 🟢 Table 7 also prints `TSO_RNAhybrid` as `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG` and the 10x PCR primer with terminal phosphorothioates as `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C`; Table 10 prints barcoded bead oligos with `/5AzideN/`, `/5Phos/` and `/3AzideN/` modifications. These are specific to the published back ends and should not be merged into a single synthetic library.
- 🟢 The 10x branch's supplementary Methods cite CG000315 Rev B for 3′ gene expression and CG000317 Rev B for PLA/feature library construction. Both [CG000315 Rev B](https://assets.ctfassets.net/an68im79xiti/1C16trEdzy1Folq5xbOijE/7e6fb1f504e130bd561d898384da99d9/CG000315_ChromiumNextGEMSingleCell3-_GeneExpression_v3.1_DualIndex__RevB.pdf) and [CG000317 Rev B](https://assets.ctfassets.net/an68im79xiti/1W6c5yj7FAiBQdOKabc6sc/06e962ecdff018fb20b6979859cf1f5e/CG000317_ChromiumNextGEMSingleCell3-v3.1_CellSurfaceProtein_RevB.pdf) are available as original guides. Their appendices print the fixed gel-bead oligo and library shells; the Prox-seq-specific reporter remains defined by its own supplementary Tables 7–10.
- 🟡 The method supports several single-cell back ends. The current schematic represents invariant reporter roles and does not select one of the printed readout branches; its generic sequencing-primer limitation concerns that unspecialized schematic, not an absence of published branch-specific primer sequences.
