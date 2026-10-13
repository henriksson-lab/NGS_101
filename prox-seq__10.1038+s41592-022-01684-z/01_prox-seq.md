# Prox-seq — research notes

Defining source: Vistain *et al.*, *Nature Methods* (2022), [10.1038/s41592-022-01684-z](https://doi.org/10.1038/s41592-022-01684-z) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC11289786/)). Detailed protocol: Vistain *et al.*, *Nature Protocols* (2024), [10.1038/s41596-024-01030-x](https://doi.org/10.1038/s41596-024-01030-x).

## 1. What it is

Proximity sequencing (Prox-seq) measures extracellular proteins, pairwise protein proximity and mRNA in the same single cells. Two antibody-conjugated DNA probes produce a sequenceable product only when a universal connector can bridge and ligate them.

## Sources read

- 🟢 The defining paper's probe design, reporter architecture, ligation and single-cell library sections were read in the open manuscript.
- 🟢 The later authoritative protocol confirms probe-conjugation and branch-amplification order.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1038/s41592-022-01684-z`.

## Reporter architecture

- 🟢 Every protein target has probe A and probe B, each an antibody conjugated to a single-stranded oligo.
- 🟢 A universal connector allows any A probe to join any B probe. The complete PLA product is 119 bases; a 20-base region spanning the joined probe ends is hybridized to the connector during ligation. The connector is a splint, not inserted sequence.
- 🟢 The product contains protein barcode A, protein barcode B, a UMI, a PCR-primer site and a 3′ poly(A) tail. Both probe halves are required for an amplifiable product.
- 🟢 The poly(A) tail makes the reporter compatible with the same oligo(dT) capture used for cellular mRNA.

## Construct-changing steps

1. Bind antibody–oligo probe pairs to cell-surface proteins.
2. 🟢 Hybridize connector and ligate proximal A/B probes. This is the defining covalent junction and is explicitly marked on the schematic.
3. 🟢 Process cells by Drop-seq, Smart-seq2 or 10x. The single-cell barcode joins the PLA and transcriptome measurements computationally.
4. 🟢 PCR the PLA and cDNA branches and add platform indexes.

## Final library and uncertainty

- 🟢 Each PLA read identifies both targeted proteins and its UMI; reporter counts encode protein abundance and proximity.
- 🟡 The method deliberately supports several single-cell back ends. There is no one source-supported universal cell-barcode length, read cycle count or sequencing-primer geometry, so the page models invariant reporter regions and declares the platform-dependent endpoint.
