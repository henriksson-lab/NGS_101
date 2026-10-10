# Pore-C — research notes

## Defining source

Deshpande AS, Ulahannan N, Pendleton M, et al. *Identifying synergistic high-order 3D chromatin conformations from genome-scale nanopore concatemer sequencing.* Nature Biotechnology 2022. [doi:10.1038/s41587-022-01289-z](https://doi.org/10.1038/s41587-022-01289-z)

Operational source read: Oxford Nanopore Technologies, *Restriction Enzyme Pore-C protocol* and accompanying information sheet. This is the openly available wet-lab procedure underlying the paper's method class.

## Evidence record

- 🟢 Formaldehyde fixes spatially proximal chromatin within intact nuclei.
- 🟢 In-situ restriction digestion exposes cohesive ends. The source workflow supports NlaIII and DpnII and the defining Pore-C analysis uses restriction-fragment boundaries; this page draws the common NlaIII/CATG implementation.
- 🟢 T4 DNA ligase diffuses into the crosslinked nuclear scaffold and joins nearby restriction monomers into chimeric Pore-C polymers. The ONT procedure uses 16 °C for 6 h and warns that prolonged ligation can increase trans contacts.
- 🟢 Proteinase K digestion reverses crosslinks and releases the covalent dsDNA polymers. Phenol/chloroform extraction, ethanol precipitation and long-fragment selection preserve those polymers.
- 🟢 DNA ends are repaired and dA-tailed before ONT ligation-library preparation. The resulting long read contains a variable-order concatemer of restriction fragments and can therefore encode a multiway contact directly.
- 🟢 The method is PCR-free when DNA yield permits, retaining native methylation.
- 🟡 A three-locus concatemer is the smallest visually useful multiway example. It is not a claim that Pore-C molecules contain exactly three fragments.
- 🔴 Exact ONT adapter bases and motor attachment are proprietary and are role-labelled.

## Distinction from Hi-C

Hi-C commonly shears a proximity-ligation product and reads two ends, reducing each sequenced molecule to a pair. Pore-C instead keeps a long ligated polymer intact and reads through multiple junctions. Pairwise contacts can be computed later, but the original molecule's multiway membership remains available.
