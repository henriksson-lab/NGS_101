# Pore-C — research notes

## Defining source

Deshpande AS, Ulahannan N, Pendleton M, et al. *Identifying synergistic high-order 3D chromatin conformations from genome-scale nanopore concatemer sequencing.* Nature Biotechnology 2022. [doi:10.1038/s41587-022-01289-z](https://doi.org/10.1038/s41587-022-01289-z)

Operational sources read: Oxford Nanopore Technologies, *Restriction Enzyme Pore-C protocol* (v5) and its linked information sheet. These document a later vendor implementation of the method class. The defining article and the cited 2019 preprint remain inaccessible for this audit, so the original paper's sequencing-kit revision is unresolved.

## Evidence record

- 🟢 Formaldehyde fixes spatially proximal chromatin within intact nuclei.
- 🟢 In-situ restriction digestion exposes cohesive ends. The ONT procedure uses NlaIII as an example and its information sheet discusses NlaIII and DpnII; this page draws the NlaIII/CATG branch.
- 🟢 T4 DNA ligase diffuses into the crosslinked nuclear scaffold and joins nearby restriction monomers into chimeric Pore-C polymers. The ONT procedure uses 16 °C for 6 h and warns that prolonged ligation can increase trans contacts.
- 🟢 Proteinase K digestion reverses crosslinks and releases the covalent dsDNA polymers. Phenol/chloroform extraction, ethanol precipitation and long-fragment selection preserve those polymers.
- 🟢 The linked ONT information sheet recommends 0.85X SPRI size selection when extract yield exceeds 1 µg, but sends lower-yield extracts to library preparation without selection. It recommends a native ligation kit when input permits and offers a Rapid PCR Barcoding branch for lower yields. The page draws the native, size-selected branch.
- 🟢 The later ONT SQK-LSK114 kit documentation repairs and dA-tails DNA before ligating a dT-tailed adapter. This establishes the displayed end geometry for that later kit, not the historical kit used in the defining paper. Native preparation is PCR-free and retains DNA modifications.
- 🟡 A three-locus concatemer is the smallest visually useful multiway example. It is not a claim that Pore-C molecules contain exactly three fragments.
- 🔴 The exact sequencing kit and adapter used in the defining 2022 article remain unresolved. Later ONT SQK-LSK114 documentation publishes an adapter duplex, so the page does not claim that ONT adapter bases are universally unpublished. Its bracketed adapter represents the still-unverified historical configuration.

The maintainer source checklist in `catalogue/source_reviews/pore-c__10.1038+s41587-022-01289-z.json` records the inaccessible article and preprint, two independently checked publisher supplements, independently checked later ONT originals, and each linked-information disposition. The historical kit identity remains blocked; no missing-information claim is accepted.

**Source discrepancy:** the current ONT information-sheet overview prints document version `RPC_S1015_v1_revL_12Nov2019`, while its library-preparation section recommends SQK-LSK114 and calls SQK-LSK109/110 older alternatives. We treat that section as later kit guidance; the metadata does not establish which adapter the 2022 article used.

## Distinction from Hi-C

Hi-C commonly shears a proximity-ligation product and reads two ends, reducing each sequenced molecule to a pair. Pore-C instead keeps a long ligated polymer intact and reads through multiple junctions. Pairwise contacts can be computed later, but the original molecule's multiway membership remains available.
