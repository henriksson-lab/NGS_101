# SMART-seq3

Hagemann-Jensen et al. (2020), [doi:10.1038/s41587-020-0497-0](https://doi.org/10.1038/s41587-020-0497-0).

The current page follows the authors' **Version 3**, printed DOI [10.17504/protocols.io.bcq4ivyw](https://doi.org/10.17504/protocols.io.bcq4ivyw), with published Supplementary Table 1 as oligo authority. The printed revision was checked. The defining article HTML provided its identity/access boundary; full published Methods were not available in that view. Superseded preprints do not supply current chemistry.

## Oligos and identifiers

🟢 Native dT primer is `/5Biosg/ACGAGCATCAGCAGCATACGA(T)30VN`; native TSO is `/5Biosg/AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG`. The eight degenerate bases carry the transcript UMI. The three terminal G residues are RNA; the 5-prime biotin belongs to the ordering reagent and is not copied as a base modification.

🟢 Fwd PCR primer is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G`; Rev is `ACGAGCATCAGCAGCATAC*G*A`. The asterisks are two terminal phosphorothioate linkages. Published Table 1 and V.3 agree on the corresponding bases.

## Native molecular route

RT creates first-strand and then template-switched single-strand cDNA; asymmetric KAPA amplification consumes that actual substrate. The selected optional 0.6:1 PEG bead cleanup precedes tagmentation. Nextera XT tagmentation is followed by SDS stripping, the first Phusion extension, twelve index PCR cycles, pooling and 0.6:1 PEG bead purification.

Retained-TSO molecules carry the 11-base tag and original eight-base UMI. Internal molecules are a separate library class without that UMI. Both classes are co-purified; they are not physical fractions isolated from one another. Other native end classes retain unresolved recovery rather than being discarded merely for nonamplification.

The author protocol recommends custom 10-base or 8-base indexes. A particular native member or width is unassigned here. The endpoint retains that uncertainty; generic Nextera eight-base primers are not treated as native reagent identity. V.3 permits single-end or paired-end acquisition; run-specific read lengths/index reads remain unassigned.

## Deferred source and model work

The exact native transferred adapter geometry, first-cycle PCR relations, custom index member and run assignment need follow-up. Optional final gel selection is outside the chosen branch. This capped pass does not establish a complete chemistry review or accepted missing-information claim.

## Critical implementation follow-up (2026-10-18)

Retained the independent RNA template beside first-strand and switched cDNA. The PCR output includes an RNA-derived molecular population with unresolved thermal fate and later cleanup recovery; this conceptual class distinction does not assert physical isolation. Copied TSO orientation was rechecked against complement_segments: it already reverses order and complements bases while preserving UMI feature identity.
