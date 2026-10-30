# NICHE-seq — naive lymph-node niche selection

Medaglia et al., *Spatial reconstruction of immune niches by combining photoactivatable reporters and scRNA-seq*, [DOI10.1126/science.aao4277](https://doi.org/10.1126/science.aao4277).

## Inspected source and boundary

🟢 This capped pass inspected the associated author manuscript PMC7234837/EMS86422: identity, main47/49/51, references7/17 and supplement identities. Its header explicitly identifies a Funders Author Manuscript. It is not asserted fully equivalent to the final published Methods. The mandatory scraper ran after conversion of this original; no oligo hits supply any sequencing-adapter evidence. Zero hits are not an absence audit.

🟢 NICHE-seq combines protein-reporter photoactivation, two-photon microscopy, physical cell isolation and MARS-seq. The changed physical preparation boundary is recovery of cells from an optically selected tissue niche before library preparation. It is distinct from activating an oligo barcode or merely selecting cells computationally after sequencing.

The private partial graph selects naive inguinal lymph-node niches. It starts with prepared PA-GFP host tissue containing transferred TdTomato+ T-cell and CFP+ B-cell landmarks. Upstream animal preparation and transfer operations are outside this graph; it does not model unrelated infection, spleen or tumor application cohorts. The source's B-cell-follicle and T-cell-area selections remain alternative regions within this method, not separate protocol entries.

## Source-qualified physical route

🟢 Main49 identifies one wavelength for local photoactivation and a different subsequent imaging wavelength. The model declares and checks both literal source values. This is photochemical activation of a protein reporter, not incorporation of a spatial DNA index. Duration, dose, exact optical geometry and perfect confinement of activation are unassembled; the source's measured enrichment is not a guarantee for every cell.

🟢 Main47/51 establishes tissue dissociation followed by fluorescence-based FACS collection of photoactivated single cells. Dissociation is an opaque boundary: the exact enzymes, physical fractions, viability gates and recovery are not borrowed from another protocol. The graph preserves a physical unselected-cell remainder with an unassigned later fate. Sequencing filters are not represented as removal of cells or molecules.

🟢 MARS-seq follows physical cell collection. Main47 explicitly identifies the method and reference7 identifies Jaitin2014, DOI10.1126/science.1247651. The graph represents that delegated preparation as an opaque component with contributing molecular stocks. The original MARS paper and detailed local supplementary Methods were not read here; no constructor, current version, oligos, UMI pattern, index/read cycles or complete library molecule is inherited from similarity or citation alone.

## Open scope

🟢 Reference17 explicitly locates Materials and methods in the supplementary material. The actual linked Supplementary Materials incl Figures PDF is gated through the standard PMC route and unavailable as a Cloud TDM object; it remains needed for detailed preparation. The separate Supplementary Tables workbook is also gated, but its chemistry relevance is unreviewed: it is not assumed to be a primer table or requested as an invented missing disclosure. Full final-version applicability remains unverified.

One original was selected before opening and inspected; the source investigation is closed. This is a source-positive frontend description and private partial graph, not a finished schematic, complete library validation or accepted missing-information claim. The selected independent implementation recheck is recorded separately from the two-reader absence gate.
