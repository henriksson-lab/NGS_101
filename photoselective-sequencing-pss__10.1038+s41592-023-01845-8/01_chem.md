# Photoselective sequencing: source components

Mangiameli et al. (2023), **Photoselective sequencing: microscopically guided genomic measurements with subcellular resolution**, DOI10.1038/s41592-023-01845-8.

## Source scope and actual access limit

🟢 Published Nature main HTML identity/abstract inspected: the received page offers subscription access rather than the core Methods. Published Supplementary Information MOESM1 inspected: original p2 oligo table, ExtendedDataFig1A caption, Methods156–405 selected transposase purification, photoselection and blocking-control sections. Two originals, third reserved; source investigation closed. Mandatory scraper ran. The automatically attempted 2022 preprint is superseded and unused. This is not a full source or missing-information audit.

## Positive preparation boundary

🟢 The published abstract describes an in-situ photocaged genomic fragment library and ROI illumination that selects sequencing-competent material. SI PSS1/2 are photocleavable tagmentation adapters. PSS3 is a shortened, phosphorothioate-containing 5′-phosphorylated mosaic-end oligo ending in ddC. PSS4/5 and PSS6/7 are ligation handles and blocked splints; PSS8/9 are library-amplification primers. This preparation gates genomic fragments rather than RT-derived RNA/cDNA substrates. SI distinguishes unbiased cell libraries and chromatin-accessibility libraries in tissues; the core recipe that establishes each route is still needed.

## Verbatim whole orderings, SI p2

🟢 PSS1: `/5Alex546N//iSpPC/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`.
PSS2: `/5Alex546N//iSpPC/GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`.
PSS3: `/5phos/C*T*G*T*C*T*C*T*T*A*T*A*C*A*/3ddC/`.
PSS4: `AATGATACGGCGACCACCGAGATCTACACCTCTCTAT`.
PSS5-indN: `CAAGCAGAAGACGGCATACGAGATBBBBBBBBCGACTGACTCT`.
PSS6: `TGCCGACGAATAGAGAGGTGTAG/3ddC/`.
PSS7: `CCCACGAGACAGAGTCAGTCG/3ddC/`.
PSS8: `AATGATACGGCGACCACCGAGAT`.
PSS9: `CAAGCAGAAGACGGCATACGAGAT`.

🟢 The caption defines B as a known eight-base sample barcode, rather than an IUPAC exclusion class. It directs readers to SupplementaryTable2 for the full set of indexed P7 adapters. That workbook is downloaded but uninspected; no individual index is selected. Terminal `/3ddC/` is an additional modified C nucleotide, retained in the molecular inventory. The first fourteen bases of PSS3 carry the printed phosphorothioate marks; the terminal ddC is the fifteenth paired base.

🟡 Local Scene pairing aligns the fifteen-base modified short strand to the source PSS1 or PSS2 mosaic-end suffix. This is an inventory pairing, not proof of loaded-transposase stoichiometry, insertion chemistry, optical yield or full ligation geometry. The separate `ligation_interface()` helper pairs the source PSS6/PSS7 splint arms to cleaved PSS1/PSS2 handles and PSS4/PSS5 tails: 9+14 or 10+11 paired bases, respectively. It computes an adjacent nick without constructing a sealed product. The final modified ddC is retained outside the asserted paired arms; upstream adapter-prefix/barcode portions are explicitly unpaired in this local view, and no complete fraying/conformation is asserted. Source handles are imported from canonical sequence constants only after whole source-ordering fixtures.

## Open scope

Actual published core Methods/PDF is needed for ordered preparation, enzyme/reagent allocations, end handling and physical fractions. Therefore no whole protocol graph or sequencing-ready product is invented. Exact barcode workbook, whole sealed library/read layout, linked originals and kit details remain unreviewed, not confirmed absent. The blocking controls in SI357–405 are distinct analytical experiments and are not silently substituted for the production recipe. A private component helper and source-ordering tests accompany these notes.
