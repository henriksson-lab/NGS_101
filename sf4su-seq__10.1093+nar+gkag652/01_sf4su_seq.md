# sf4sU-Seq

**Defining paper:** Gerlevik et al., *Nucleic Acids Research* (2026), doi:10.1093/nar/gkag652.

## Molecular boundary

🟢 Pulse-labelled 4sU RNA is thiol-biotinylated and affinity-purified. The enriched RNA is resolved on a denaturing gel and an approximately 100-nt fraction is recovered.

🟢 Directional sequencing preserves both native ends of the selected short RNA population, relating initiation to promoter-proximal pausing.

## Evidence boundary

🟡 sf4sU-Seq and SNU-Seq share a defining paper but remain distinct protocol pages because size selection and bacterial polyadenylation create different molecular endpoints.

## Oligo source audit

🟢 The defining article (`PMC13309788_PMC13309788.xml.txt`, line 124) and its Supplementary Methods (`PMC13309788_supp_gkag652_supplemental_file.pdf.txt`, lines 440–452) specify 5′ decapping with RNA 5′ pyrophosphohydrolase followed by the NEBNext Small RNA Library kit for Illumina. The Results section places 3′ linker ligation before 5′ linker ligation. The supplement's Kits and Reagents table identifies **NEB E7300S** (`PMC13309788_supp_gkag652_supplemental_file.pdf.txt`, lines 725–729).

🟢 The [official NEB E7300 manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7300_e7330_e7560_e7580.pdf), version 9.0_07/25, page 16 (`manuale7300_e7330_e7560_e7580.pdf.txt`, lines 1041–1067), gives the kit molecules verbatim:

| Reagent | Oligo as written, 5′→3′ |
|---|---|
| 3′ SR adapter | `rAppAGATCGGAAGAGCACACGTCT-NH2` |
| 5′ SR adapter | `rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC` |
| SR RT primer | `AGACGTGTGCTCTTCCGATCT` |
| SR PCR primer | `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A` |
| indexed PCR primer family | `CAAGCAGAAGACGGCATACGAGAT[6-nt index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T` |

🟡 The manual establishes the E7300 reagent family, including the 5′ activated and 3′ blocked ends of the 3′ adapter. The paper does not identify its individual index tube or lot, so the six index bases stay unknown. The graph now declares RppH decapping and the separate 3′ and 5′ single-strand ligations; the final amplicon and sequencing-primer landing sites remain role-level until they are independently reconciled base by base.
