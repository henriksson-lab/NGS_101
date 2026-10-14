# sf4sU-Seq

**Defining paper:** Gerlevik et al., *Nucleic Acids Research* (2026), doi:10.1093/nar/gkag652.

## Molecular boundary

🟢 Pulse-labelled 4sU RNA is thiol-biotinylated and affinity-purified. The enriched RNA is resolved on a denaturing gel and an approximately 100-nt fraction is recovered.

🟢 Directional sequencing preserves both native ends of the selected short RNA population, relating initiation to promoter-proximal pausing.

## Evidence boundary

🟡 sf4sU-Seq and SNU-Seq share a defining paper but remain distinct protocol pages because size selection and bacterial polyadenylation create different molecular endpoints.

## Oligo source audit

🟢 The defining article (`PMC13309788_PMC13309788.xml.txt`, line 124) and its Supplementary Methods (`PMC13309788_supp_gkag652_supplemental_file.pdf.txt`, lines 440–452) specify 5′ decapping with RNA 5′ pyrophosphohydrolase followed by the NEBNext Small RNA Library kit for Illumina. The Results section places 3′ linker ligation before 5′ linker ligation. The supplement prints no kit oligo bases, so the 3′/5′ linkers, RT primer, PCR/index primers and sequencing entry remain kit-specific rather than source-verified bases in this page.
