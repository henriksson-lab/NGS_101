# snHi-C

## 1. What it is

Flyamer et al.'s single-nucleus Hi-C workflow, designed around whole-genome amplification rather than biotin-junction capture.

## 2. Sources and evidence

- 🟢 Flyamer et al., *Nature* 2017, DOI [10.1038/nature21711](https://doi.org/10.1038/nature21711); author manuscript is available as [PMC5639698](https://pmc.ncbi.nlm.nih.gov/articles/PMC5639698/).
- 🟢 The method isolates nuclei before fixation, digests with DpnII, performs proximity ligation without biotin enrichment, and uses phi29 multiple-displacement amplification.
- 🟡 The outer Illumina library shell is reconstructed because exact kit oligos are not printed.
- 🟢 The Supplementary Information (`nature21711_supp_NIHMS852755-supplement-supp_info.pdf.txt`)
  says the embryonic/zygote branch sends sonicated MDA DNA to a sequencing service
  for library preparation without naming its adapter or PCR oligos. Its separate
  K562 branch describes an indexed Illumina adapter and Illumina primers after
  end repair and dA attachment, but likewise prints no sequences. These branch
  descriptions cannot establish that one adapter set was shared by both.

## 3. Distinguishing chemistry

There is no biotin fill or streptavidin bottleneck. All purified ligation products enter MDA, so the schematic keeps contact formation and the later sequencing-library conversion as separate transformations.

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC5639698 JATS title/abstract; original NIHMS852755 supp_info PDF text18–68.

Selected isolated-pronucleus branch; separated whole-zygote branch. Added decrosslinking condition within MDA workflow without inventing internal order, post-MDA cleanup, sonication, second cleanup and paired 125-base reads. Library boundary names VBCF preparation and 300–1300 bp selection; withdrew unverified absence claims.

Deferred: Cited GenomiPhi method31, manufacturer manual and VBCF detailed preparation; exact adapter assignment. This is a capped check, not a full source audit or confirmation of missing information.
