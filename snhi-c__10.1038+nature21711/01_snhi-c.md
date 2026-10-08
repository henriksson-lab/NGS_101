# snHi-C

## 1. What it is

Flyamer et al.'s single-nucleus Hi-C workflow, designed around whole-genome amplification rather than biotin-junction capture.

## 2. Sources and evidence

- 🟢 Flyamer et al., *Nature* 2017, DOI [10.1038/nature21711](https://doi.org/10.1038/nature21711); author manuscript is available as [PMC5639698](https://pmc.ncbi.nlm.nih.gov/articles/PMC5639698/).
- 🟢 The method isolates nuclei before fixation, digests with DpnII, performs proximity ligation without biotin enrichment, and uses phi29 multiple-displacement amplification.
- 🟡 The outer Illumina library shell is reconstructed because exact kit oligos are not printed.

## 3. Distinguishing chemistry

There is no biotin fill or streptavidin bottleneck. All purified ligation products enter MDA, so the schematic keeps contact formation and the later sequencing-library conversion as separate transformations.
