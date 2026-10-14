# HiChIP — research notes

## 1. What it is

HiChIP combines in situ proximity ligation, protein-directed immunoprecipitation and
on-bead Tn5 library construction to enrich chromatin contacts involving a chosen target.

## 2. Sources read

- 🟢 Mumbach et al. 2016, DOI [10.1038/nmeth.3999](https://doi.org/10.1038/nmeth.3999),
  complete author manuscript and Online Methods (PMC5501173).
- 🟢 The defining paper's supplementary HiChIP protocol gives Nextera Ad1.1 (universal) and Ad2.x (barcoded) as the first PCR primers; the Online Methods calls them Ad1_noMX and Ad2.X. Supplementary Table 4 publishes the full Ad1_noMX `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG` and all 24 Ad2 indexed members; for example Ad2.1_TAAGGCGA is `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT`. The Online Methods also prints the optional post-pooling short primers verbatim: Nextera_i7Short `CAAGCAGAAGACGGCATACGAGAT` and Nextera_i5Short `AATGATACGGCGACCACCGAGATCTACAC`. These short primers perform additional amplification after indexed primer incorporation; they do not introduce the indexes.

## 3. Molecular order

🟢 Crosslinked chromatin is restriction-digested in nuclei, ends are filled with biotin,
and spatially adjacent ends are ligated. The chromatin is then sheared and
immunoprecipitated. 🟢 Tn5 tagmentation occurs while the enriched chromatin remains on
beads. Streptavidin capture retains biotin-containing contact products before PCR.

## 4. Evidence boundary

🟢 The defining implementation uses MboI; the page derives its filled GATCGATC contact
scar and top-strand biotin coordinate from the shared restriction model. The final adapter
architecture is the canonical Nextera molecule implied by Tn5 library construction.
🟢 Supplementary Table 4 lists the complete Ad2.1–Ad2.24 primer family. The selected
Ad2.X member varies by sample; the published short-primer pair is fixed across samples.
