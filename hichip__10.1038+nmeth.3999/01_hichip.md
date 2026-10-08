# HiChIP — research notes

## 1. What it is

HiChIP combines in situ proximity ligation, protein-directed immunoprecipitation and
on-bead Tn5 library construction to enrich chromatin contacts involving a chosen target.

## 2. Sources read

- 🟢 Mumbach et al. 2016, DOI [10.1038/nmeth.3999](https://doi.org/10.1038/nmeth.3999),
  complete author manuscript and Online Methods (PMC5501173).

## 3. Molecular order

🟢 Crosslinked chromatin is restriction-digested in nuclei, ends are filled with biotin,
and spatially adjacent ends are ligated. The chromatin is then sheared and
immunoprecipitated. 🟢 Tn5 tagmentation occurs while the enriched chromatin remains on
beads. Streptavidin capture retains biotin-containing contact products before PCR.

## 4. Evidence boundary

🟢 The defining implementation uses MboI; the page derives its filled GATCGATC contact
scar and top-strand biotin coordinate from the shared restriction model. The final adapter
architecture is the canonical Nextera molecule implied by Tn5 library construction.
