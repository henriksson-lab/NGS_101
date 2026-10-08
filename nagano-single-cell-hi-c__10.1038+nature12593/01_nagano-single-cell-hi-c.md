# Nagano single-cell Hi-C

## 1. What it is

The foundational single-cell Hi-C implementation: proximity ligation is performed in a population of fixed nuclei, after which individual nuclei are isolated and their contact libraries are prepared separately.

## 2. Sources and evidence

- 🟢 Nagano et al., *Nature* 2013, DOI [10.1038/nature12593](https://doi.org/10.1038/nature12593), defining paper.
- 🟢 Nagano et al., *Nature Protocols* 2015, DOI [10.1038/nprot.2015.127](https://doi.org/10.1038/nprot.2015.127), detailed implementation used to resolve BglII, biotin-14-dATP, streptavidin capture and PCR ordering.
- 🟢 The 2013 supplementary methods print twelve customized Illumina adapters, each with
  a 3-bp identification tag, and both PCR primers. The schematic uses the printed CAA
  adapter member and the historical paired-end primer sites; there is no invented i7 read.

## 3. Distinguishing chemistry

🟢 BglII digestion and biotinylated fill-in precede proximity ligation. 🟢 Single nuclei are separated after the bulk ligation reaction. 🟢 Following crosslink reversal and fragmentation, streptavidin enrichment selects marked contact junctions for PCR library construction.
