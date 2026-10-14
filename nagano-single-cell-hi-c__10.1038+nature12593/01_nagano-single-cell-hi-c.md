# Nagano single-cell Hi-C

## 1. What it is

The foundational single-cell Hi-C implementation: proximity ligation is performed in a population of fixed nuclei, after which individual nuclei are isolated and their contact libraries are prepared separately.

## 2. Sources and evidence

- 🟢 Nagano et al., *Nature* 2013, DOI [10.1038/nature12593](https://doi.org/10.1038/nature12593), defining paper.
- 🟢 Nagano et al., *Nature Protocols* 2015, DOI [10.1038/nprot.2015.127](https://doi.org/10.1038/nprot.2015.127), detailed implementation used to resolve BglII, biotin-14-dATP, streptavidin capture and PCR ordering.
- 🟢 The [2013 supplementary methods](https://www.nature.com/nature/journal/v502/n7469/extref/nature12593-s1.pdf),
  pp. 3–5, print twelve customized Illumina adapters, each with a 3-bp identification
  tag, and both PCR primers. The schematic uses the printed CAA adapter member and the
  historical paired-end primer sites; there is no invented i7 read.

The ordered CAA adapter is 🟢 `5′-ACACTCTTTCCCTACACGACGCTCTTCCGATCTCAA*T-3′`
paired with `5′-pTTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG-3′`. Here `p` is a 5′
phosphate and `*` is a terminal phosphorothioate linkage, as defined by the supplement.
These modifications belong to the ligation oligos; the final PCR product is drawn as
ordinary DNA. The two printed library-amplification primers also end with a terminal
phosphorothioate linkage. The paper reports 2 × 40 bp reads beginning with the inline
identification tag, with no separate index read.

## 3. Distinguishing chemistry

🟢 BglII digestion and biotinylated fill-in precede proximity ligation. 🟢 Single nuclei are separated after the bulk ligation reaction. 🟢 Following crosslink reversal and fragmentation, streptavidin enrichment selects marked contact junctions for PCR library construction.
