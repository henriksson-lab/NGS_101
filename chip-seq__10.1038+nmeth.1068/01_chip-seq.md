# Conventional ChIP-seq — research notes

## 1. What it is

ChIP-seq enriches crosslinked chromatin fragments associated with an antibody target and
turns the recovered DNA into a conventional end-repaired, dA-tailed, adapter-ligated
sequencing library.

## 2. Sources read

- 🟢 Robertson et al. 2007, DOI [10.1038/nmeth.1068](https://doi.org/10.1038/nmeth.1068),
  the defining genome-wide ChIP sequencing report.
- 🟢 Park 2009 protocol/review, DOI `10.1038/nprot.2009.12`, complete text
  (PMC4052679), for the explicit end-repair, dA and adapter-ligation sequence.

## 3. Molecular order

🟢 Formaldehyde-crosslinked chromatin is sonicated and immunoprecipitated with a
target-specific antibody. Crosslinks are reversed and DNA recovered. The sequencing
workflow repairs ends, phosphorylates them, adds a single 3′ dA, ligates an asymmetric
Illumina adapter, size-selects and PCR-amplifies.

## 4. Evidence boundary

🟢 Reaction order is explicit in the protocol. 🟡 The old proprietary Illumina adapter is
shown as a dotted canonical paired-end TruSeq layout so primer geometry is visible without
claiming that current adapter bases were printed in the 2007 paper.
