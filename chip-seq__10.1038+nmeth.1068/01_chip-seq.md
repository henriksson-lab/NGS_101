# Conventional ChIP-seq — research notes

## 1. What it is

ChIP-seq enriches crosslinked chromatin fragments associated with an antibody target and
turns the recovered DNA into a conventional end-repaired, dA-tailed, adapter-ligated
sequencing library.

## 2. Sources read

- 🟢 Robertson et al. 2007, DOI [10.1038/nmeth.1068](https://doi.org/10.1038/nmeth.1068),
  the defining genome-wide ChIP sequencing report.
- 🟢 Schmidt et al. 2009, *Methods* 48:240–248,
  DOI [`10.1016/j.ymeth.2009.03.001`](https://doi.org/10.1016/j.ymeth.2009.03.001)
  (PMC4052679), for the explicit end-repair, dA and adapter-ligation sequence.

## 3. Molecular order

🟢 Formaldehyde-crosslinked chromatin is sonicated and immunoprecipitated with a
target-specific antibody. Crosslinks are reversed and DNA recovered. The sequencing
workflow repairs ends, phosphorylates them, adds a single 3′ dA, ligates an asymmetric
Illumina adapter, size-selects and PCR-amplifies.

## 4. Historical library identity

🟢 The contemporary protocol identifies Illumina Genomic Adapter Oligo Mix (part
1000521) and Genomic PCR primers 1.1 and 2.1. Illumina's archived Genomic DNA Sample
Prep documentation and adapter-sequence guide publish the corresponding single-read
oligos. The 2007 experiment produced unindexed single-end reads on the Solexa Genome
Analyzer; representing it as paired-end, indexed TruSeq would be anachronistic.

## Bounded curation round 2 (2026-10-16)

Inspected: Local source bundle MANIFEST.tsv inspected; existing get_sources --slugs chip-seq__10.1038+nmeth.1068 attempted:0 originals returned.

No correction made: defining original and historical vendor originals were not available in this bundle. No source claim was re-verified.

Deferred: Defining Robertson2007 article, Schmidt2009 protocol and archived Illumina original11257047/adapter guide; retrieval yielded no source and remains deferred. This capped pass does not establish a complete source review or publication-wide absence.
