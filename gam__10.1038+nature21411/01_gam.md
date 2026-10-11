# GAM — research notes

Beagrie et al., *Nature* (2017),
[doi:10.1038/nature21411](https://doi.org/10.1038/nature21411).

## Sources read

- 🟢 Defining article, Methods and Supplementary Guide (PMCID: PMC5366070).
- 🟢 Primary sequencing deposit GSE64881.
- 🟢 The later Multiplex-GAM paper
  ([doi:10.1038/s41592-023-01903-1](https://doi.org/10.1038/s41592-023-01903-1))
  was read only to separate the original one-profile-per-library protocol from the later
  multi-profile variant.

## What constitutes one GAM observation

GAM does not ligate spatial neighbours. A random approximately 220-nm plane through one
fixed nucleus selects many physically separate pieces of chromosomes. One column of the
eventual segregation matrix is the set of genomic windows detected in that nuclear profile
(NP); contacts are inferred because nearby loci co-segregate across many independent NPs.

1. 🟢 Mouse ESC nuclear architecture is strongly fixed, embedded in 2.1 M sucrose, frozen
   and cut into approximately 220-nm ultrathin cryosections.
2. 🟢 Sections are placed on PEN membrane slides. One visually identified NP is cut and
   laser-catapulted into an adhesive PCR cap. One cap in each eight-cap strip is empty and
   follows the entire WGA procedure as a contamination control.
3. 🟢 DNA is extracted and whole-genome amplified with the Sigma WGA4 single-cell kit.
   The defining methods disclose the temperatures, times and reagent roles but not the
   proprietary adapter bases.
4. 🟢 WGA material is purified. A TruSeq DNA HT library uses up to 1.1 µg and is Pippin
   selected to 300–500 nt; the alternative TruSeq Nano DNA HT workflow uses up to 200 ng
   and a bead selection centred near 350 nt.
5. 🟢 Ninety-six indexed NP libraries are pooled. The original mESC data comprise 408
   passing NPs after excluding samples with less than 15% uniquely mapped mouse reads.

## Final molecule and read entry

- 🟢 WGA leaves a 30-base amplification adapter on both ends of each genomic piece.
- 🟢 The final sequencing shell is the named TruSeq DNA HT/Nano chemistry. Its sample
  indexes identify an NP library, not a molecular contact or cell barcode.
- 🟢 Pools were sequenced single-end on HiSeq in rapid-run mode. The first 30 Read 1
  cycles traverse WGA adapter and were deliberately not imaged (dark cycles), followed by
  100 imaged bases entering genomic DNA. There is no Read 2 observation in this defining
  run.
- 🔴 WGA4 adapter base sequences are proprietary. The schematic preserves their published
  30-nt length and marks the bases inferred/unknown rather than inventing an oligo.

## Protocol/version boundary

The later Multiplex-GAM protocol can put several NPs into one library and uses Nextera;
later GAM applications also introduced a public GAT-7N/GAM-COM DeepVent WGA. Those are
useful developments, but neither is silently substituted for the original 2017 GAM page.
