# sn-m3C-seq

## 1. What it is

A single-nucleus multi-omic protocol that reads chromosome contacts and cytosine methylation from the same bisulfite-converted molecules.

## 2. Sources and evidence

- 🟢 Lee et al., *Nature Methods* 2019, DOI [10.1038/s41592-019-0547-z](https://doi.org/10.1038/s41592-019-0547-z); full text [PMC6765423](https://pmc.ncbi.nlm.nih.gov/articles/PMC6765423/).
- 🟢 Proximity-ligation products are carried into bisulfite conversion without a biotin-junction pull-down.
- 🟢 Library construction follows the snmC-seq2 strategy: indexed P5-bearing random priming, adapter tagging and PCR.
- 🟢 The cited snmC-seq2 source (Luo et al., DOI
  [10.1038/s41467-018-06355-2](https://doi.org/10.1038/s41467-018-06355-2))
  prints eight random-primer designs in its Supplementary Methods and 768 paired
  P5/i5 and P7/i7 amplification-primer pairs in Supplementary Table 2. Those are
  available in the separate [snmC-seq2 note](../snmc-seq2__10.1038+s41467-018-06355-2/01_snmc-seq2.md).
- 🔴 The defining sn-m3C source does not print a complete outer adapter/primer set. The
  paper does not specify which predecessor primer/barcode members were used for each
  sn-m3C library. The page keeps those arms unresolved rather than assigning canonical TruSeq; converted
  insert bases are sample-dependent and are not falsely written as a fixed sequence.
- 🟢 The authors trim the first 25 bases and final 3 bases of each read to remove random-primer sequence and the low-complexity Adaptase tail; both regions remain explicit in the model.

## 3. Distinguishing chemistry

Bisulfite conversion happens after the contact molecule has been created. Consequently each paired read supports both a ligation contact and local methylation inference, unlike a conventional Hi-C read pair.
