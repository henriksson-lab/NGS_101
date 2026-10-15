# sn-m3C-seq

## 1. What it is

A single-nucleus multi-omic protocol that reads chromosome contacts and cytosine methylation from the same bisulfite-converted molecules.

## 2. Sources and evidence

- 🟢 Lee et al., *Nature Methods* 2019, DOI [10.1038/s41592-019-0547-z](https://doi.org/10.1038/s41592-019-0547-z); full text [PMC6765423](https://pmc.ncbi.nlm.nih.gov/articles/PMC6765423/).
- 🟢 Proximity-ligation products are carried into bisulfite conversion without a biotin-junction pull-down.
- 🟢 Library construction follows the snmC-seq2 strategy: random priming with a
  partial Read 1 tail, inline index and 3′ H9 arm, followed by adapter tagging and
  PCR. The inherited oligo is printed in Luo et al. Supplementary Methods under
  “Random primers” (for example, P5L_AD002_H); the P5/i5 and P7/i7 arms are added
  by the later indexed PCR primers. The sn-m3C article identifies this predecessor
  procedure in Methods under “Bulk and single-cell m3C-seq library preparation.”
- 🟢 The cited snmC-seq2 source (Luo et al., DOI
  [10.1038/s41467-018-06355-2](https://doi.org/10.1038/s41467-018-06355-2))
  prints eight random-primer designs in its Supplementary Methods and 768 paired
  P5/i5 and P7/i7 amplification-primer pairs in Supplementary Table 2. Those are
  available in the separate [snmC-seq2 note](../snmc-seq2__10.1038+s41467-018-06355-2/01_snmc-seq2.md).
- 🟡 The page draws a representative inherited snmC-seq2 endpoint using the
  source-printed P5L_AD002_H inline barcode and dual-index PCR pair 1. The
  sn-m3C paper does not specify which members of those published families were
  used for each library. The representative must not be read as that sample's
  actual barcode or index identity. 🔴 The proprietary Adaptase-added bases
  remain unpublished; converted contact insert bases are sample-dependent.
- 🟢 The authors trim the first 25 bases and final 3 bases of each read to remove
  random-primer sequence and the low-complexity Adaptase tail. These are read
  processing lengths; the drawn random H9 and schematic Adaptase span must not
  be read as the full lengths of those trimmed read regions.

## 3. Distinguishing chemistry

Bisulfite conversion happens after the contact molecule has been created. Consequently each paired read supports both a ligation contact and local methylation inference, unlike a conventional Hi-C read pair.
