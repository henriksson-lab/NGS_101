# Stevens single-cell Hi-C

## 1. What it is

A single-cell Hi-C workflow in which nuclei are isolated before their individual restriction, ligation and amplification reactions.

## 2. Sources and evidence

- 🟢 Stevens et al., *Nature* 2017, DOI [10.1038/nature21429](https://doi.org/10.1038/nature21429), defining paper and methods.
- 🟢 The contact digest is **MboI**, followed by biotin-14-dATP fill-in and proximity
  ligation. AluI is a later fragmentation step on captured contact DNA.
- 🟢 The supplement prints the customized 3-bp-tagged Illumina adapters and both PCR
  primers. It also reports Nextera XT dual-index preparation as an alternative route.

## 3. Distinguishing chemistry

The essential distinction from Nagano 2013 is early isolation of individual nuclei. The
contact junction is the biotin-filled MboI junction; AluI must not be drawn there.
