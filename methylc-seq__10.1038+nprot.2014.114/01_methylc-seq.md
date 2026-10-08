# MethylC-seq / conventional WGBS — research notes

## 1. What it is

MethylC-seq is ligation-first whole-genome bisulfite sequencing: completed, protected
adapter-tagged genomic fragments undergo sodium-bisulfite conversion and limited PCR.

## 2. Sources read

- 🟢 Urich et al. 2015, DOI [10.1038/nprot.2014.114](https://doi.org/10.1038/nprot.2014.114),
  complete protocol text (PMC4465251).

## 3. Molecular order

🟢 Genomic DNA is sonicated to about 200 bp, end-repaired and 5′ phosphorylated, given 3′
dAMP overhangs, and ligated to adapters in which every cytosine is methylated. 🟢 The
library is denatured; bisulfite converts unprotected C to U and leaves 5mC as C. The two
converted strands are no longer complementary. A uracil-tolerant polymerase restores a
duplex library during limited PCR.

## 4. Evidence boundary

🟢 The protocol establishes the methylated Y-adapter requirement and reaction order.
🟡 The page uses dotted canonical TruSeq arms to expose run geometry because the adapter
product, not an exact current index-kit sequence, is the relevant invariant.
