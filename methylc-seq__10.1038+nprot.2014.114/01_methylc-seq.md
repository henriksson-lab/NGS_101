# MethylC-seq / conventional WGBS — research notes

## 1. What it is

MethylC-seq is ligation-first whole-genome bisulfite sequencing: completed, protected
adapter-tagged genomic fragments undergo sodium-bisulfite conversion and limited PCR.

## 2. Sources read

- 🟢 Urich et al. 2015, DOI [10.1038/nprot.2014.114](https://doi.org/10.1038/nprot.2014.114),
  complete protocol text (PMC4465251).
- 🟢 Bioo Scientific / NEXTflex product documentation for Bisulfite-Seq Barcodes,
  including the single 6-nt index architecture and matching PCR-primer sequences.
- 🟢 Illumina Adapter Sequences, document 1000000002694, for the canonical sequencing
  primer sequences used to identify the compatible sites computationally.

## 3. Molecular order

🟢 Genomic DNA is sonicated to about 200 bp, end-repaired and 5′ phosphorylated, given 3′
dAMP overhangs, and ligated to adapters in which every cytosine is methylated. 🟢 The
library is denatured; bisulfite converts unprotected C to U and leaves 5mC as C. The two
converted strands are no longer complementary. A uracil-tolerant polymerase restores a
duplex library during limited PCR.

## 4. Evidence boundary

🟢 The protocol specifies **NEXTflex Bisulfite-Seq Barcodes–12, Bioo Scientific cat.
no. 511912**, including the corresponding PCR primers. The product is a methylated,
single-index Illumina library with a 6-nt i7 barcode. Its final construct contains the
standard TruSeq-compatible Read 1, Index 1 and Read 2 primer sites; it has no i5 index.

The vendor name matters here: the amplification oligos were the matching **NEXTflex PCR
primers**, not an Illumina-branded TruSeq kit. At the sequence/flow-cell level, however,
the resulting endpoint uses the ordinary TruSeq primer geometry, so standard TruSeq
sequencing primers are the correct ones to draw.
