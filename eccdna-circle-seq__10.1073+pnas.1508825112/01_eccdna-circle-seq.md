# eccDNA Circle-Seq — research notes

## Sources read

- 🟢 Møller et al. 2015, DOI [10.1073/pnas.1508825112](https://doi.org/10.1073/pnas.1508825112), PMC4475933, full text and Methods.

## Molecular order

🟢 Circular DNA is purified from yeast. NotI cuts chromosomal DNA to expose more ends,
then Plasmid-Safe ATP-dependent DNase is replenished over 130 hours to exhaustively remove
linear double-stranded DNA. One percent of the enriched circular fraction is amplified for
16 hours with phi29 polymerase. The rolling-circle product is sonicated to about 300 nt,
indexed on an Apollo 324 system and sequenced single-end on HiSeq 2000.

## Boundary and uncertainty

The circles are biological inputs, not library-prep intermediates. 🔴 The paper names the
Apollo 324 adapter/index workflow but does not print its adapter sequences; the public
page marks conventional Illumina-arm geometry as inferred and does not fabricate bases.
This is separate from CRISPR CIRCLE-seq.
