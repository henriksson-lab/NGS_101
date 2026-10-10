# GRO-seq — research notes

Defining source: Core, Waterfall and Lis, *Science* (2008), [10.1126/science.1162228](https://doi.org/10.1126/science.1162228) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC2833333/)).

## 1. What it is

Global run-on sequencing (GRO-seq) maps the positions and orientations of transcriptionally engaged RNA polymerases. It is a nascent-RNA selection chemistry, not ordinary abundance RNA-seq.

## Sources read

- 🟢 The defining article's main text and Methods/Supplementary Methods pointer were read through the open author manuscript.
- 🟢 A later stepwise account, *Global Run-On sequencing (GRO-seq)* ([PMC5880289](https://pmc.ncbi.nlm.nih.gov/articles/PMC5880289/)), was used only to clarify library order where the short defining article points to its supplement.
- The repository source fetch records the publisher copy as manual: `python3 tools/get_sources.py --doi 10.1126/science.1162228`.

## Construct-changing steps

1. 🟢 Isolated nuclei undergo nuclear run-on under conditions that prevent new initiation. Already-engaged polymerases extend their attached RNAs with BrUTP.
2. 🟢 Run-on RNA is kept short and chemically hydrolyzed to about 100 bases.
3. 🟢 BrU-containing RNA is selected three times with an antibody to the nucleotide analogue.
4. 🟢 The original library represents the 5′ end of the fragmented BrU-incorporated RNA. The detailed later protocol polyadenylates RNA, primes first-strand cDNA with an oligo(dT)-containing library adapter, and completes the sequencing library.

The BrU is therefore at the in-vitro-extended RNA 3′ end; it is a selection handle, not a cell barcode or UMI.

## Final library and read layout

- 🟢 The defining experiment used 33-base single-end Illumina reads.
- 🟢 Strand orientation reports polymerase direction.
- 🟡 The public schematic labels the historical sequencing arms by role. It does not replace them with a modern TruSeq sequence, because the defining paper does not establish that identity.

## Open questions

- 🔴 No single current sequencing-primer set can be asserted for every implementation called GRO-seq. A later named revision should be a separate page if its library architecture changes.
