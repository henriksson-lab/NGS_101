# UDiTaS — research notes

Giannoukos et al., *BMC Genomics* (2018), [doi:10.1186/s12864-018-4561-9](https://doi.org/10.1186/s12864-018-4561-9), PMCID: PMC5861650.

## Sources read

- 🟢 Full open-access article and Methods.
- 🟢 Additional file 2, the detailed UDiTaS protocol and oligonucleotide table.

## Construct logic

1. 🟢 A custom hyperactive Tn5 transposome carries P5, one of eight i5 pooling barcodes, a 10-base random UMI, Nextera S5 and the mosaic end. Its bottom oligo is 5′ phosphorylated and ends in ddC.
2. 🟢 Tagmentation simultaneously fragments genomic DNA and marks random endpoints. Unlike a fixed pair of locus primers, the random endpoint lets molecules with short indels, long deletions, inversions and translocation junctions enter the same assay.
3. 🟢 Round-one PCR pairs the common P5 primer with one target-specific anchor carrying the full TruSeq Read 2 sequence at its 5′ end.
4. 🟢 Round-two PCR uses P5 again and a P7–i7 primer. The final molecule is dual-indexed: i7 is a conventional sample index, while the opposite index read traverses the 10-nt UMI and then the eight-base i5 pooling barcode.
5. 🟢 The defining experiments were sequenced on MiSeq. The protocol run sheet specifies Read 1 at 300 bases, Index 1 at eight bases and Index 2 at 18 bases, with no Read 2. The repository derives barcode/UMI cycles from the final construct and primer directions rather than copying them into a separate layout.

## Represented branch

The page uses the published N501 Tn5 top oligo, N701 P7 primer and CEP290 guide-323 anchor as concrete examples. Supplementary Table 1 publishes eight Tn5 pooling-index variants, multiple locus anchors and a set of i7 primers. These alternatives change index or target sequence, not the library architecture.
