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

## Later published GRO-seq library oligos

🟢 The later author protocol PMCID `PMC5880289`, **Materials, Oligonucleotides**
(JATS text lines 107–115), prints the following ordered oligos. The internal
`/idSp/` is a stable abasic spacer; `V` and `N` are degenerate bases.

| Role | 5′→3′ sequence as printed |
|---|---|
| NTI223 poly(dT) first-strand/adapter oligo | `/5Phos/GAT CGT CGG ACT GTA GAA CTC T/idSp/CA AGC AGA AGA CGG CAT ACG ATT TTT TTT TTT TTT TTT TTT VN` |
| NTI200 PCR primer | `CAA GCA GAA GAC GGC ATA` |
| NTI201 PCR primer | `AAT GAT ACG GCG ACC ACC GAC AGG TTC AGA GTT CTA CAG TCC GAC G` |
| NTI202 sequencing primer | `CGACAGGTTCAGAGTTCTACAGTCCGACGATC` |

🟡 This is a source-verified **later implementation**, not an identity claim
for the original 2008 Science library. The short defining author manuscript
does not print its library oligos, and its Science supplement was not available
through the source fetcher. Its historical 33-base run therefore remains
role-level on the public page. The later 2018 primer family is a concrete
starting point only if that implementation receives its own schematic.

## Final library and read layout

- 🟢 The defining experiment used 33-base single-end Illumina reads.
- 🟢 Strand orientation reports polymerase direction.
- 🟡 The public schematic labels the historical sequencing arms by role. It does not replace them with a modern TruSeq sequence, because the defining paper does not establish that identity.

## Open questions

- 🔴 No single current sequencing-primer set can be asserted for every implementation called GRO-seq. A later named revision should be a separate page if its library architecture changes.
