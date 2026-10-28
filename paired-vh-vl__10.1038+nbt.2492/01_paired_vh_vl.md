# Paired VH:VL — physical linkage of antibody-chain amplicons

DeKosky et al. (2013), *High-throughput sequencing of the paired human immunoglobulin heavy and light chain repertoire*, [10.1038/nbt.2492](https://doi.org/10.1038/nbt.2492).

## Bounded sources inspected

Two originals were read: the published-version text available as
`nbt.2492_PMC3910347.xml.txt`, Results overview (paragraph 71), and Online Methods,
Analysis of memory B-cell VH:VL pairings in response to seasonal influenza
vaccination (paragraphs 103–107); and the published supplement
`nbt.2492_41587_2013_BFnbt2492_MOESM35_ESM.pdf.txt`, Supplementary Figure 1 and
Tables 5–7; the rendered original Supplementary Figure 1 (PDF page 2, printed
page 10) was inspected as well. The available article XML is a **final-version author manuscript**;
this provenance is not equated to a preprint. Its Methods was compared with the
published supplement for the reported linkage architecture. The source-bundle
oligo scraper was run. Large supplementary sequencing-result tables are data,
not experimental primer authority; one >10 MB result table was automatically
skipped by the scraper, not judged irrelevant by an absence audit.

## Source-defined native pairing

🟢 Cells settle into 125-pL microwells. Poly(dT) magnetic beads capture mRNA after
lysis under a dialysis membrane. Beads are recovered and washed before they are
emulsified with Quanta OneStep Fast RT-PCR mix and the Table 5 primer pools.
RT is followed by overlap-extension PCR, linking heavy and light chains from the
captured material. The emulsion is broken, cDNA recovered on a silica column,
then amplified with Table 6 nested primers and gel-purified.

🟢 The approximately 850-bp product runs from N-terminal CH1 through VH, the
linker, VL and the N-terminal Cκ/Cλ region. Paired-end MiSeq 2×250 sequencing
interrogates the terminal CDR-H3/CDR-L3 portions. Separate VH/VL amplifications
(Table 7) are an optional branch for upstream variable-region information; they
must not be substituted for the linked-repertoire branch. This is physical
chain linkage, rather than separate chains associated only by a cell barcode.

## Complementary overlap tails

🟢 Table 5 heavy-variable primers start with `TATTCCCATGGCGCGCC`; light-variable
primers start with `GGCGCGCCATGGGAATA`. Their reverse-complement relationship is
computed in the local model. Supplementary Figure 1 describes priming of first
strands, incorporation of complementary tails, denaturation, annealing and
extension to generate the linked product. No biological V(D)J sequence is
invented to demonstrate a complete source-derived receptor.

🟢 Representative Table 5 oligos (written 5′→3′):

| Primer | Sequence |
|---|---|
| CHrev-AHX89 | CGCAGTAGCGGTAAACGGC |
| CLrev-BRH06 | GCGGATAACAATTTCACACAGG |
| hVH4-DP63-fwd-OE | TATTCCCATGGCGCGCCCAGGTGCAGCTACAGCAGTGGG |
| hVK5-fwd-OE | GGCGCGCCATGGGAATAGCCGAAACGACACTCACGCAGTCTC |

These representatives do not replace the full degenerate primer pools. Many
other Table 5 sequences contain IUPAC positions that must remain degenerate if
modelled. A complete primer-table transcription and family-site coverage audit
remain deferred; inspected sequences are not an absence claim.

## Current model scope and deferred targets

The local representation models bead capture, emulsion RT/linkage, recovery,
nested PCR and gel selection. A primer-tail overlap can be checked without
claiming transcript-dependent amplification specificity or correct single-cell
pairing. Exact chain templates, nested primer binding and sequencing adapter
construction remain unassigned in this capped pass. A source-limited handoff to
sequencing is explicit; no generic ligation/tagmentation kit is invented.

Deferred: full rendered Tables 5–7 comparison; remaining
article/supplement sections; underlying cited repertoire primer methods
(references 3 and 25) and exact inheritance scope; complete sequencing-library
adapter/primer procedure; optional separate-chain and expression-validation
branches; independent review. PMC manual-download rows are not evidence that a
published counterpart is absent: the published MOESM35 PDF was obtained directly.
No missing-information claim is accepted.

The rendered Supplementary Figure 1c explicitly shows a complementary-strand
**by-product**. The graph retains this material class with an unresolved terminal
fate; it does not silently declare it discarded by the gel or incapable of
amplification. This limitation is displayed by the graph renderer. The local
linkage reference checks the two primer tails and retains one overlap, while
chain bodies/lengths are schematic placeholders, not a source-defined receptor.
