# DISCOVER-Seq — research notes

Wienert et al., *Science* (2019),
[doi:10.1126/science.aav9023](https://doi.org/10.1126/science.aav9023).

## Assay target

The sequenced material is human chromatin fragments enriched by MRE11 occupancy at nuclease-induced double-strand breaks. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a stock ChIP-seq library kit). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining Science author manuscript and associated description of MRE11-bound
  read endpoints.
- 🟢 Wienert et al., *Nature Protocols* (2020),
  [doi:10.1038/s41596-020-0309-5](https://doi.org/10.1038/s41596-020-0309-5),
  the authoritative stepwise DISCOVER-Seq procedure.
- 🟢 Illumina, *TruSeq ChIP Sample Preparation Guide*, part 15023092 Rev. B, the manual
  explicitly named by the detailed protocol.

## Construct logic

1. 🟢 CRISPR nuclease is delivered to living cells or tissue and MRE11 is allowed to bind
   newly made double-strand breaks. The appropriate capture time depends on delivery;
   the detailed protocol suggests roughly 8–16 h for RNP editing.
2. 🟢 Cells are formaldehyde-crosslinked, lysed and sonicated to 150–500-bp chromatin
   fragments (Nature Protocols workflow and lines 1233–1237, 1430–1448 of the extracted
   article). This is ordinary mechanical fragmentation around a biologically created end.
3. 🟢 MRE11-bound chromatin is immunoprecipitated with an anti-MRE11 antibody prebound to
   Protein A magnetic beads. Reverse crosslinking, RNase A and proteinase K release DNA.
   ChIP-qPCR at the intended target is a checkpoint before library preparation.
4. 🟢 The detailed procedure delegates the remaining chemistry to the Illumina TruSeq
   ChIP library kit (Nature Protocols lines 799–845). The kit manual specifies end repair,
   3′ adenylation, indexed-adapter ligation, size selection and PCR enrichment.
5. 🟢 The historical kit is single-indexed. Its guide says to record six index bases;
   AR001 is `ATCACG`, with the parenthesized seventh base explicitly excluded from the
   sample-sheet index (TruSeq guide lines 2109–2186).
6. 🟢 Libraries are sequenced as 100-bp paired-end reads at about 30–40 million read pairs
   per ChIP sample (Nature Protocols lines 843–849). The BLENDER analysis uses the pileup
   and orientation of fragment endpoints around a cut rather than an assay-specific UMI.

## Why the end is informative

One end of an enriched fragment can be the cellular nuclease break occupied by MRE11;
the other is normally a random sonication boundary. A real cut therefore produces many
MRE11-enriched fragments whose reads start or end at a common genomic coordinate, while
their distal endpoints vary. The paired-end construct on the public page preserves that
asymmetry without pretending the genomic sequence is fixed.

## Commercial-sequence boundary

The six-base adapter indexes and the final TruSeq sequencing architecture are public.
The discontinued ChIP kit guide does not provide every supplied reagent as a complete
ordering oligo, so the page does not invent such sequences or call the library dual-indexed.
