# scGESTALT — research notes

Primary protocol: Raj and Gagnon, *Nature Protocols* (2018), [10.1038/s41596-018-0058-x](https://doi.org/10.1038/s41596-018-0058-x) ([open manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279253/)). Defining study: Raj *et al.*, *Nature Biotechnology* (2018), [10.1038/nbt.4103](https://doi.org/10.1038/nbt.4103).

## 1. What it is

Single-cell GESTALT (scGESTALT) writes cumulative CRISPR–Cas9 indels into an integrated target array at multiple developmental times, expresses the edited array as RNA, and reads it with the same cell identifier as an inDrops transcriptome.

## Sources read

- 🟢 The protocol's experimental design, target-array capture, inDrops amplification, nested lineage PCR and run-configuration sections were read.
- The source fetch command for the defining study is `python3 tools/get_sources.py --doi 10.1038/nbt.4103`; the protocol above is the primary chemistry authority for this page.

## Construct-changing steps

1. 🟢 A Tol2-integrated transgene expresses a DsRed-associated CRISPR target array. Early guide injection and later inducible Cas9/guide expression create cumulative indels.
2. 🟢 The edited cassette is transcribed and polyadenylated, allowing it to be captured alongside cellular mRNA in inDrops.
3. 🟢 Reverse-transcribed molecules receive the inDrops cell identifier, undergo second-strand synthesis and are linearly amplified by IVT.
4. 🟢 The amplified RNA is split. Fragmented material enters the standard transcriptome library; unfragmented material preserves the lineage array.
5. 🟢 Random-hexamer RT and two nested PCRs enrich the expressed lineage barcode while retaining the cell identifier and UMI. Matched Illumina multiplexing indexes are added.

## Identifier structure and read layout

- 🟢 The lineage library run is: Read 1, 250 cycles (edited barcode array); Index 1, 8 cycles (cell-barcode part 1); Index 2, 8 cycles (library index); Read 2, 14 cycles (8-base cell-barcode part 2 followed by 6-base UMI).
- 🟢 The two cell-barcode parts form one combinatorial cell identifier. The UMI identifies lineage-transcript molecules; the library index is a sample/library index.
- 🟡 Exact inDrops primer bases remain role-labelled because the accessible protocol text establishes run geometry more directly than a single assembled final-library strand.

## Boundary

The public page documents scGESTALT's expressed lineage-cassette branch. GESTALT genomic-DNA recovery and later optimized scGESTALT cassettes are related protocols, not silently merged chemistry.
