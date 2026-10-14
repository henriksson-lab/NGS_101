# SITE-Seq — research notes

Cameron et al., *Nature Methods* (2017), [doi:10.1038/nmeth.4284](https://doi.org/10.1038/nmeth.4284).

## Assay target

The sequenced material is purified human genomic DNA captured at in-vitro Cas9 cleavage sites. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a biotinylated end-capture adapter and library primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The published SITE-Seq Supplementary Protocol (March 2017), including its oligonucleotide table and reaction instructions.
- 🟢 Supplementary Table 26, an independent spreadsheet transcription of the SITE-Seq oligos and modification key.

## Construct logic

1. 🟢 High-molecular-weight purified genomic DNA is incubated with Cas9–sgRNA at a concentration series. SITE-Seq therefore measures biochemical cleavage susceptibility in naked DNA, not cleavage inside cells.
2. 🟢 Cas9-treated DNA is 3′ dA-tailed and ligated to Adapter 1. Its forward strand carries a 5′ biotin and a phosphorothioate bond before the terminal T. Because this happens before fragmentation, Adapter 1 uniquely marks a Cas9-created end.
3. 🟢 Adapter-1-labelled DNA is fragmented with NEBNext dsDNA Fragmentase to approximately 200–800 bp. End repair then prepares the new distal ends.
4. 🟢 Adapter 2 is assembled as an equimolar mixture of forward strands containing five, six, or seven random bases and a common reverse strand. The source does not call these random bases a UMI or use them for molecule counting; the schematic conservatively labels them a length-diversity spacer.
5. 🟢 After Adapter 2 ligation, streptavidin beads enrich fragments carrying biotinylated Adapter 1. Recovery PCR uses the published RF/RR pair, and a second PCR installs P5/P7 plus eight variable bases on each side.
6. 🟢 The supplementary protocol recommends 2–3 million 150-bp single-end reads per sample. Read 1 enters genomic DNA from the Adapter 1/Cas9-cleavage end. The two eight-base sample indexes are derived from the final construct and primer directions.

## Exact oligos and modifications

The public page reproduces all ten oligos from the supplementary protocol. `/5Biosg/` is the 5′ biotin modification, `/5Phos/` is a 5′ phosphate, and `*` is the published phosphorothioate linkage. The displayed N7 construct is one of three Adapter 2 branches; N5 and N6 differ only in spacer length.

## Interpretation boundary

SITE-Seq candidate sites are supported by a pileup of reads beginning at the same in-vitro cleavage coordinate. Biological editing still requires a cell-based validation assay; the library chemistry alone does not establish that a candidate site is cut in vivo.
