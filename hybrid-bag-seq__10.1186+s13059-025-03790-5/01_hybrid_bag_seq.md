# Hybrid BAG-seq — joint single-nucleus DNA and RNA sequencing

## Sources read

- 🟢 Li *et al.*, “Joint single-cell DNA and RNA sequencing with hybrid BAG-seq,” *Genome Biology* (2025), DOI `10.1186/s13059-025-03790-5`, PMCID `PMC12476631` (defining paper and methods).
- 🟢 Li *et al.*, “Single-cell DNA sequencing based on hydrogel whole-genome amplification,” *Genome Research* (2020), DOI `10.1101/gr.253047.119`, PMCID `PMC6961581` (original BAG platform and printed capture oligos).
- 🟢 Stepwise hybrid protocol deposited as protocols.io DOI `10.17504/protocols.io.rm7vzqr92vx1/v1`.
- 🟢 The protocols.io PDF's **Primers** section and its attached workbook print
  the complete three-round barcode families, indexing primers and custom Read 1
  primer. They were retrieved through `tools/get_sources.py` and read as text.

## Molecular interpretation

- 🟢 A nucleus is co-encapsulated with acrylamide monomer, a 5′-Acrydite random-TG DNA primer and a 5′-Acrydite poly(dT) RNA primer. Both primers become covalently embedded in the gel.
- 🟢 Hybrid BAG formation uses 85 °C for 5 min rather than the DNA-only protocol's 95 °C for 12 min, avoiding severe RNA destruction.
- 🟢 Reverse transcriptase and DNA polymerase make gel-bound RNA- and genomic-DNA-derived copies. A TSO gives cDNA a double-stranded NlaIII site; Klenow extends captured genomic DNA.
- 🟢 NlaIII cleavage exposes a sticky end. Round 1 ligation-extension adds barcode 1 plus varietal tag 1; round 2 hybridization-extension adds barcode 2 plus varietal tag 2 and a common adapter; split PCR adds barcode 3.
- 🟢 The complete BAG barcode is 17 nt (5 + 6 + 6). The 12-nt varietal tag comprises 4 + 4 random bases plus four genomic bases two bases from the CATG cut site.
- 🟢 Pooled amplicons are tagmented and sequenced 2 × 150 bp on NovaSeq 6000. Read 1 begins with an 80-base tag/adapter region; Read 2 starts from a random tagmentation end.
- 🟢 The author-deposited protocol supplies exact barcode-family and primer
  sequences, though the current schematic has not yet assembled them into a
  final base-level construct.

## Verbatim capture oligos

- 🟢 DNA: `/5ACryd//iSp18/TGTGTTGGGTGTGTTTGGKKKKKKKGKKKKKKKKNN`
- 🟢 RNA: `/5ACryd//iSp18/AAGCAGTGGTATCAACGCAGAGTNNWNNNSTTTTTTTTTTTTTTTTTTTTTTTTTTTT`

These are retained with modification syntax and IUPAC symbols exactly as printed by the original BAG paper.

## Hybrid-specific oligos recovered from the author protocol

The stepwise hybrid protocol's **Primers** section (PDF text lines 77–89) prints
the following oligos. The attached workbook has 96 rows for each barcoding
round; its first rows are shown as concrete representatives. Spaces separate
source columns and are not nucleotide symbols. `+` marks locked nucleotides,
`*` phosphorothioate bonds, and `N`/`B`/`K`/`D` are IUPAC or degenerate bases.

| Role | Published 5′→3′ oligo or representative |
|---|---|
| Hybrid DNA capture | `/5ACryd//iSp18/TGTGTTGGGTGTGTTTGG NBBBB KKKKKKKGKKKKKKKKNN` |
| Hybrid RNA capture | `/5ACryd//iSp18/TGTGTTGGGTGTGTTTGG DDDDDNNN T(30)` |
| Template switching | `AAGGCTCTTGACATGAAGGATAGrGrGrG` |
| Round 1, member 1 | `ACGCAGAGTCCGCTC CTAAG NNNN CA+T+G` |
| Round 2, member 1 | `GACCGACTCGCATTACCCTAT AATTTC NNNN ACGCAGAGTCCGCTC/3InvdT/` |
| Round 3, 2N member 1 | `AAGGCTCTTGACACGAAGGATAG TACCGTNN GACCGACTCGCATTACCCT*A*T` |
| Final P5-side PCR primer | `AATGATACGGCGACCACCGAGATCTACAC CGCTGTCC AAGGCTCTTGACACGAAGGAT*A*G` |
| Final P7-side N701 primer | `CAAGCAGAAGACGGCATACGAGAT TCGCCTTA GTCTCGTGGGCTCGG` |
| Custom Read 1 primer | `CGCTGTCC AAGGCTCTTGACACGAAGGATAG` |

🟢 The paper's analysis Methods additionally print the Read 1 and Read 2
trimming motifs `CCAAACACACCCAA` and `GAGCGGACTCTGCG`, respectively.
These are observed adapter-related read motifs, not separately ordered
oligos. 🟡 Exact pairing, strand orientation and sequencing-primer landing
should be computed from the published oligos in a molecular model rather
than copied from the trimming rules. 🔴 The current shared factory still
draws only roles for the final library.
