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
- 🟢 Pooled amplicons are tagmented, then amplified with P5-Adapter3 and one Nextera_N70x indexed primer (author protocol step 38). The libraries are sequenced 2 × 150 bp on NovaSeq 6000. Read 1 begins with an 80-base tag/adapter region; Read 2 starts from a random tagmentation end.
- 🟢 The author-deposited protocol supplies exact barcode-family and primer
  sequences. Since the 2026-10-17 bounded pass the schematic assembles the
  tag-bearing end from their common arms; barcodes, varietal tags and the
  round-3 stagger remain placeholders and the Tn5/N7/P7 end remains a role.
- 🟢 The original BAG predecessor's linked Supplemental Methods also print
  its earlier Acrydite capture primer and adapter arms. Those differ from the
  hybrid-specific primers and do not fill any hybrid source gap.

## Predecessor BAG capture oligos

- 🟢 DNA: `/5ACryd//iSp18/TGTGTTGGGTGTGTTTGGKKKKKKKGKKKKKKKKNN`
- 🟢 RNA: `/5ACryd//iSp18/AAGCAGTGGTATCAACGCAGAGTNNWNNNSTTTTTTTTTTTTTTTTTTTTTTTTTTTT`

These are retained with modification syntax and IUPAC symbols exactly as printed by the original BAG paper. They describe its DNA-only and RNA-only branches. The hybrid author protocol prints different capture primers below, which are the ones displayed on the hybrid schematic.

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

🟢 In the author protocol, round 1 adds Adapter1 using Klenow exo− plus Quick Ligase (steps 25–26). Round 2 hybridizes Adapter2 and extends it with Bst 2.0 DNA polymerase (steps 31–32); it is not a second ligation. The public reaction graph follows these separate operations.

🟢 The paper's analysis Methods additionally print the Read 1 and Read 2
trimming motifs `CCAAACACACCCAA` and `GAGCGGACTCTGCG`, respectively.
These are observed adapter-related read motifs, not separately ordered
oligos. 🟡 Exact pairing, strand orientation and sequencing-primer landing
should be computed from the published oligos in a molecular model rather
than copied from the trimming rules. 🔴 The local protocol model still draws only roles for the final library.

## Bench operations added in the 2026-10-17 bounded pass

🟢 From the protocols.io v1 PDF (text lines 340–840, steps 14–41):

- RT and DNA extension in one tube (step 15): Maxima H− RT 15 min at room temperature, TSO added for 15 min more, 42 °C for 1 h; then Klenow (M0210M) at 42 °C for 30 min.
- Exonuclease I plus RNase H, 37 °C for 1.5 h, then STOP-buffer washes (steps 17–19). The protocol does not state which molecules these enzymes remove; the graph shows the step without assigning targets.
- Two sequential NlaIII reactions (steps 21–22).
- Round 1 (steps 25–26): Adapter1-BC1-vt1-X hybridized in Quick Ligase buffer, then Klenow exo− plus Quick Ligase.
- NaOH/Brij-35 denaturation (three washes) and Tris/NaCl/EDTA/Tween neutralization before round 2 (steps 29–30). 🟡 Which strands leave is inferred from the covalent Acrydite anchor.
- Round 2 (steps 31–32): Adapter2-BC2-vt2-Adapter1-Y hybridized at 75 °C then 57 °C, then Bst 2.0 extension at 57 °C for 40 min.
- 150 µm strainer removes fused BAGs (step 34). Round 3 (step 35): TG primer `TGTGTTGGGTGTGTTT*G*G` plus Adapter3-BC3-Adapter2-Z with KAPA HiFi, 10 cycles.
- Pooled 0.8× AMPure XP cleanups, twice (step 36); about 1.5 ng tagmented with Amplicon Tagment enzyme at 55 °C for 5 min (step 38). Then 10 cycles with P5-Adapter3, Nextera_N70x and Nextera PCR mix, 0.8× AMPure (step 39). The run uses the custom Read1-Adapter3 primer (step 41).

## Read 1 order (model correction)

🟡 Derived from oligo polarity: Adapter3's 3′ arm equals Adapter2's 5′ arm, and Adapter2's 3′ arm (`ACGCAGAGTCCGCTC`) equals Adapter1's 5′ arm. Read 1 from Read1-Adapter3 therefore reads BAG barcode 3, the 2N/3N stagger, Adapter2, barcode 2 and varietal tag 2, Adapter1, barcode 1 and varietal tag 1, then CATG. The article (PMC12476631 Step 2 of data processing) places the four template varietal bases two bases *beyond* CATG. The earlier role-only model put those bases before CATG and omitted both common arms. This was a repository modelling error, now corrected. The Read 1 and Read 2 trimming motifs are consistent with this order: they are the reverse complements of UP1 and Adapter1.

## Open questions

- 🟡 For cDNA, a cut at the TSO-encoded CATG would place TSO bases (`AAGGATAG` then the rG-derived Gs) after CATG, not transcript bases. The article's description of the varietal bases as "genomic" was not reconciled with this inside the bounded pass. This is a sequence comparison, not a source statement.
- The Nextera Tn5 adapter and Read 2/index primers are not printed in the hybrid sources, and the index-read cycle count was not found in the passages inspected.
