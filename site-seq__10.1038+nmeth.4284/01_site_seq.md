# SITE-Seq — research notes

Cameron et al., *Nature Methods* (2017), [doi:10.1038/nmeth.4284](https://doi.org/10.1038/nmeth.4284).

## Assay target

The sequenced material is purified human genomic DNA captured at in-vitro Cas9 cleavage sites. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a biotinylated end-capture adapter and library primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The published SITE-Seq Supplementary Protocol (March 2017, Springer `MOESM94`), read in full: reagent list and oligonucleotide table (lines 33–105), section C Cas9 digestion (lines 163–202), D–H SPRI cleanups, dA-tailing and Adapter 1 ligation (lines 203–287), G fragmentation (lines 288–317), I end repair and Adapter 2 ligation (lines 341–383), J streptavidin purification (lines 384–411), K recovery and indexing PCR (lines 412–463), L–O pooling, cleanup, QC and sequencing (lines 464–515).
- 🟢 Supplementary Table 26 (`MOESM93`, sheet `Table26_SITE-Seq Oligos`, lines 12–27), an independent spreadsheet transcription of the ten oligos AS1A/AS1B, AS2A–D, RF/RR, IF/IR and the modification key.
- The defining *Nature Methods* article full text is not in the local bundle (publisher page only; the author-hosted mirror returns HTTP 403). The chemistry below therefore rests on the two published supplements.

## Construct logic

1. 🟢 High-molecular-weight purified genomic DNA is incubated with Cas9–sgRNA at a concentration series. SITE-Seq therefore measures biochemical cleavage susceptibility in naked DNA, not cleavage inside cells.
2. 🟢 The 4 h, 37 °C cleavage reaction (RNP 0.25–1024 nM, guide:Cas9 ≥ 3:1) is terminated with a Proteinase K / RNase A mix (37 °C 20 min, 55 °C 20 min) and cleaned with 1× SPRIselect. Cas9-treated DNA is then 3′ dA-tailed with Klenow Fragment (37 °C, 30 min) and ligated to Adapter 1. Its forward strand carries a 5′ biotin and a phosphorothioate bond before the terminal T. Ligation is with NEB Quick Ligase (20 °C 30 min, then 16 °C overnight). Because this happens before fragmentation, Adapter 1 uniquely marks a Cas9-created end. The 0.5× SPRI 2 cleanup that follows is described as removing residual Adapter 1 so it cannot form adapter dimer — a physical removal, now drawn as a waste branch.
3. 🟢 Adapter-1-labelled DNA is fragmented with NEBNext dsDNA Fragmentase (37 °C, 1 h, EDTA quench) to approximately 200–800 bp. The NEBNext Ultra End Repair / dA-tailing enzyme mix (20 °C 30 min, 65 °C 30 min) then prepares the new distal ends.
4. 🟢 Adapter 2 is assembled as an equimolar mixture of forward strands containing five, six, or seven random bases and a common reverse strand (95 °C 5 min, cooled to RT over ~45 min). The source does not call these random bases a UMI or use them for molecule counting; the schematic conservatively labels them a length-diversity spacer.
5. 🟡 Adapter 2 end geometry is derived, not printed: the 15-base Adapter 2 Reverse strand (`ACTACTAATACGAC*T`) pairs with the first 14 bases of each Forward strand, leaving its terminal phosphorothioate T unpaired as a 3′-T overhang. That matches the named Blunt/TA Ligase Master Mix and the End Repair / dA-tailing module, so the distal junction is modelled as T/A rather than blunt. The protocol's own section heading says only "end repair", so the complementary dA on the fragment is inferred from the module identity.
6. 🟢 After Adapter 2 ligation, Dynabeads M-280 streptavidin enrich fragments carrying biotinylated Adapter 1 (washes in 1× B&W, then Tris). Recovery PCR is run **on the beads** with the published RF/RR pair and Q5 Hot Start 2× master mix (12 cycles, 61 °C anneal); 30 µL of supernatant is taken off the pellet, diluted 1:50, and a second 12-cycle PCR at 60 °C installs P5/P7 plus eight variable bases on each side. The pool receives a 0.7× SPRI cleanup.
7. 🟢 The supplementary protocol recommends 2–3 million 150-bp single-end reads per sample (section O). Read 1 enters genomic DNA from the Adapter 1/Cas9-cleavage end.

## The published Index Reverse primer shifts the i7 index

🟢 Index Reverse is printed identically in both supplements as `CAAGCAGAAGACGGCATACGAGATTANNNNNNNNGACTGGAGTTCAGACGTGTGCTC`. Its 3′ 19 bases (`GGAGTTCAGACGTGTGCTC`) prime on the Adapter 2 arm, so everything 5′ of them is a tail that is copied into the library: the finished top strand runs `…GAACTCCAGTC`-[i7 index]-`TA`-P7′, where a stock TruSeq construct would read `…GAACTCCAGTCAC`-[i7 index]-P7′.

🟡 Consequence: the i7 index sits two bases closer to the Read 2 arm, and the stock TruSeq Index 1 sequencing primer (…`AGTCAC`) has no exact site on this library. The model now builds the distal end from the published primer rather than from the canonical `INDEX1_PRIMER` constant, declares the Index 1 primer as mismatching with that reason, and says so on the schematic. 🔴 No custom index-read primer and no acquired index cycle counts are published, so how the i7 index is actually read is unverified; the i5 index is read off the P5 arm in the ordinary way.

## Exact oligos and modifications

The public page reproduces all ten oligos from the supplementary protocol. `/5Biosg/` is the 5′ biotin modification, `/5Phos/` is a 5′ phosphate, and `*` is the published phosphorothioate linkage. The displayed N7 construct is one of three Adapter 2 branches; N5 and N6 differ only in spacer length.

## Open questions

- 🔴 The defining article full text (paywalled) has not been read; only the two published supplements.
- 🔴 The i7 index read primer, and the acquired index cycle counts, are not published (above).
- 🟡 The end-repair step's dA-tailing activity is inferred from the NEBNext Ultra End Repair / dA-tailing module name, not stated.
- 🔴 The sgRNA sequences are study-specific and stay unassigned in the model.

## Interpretation boundary

SITE-Seq candidate sites are supported by a pileup of reads beginning at the same in-vitro cleavage coordinate. Biological editing still requires a cell-based validation assay; the library chemistry alone does not establish that a candidate site is cut in vivo.

## Source-owned validation migration, 2026-10-18

Replaced shared dsb_capture PCR finishing with local published index-primer inputs, Q5 polymerase and source run declaration. Added paired Adapter1 and all three Adapter2 input geometries, explicit dA-tailed substrate and an adapter-derived N7 pre-PCR template. RecoveryPCR and indexingPCR now locate exact sites on their actual predecessors; source primer names are preserved, although Recovery Reverse is the construct-left extending primer in this displayed orientation. Index Reverse retains the published noncanonical i7 junction and mismatching-stock-primer notice. Only150-base single-end Read1 acquisition is declared; eight variable positions in the index primers do not establish acquired index cycles.

Source scope: Published supplementary protocol MOESM94 PDF oligo table lines64–103, sectionK recovery/indexing PCR lines412–458 and sectionO sequencing recommendation around515.

Deferred: N5/N6 molecular product branches remain represented by the N7 geometry example rather than separate graph lineages; Exact index-read primer identity/acquisition and study sample-index members remain unresolved; Full cleavage/ligation/copying product relations and final pooled-library cleanup coverage remain unresolved. This bounded check establishes neither complete chemistry validation nor an accepted missing-information claim.
