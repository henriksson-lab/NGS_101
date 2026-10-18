# BLESS — research notes

Crosetto et al., *Nature Methods* (2013), [doi:10.1038/nmeth.2408](https://doi.org/10.1038/nmeth.2408), PMCID: PMC3651036.

## Sources read

- 🟢 Author-manuscript full text from PMC.
- 🟢 Supplementary Information as published by *Nature Methods* (Springer `41592_2013_BFnmeth2408_MOESM325_ESM.pdf`), used in preference to the PMC author-manuscript copy (`NIHMS449608-supplement-1.pdf`): Supplementary Table 1 (sequencing platforms/modes, text lines 455–690) and Supplementary Table 5 (page 12, text lines 1126–1199; rendered page inspected for the colour-coded biotin, barcodes and XhoI/I-SceI sites).
- 🟢 Online Methods paragraphs read in round 2026-10-17-breadth-first-07 (PMC JATS `nmeth.2408_PMC3651036.xml.txt`): linkers obtained from Sigma and annealed in 1× T4 ligase buffer (line 125); blunting, in situ ligation, gDNA extraction and HaeIII fragmentation (line 131); streptavidin capture, distal ligation and I-SceI digestion (line 133); P1/P3 PCR and XhoI digestion (line 135); Illumina library preparation and sequencing (line 139); Illumina read filtering (line 153).

## Molecular interpretation

1. 🟢 Cells are fixed and gently lysed while nuclei remain intact. Native DSB ends are blunted and 5′ phosphorylated (Quick Blunting kit, 45 min RT).
2. 🟢 The 91-nt proximal hairpin L1 carries a biotinylated base and is ligated to DSBs in situ with T4 ligase (18–20 h, 16 °C). 🔴 Whether the annealed linker itself was 5′-phosphorylated is not stated; Supplementary Table 5 prints no phosphate. 🟡 After XhoI cutting its `TCGAGGTAGTA` remnant identifies a proximal end (XhoI C^TCGAG site + red GTAGTA barcode in Table 5).
3. 🟢 gDNA is extracted (proteinase K, isopropanol–ethanol) and fragmented with **HaeIII** (6 U per million cells, 18–20 h, 16 °C), then captured on Dynabeads MyOne C1 streptavidin. This digestion creates the free distal end; it is not a second native break.
4. 🟢 The annealed distal linker and T4 ligase are added directly to the captured beads (16–18 h, 16 °C); 🟡 no separate distal blunting step is described, and HaeIII (GG^CC) cuts blunt. The later L3 design leaves `TCGAGACGACG` as the distal tag. *(Model correction, round 2026-10-17-breadth-first-07: an earlier model drew a distal blunting step and a generic fragmentation condition.)*
5. 🟢 I-SceI (4 h, 37 °C) releases the captured construct into the supernatant; P1/P3 (Illumina) amplify it with Phusion, 18 cycles, Ta 55 °C; the product is gel-purified, cut with XhoI and gel-purified again, leaving the two short end tags. *(Model correction: the earlier graph cut with XhoI before any PCR and replaced the P1/P3 PCR with a generic "amplify and index" step.)*
6. 🟢 The paper's Methods, "Next-generation sequencing," specifies TruSeq DNA Sample Preparation Kit v2 (Illumina), without fragmentation or library size selection; the model represents this as a warned kit boundary. Reads containing the proximal tag report the native DSB; reads containing the distal tag start at the later fragmentation end.
7. 🟢 Supplementary Table 1: Illumina GAII single-end for most samples, HiSeq 2000 paired-end for N1/N2 (and gDNA controls). 🔴 Read lengths and any index read are not stated; Methods retain reads ≥ 34 nt after quality trimming. *(Model correction: the earlier page declared unsourced 75 + 6 + 75 cycles.)*
8. 🟡 The public diagram draws the final sequencing-kit shell in the familiar TruSeq orientation solely to locate the standard sequencing primers. The paper does not print the complete kit adapter oligos, so the outer shell is visibly inferred; the BLESS insert and tags are not.

## Exact oligos

Supplementary Table 5 gives L1, L3, P1 and P3. Supplementary Fig. 1 and the analysis section independently identify the retained 11-base proximal and distal tags. The published Table 5 marks the biotinylated base only by cyan highlighting (the 46th base of L1, inside the seven-T loop); it is rendered `/iBiodT/` on the page, the internal biotin-dT form being 🟡 an interpretation of that highlight. No 5′-phosphate is printed for L1 or L3, so the page displays none.

The first published experiments also used distal linker L2. The later L3 carries a separate distal barcode and is the form represented on the schematic because it is the specificity-improved BLESS construct. Supplementary Table 5 also prints P4/P5 for a Roche 454 branch; those are not omitted Illumina PCR primers. For the illustrated Illumina branch, L1/L3 and P1/P3 are complete through the BLESS insert. The downstream kit is explicitly TruSeq DNA Sample Preparation Kit v2; the paper does not print its adapter oligos or specify the selected index member.

## Open questions

- 🔴 TruSeq DNA Sample Preparation Kit v2 original (end repair/A-tailing on XhoI-cut ends, adapter oligos, PCR, index member) not inspected.
- 🔴 Linker 5′ phosphorylation and the exact biotin chemistry (internal biotin-dT vs. another linkage) are not stated.
- 🔴 Read lengths and index reads for the GAII/HiSeq runs.
- The L2 initial distal linker and the Roche 454 P4/P5 branch are not drawn.

## Source-owned validation migration, 2026-10-18

Protocol now owns its conditional TruSeq v2 endpoint instead of inheriting a six-base break-capture shell. The selected kit index member and molecular length remain unknown; exact proximal/distal tags retain their structured identities. Graph consumes the full source oligos, with their source-owned HairpinTopology retained as chemical state; public folded views remain. The sequencing plan explicitly selects the representative N1/N2 paired-end HiSeq2000 branch and distinguishes GAII single-end runs; no read/index cycle counts are assigned. TruSeq library preparation remains an explicitly unmodelled source-limited transition.

Source scope: Published article PMC3651036 Methods lines128–139 (proximal ligation, extraction/HaeIII, streptavidin/distal/I-SceI, P1/P3 PCR/XhoI and TruSeq v2 library preparation); published supplementary PDF MOESM325 Table1 lines594–626 and Table5 lines1126–1174.

Deferred: TruSeq DNA Sample Preparation Kit v2 original, exact supplied adapter/primer geometry and selected index remain deferred; I-SceI-cut remnants and complete P1/P3 amplicon geometry remain unmodelled; Acquired cycle counts and other source single-end/paired-end run variants remain unassigned. This bounded check establishes neither complete chemistry validation nor an accepted missing-information claim.
