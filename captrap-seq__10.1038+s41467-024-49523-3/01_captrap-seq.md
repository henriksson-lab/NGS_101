# CapTrap-seq — research notes

> **Evidence marks:** 🟢 printed in the published source · 🟡 derived or inferred here ·
> 🔴 not established from the sources read.

## Sources read

- 🟢 Carbonell-Sala S, Perteghella T, Lagarde J, et al., including Carninci P,
  Uszczynska-Ratajczak B and Guigó R. "CapTrap-seq: a platform-agnostic and quantitative
  approach for high-fidelity full-length RNA sequencing." *Nature Communications* 15, 5278
  (2024), doi:[10.1038/s41467-024-49523-3](https://doi.org/10.1038/s41467-024-49523-3). Published
  full text read from the PMC Cloud mirror PMC11211341 (JATS XML `.txt`), Methods
  "CapTrap-seq full-length cDNA libraries" (text lines 206–238) and "ONT and PacBio
  sequencing" (lines 252–260). This is the authority for every oligo and step below.
  Earlier repository text cited this paper as "Glinos et al."; that was a repository
  attribution error, corrected here from the JATS author list.
- 🟢 Published Springer supplementary PDFs MOESM1–3 are in the source bundle; the bounded
  2026-10-17 pass did not re-read them beyond the scraper (no oligo hits).
- The author bioRxiv preprint (doi:10.1101/2023.06.16.543444, v1) was the source of the
  earlier version of this note. It is superseded by the publication. Its Table S1 printed
  the UMI8 primer one base shorter (`…CAGAGTA`); this note follows the published `…CAGAGTAC`.
- 🔴 Not read: the full protocol cited as Nature Protocol Exchange (ref. 45), the cited
  two-step linker-ligation predecessor (ref. 46), and the ONT SQK-LSK108/109 and PacBio
  SMRTbell Express Template Prep Kit 2.0 manuals.

## Oligos (published Methods, verbatim)

🟢 `[Phos]`/`P` is a phosphate, `U` a uracil residue, `V`/`N` degenerate bases; the UMI8
primer carries an 8-nt UMI. The two linker strands are printed on separate lines without
5′/3′ labels.

- OligodT 16VN: `5′ [Phos]TTTTTTTTTTTTTTTTVN 3′`.
- 5′ linker: `GTGGTAUCAACGCAGAGUACGNNNNN-P` / `P-CACCATAGTTGCGTCTCATG-P`.
- 3′ linker: `AAAAAGCAUCGCUGTCTCUTAUACACAUCUCCGAGCCCACGAGAC-P` /
  `CGTAGCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG`.
- Second primer UMI8: `5′ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNGTGGTATCAACGCAGAGTAC 3′`.
- LA-PCR primers: forward `5′ TCGTCGGCAGCGTC 3′`, reverse `5′ GTCTCGTGGGCTCGG 3′`.

## How the oligos interlock (computed)

🟡 Each printed lower linker strand is the plain complement of its upper strand, i.e. it is
printed aligned 3′→5′. Read 5′→3′, the 5′-linker lower strand is
`{{= "CACCATAGTTGCGTCTCATG"[::-1] }}` and the 3′-linker lower strand is
`{{= "CGTAGCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG"[::-1] }}`. The model's `Scene.anneal`
checks that these pair antiparallel with the upper strands (20 and 40 columns). The
resulting geometry is:

- 5′ linker: the 3′-blocked upper strand leaves a 3′ `G` + `N5` tail beyond the duplex; it
  splints the lower strand's 5′ phosphate against the single-strand cDNA 3′ end (the end
  copied from the capped RNA 5′ end).
- 3′ linker: the upper strand leaves a 5′ `AAAAA` tail; it splints the lower strand's
  3′-OH against the 5′-phosphorylated OligodT 16VN end of the cDNA. This explains why the
  RT primer is 5′-phosphorylated.
- SAP/USER then removes the U-containing upper strands and 3′ phosphates, leaving
  single-strand handles. The UMI8 primer's 3′ 20 nt (`GTGGTATCAACGCAGAGTAC`) pair with the
  ligated 5′-linker lower strand (checked by `Scene.anneal`).
- The forward LA-PCR primer is the UMI8 primer's first 14 nt; the reverse LA-PCR primer is
  the first 15 nt of the 3′-linker lower strand. `PCRGeometry.exact` recomputes both inward
  sites on the drawn LA-PCR product.

🔴 The cDNA bases opposite the `G-N5` and `A5` splint tails (for example a cap-templated
extra C) are not stated; the model does not draw them as sequence.

## Step by step (published Methods)

1. 🟢 RT: 5 µg total RNA with OligodT 16VN, 65 °C 5 min; PrimeScript II RT 42 °C 60 min;
   RNAClean XP 1.8x.
2. 🟢 Cap trapping: NaIO4 oxidation of diols (5 min on ice; Tris pH 8.5 stop), Biotin
   (Long Arm) Hydrazide 30 min 40 °C, RNase ONE 30 min 37 °C to remove single-strand RNA and
   its 3′-end biotin, capture on M-270 streptavidin beads in LiCl buffer, three washes.
3. 🟢 Release: two 95 °C 5 min releases in RNase ONE buffer; RNase H + RNase ONE 30 min
   37 °C; AMPure XP 1.8x; speed-vac to 4 µl.
4. 🟢 5′ linker ligation with Mighty Mix (linker denatured 55 °C, sample 95 °C; 4 h 30 °C or
   16 h 16 °C); two AMPure XP 1.8x purifications to remove non-incorporated linker.
5. 🟢 3′ linker ligation the same way (linker denatured 65 °C); AMPure XP 1.8x.
6. 🟢 SAP + USER 30 min 37 °C, 5 min 95 °C; AMPure XP 1.8x.
7. 🟢 Second-strand synthesis with KAPA HiFi, DMSO and the UMI8 primer (95 °C 5 min, 55 °C
   5 min, 72 °C 30 min); Exonuclease I 30 min 37 °C; AMPure XP 1.8x then 1.4x.
8. 🟢 LA-PCR with TaKaRa LA Taq in two independent reactions, 16 cycles (95 °C 15 s,
   55 °C 15 s, 65 °C 8 min); merged, AMPure XP 1x.
9. 🟢 Platform libraries: ONT from 500–1000 ng cDNA "strictly following" SQK-LSK108 and
   SQK-LSK109 (MinION, R9.4); PacBio Sequel I/II from 500 ng with SMRTbell Express Template
   Prep Kit 2.0. No size selection before loading.

## Final library and reads

🟡 The LA-PCR product is, 5′→3′ on the mRNA-sense strand: forward-primer site, 19-nt
spacer, UMI (8 nt), 5′-linker-complementary 20 nt, cDNA, anchor complement, poly(A) copy,
then the 3′-linker-derived 40 nt. 🔴 The ONT kit end-prep, adapter bases and adapter
positions, and the PacBio SMRTbell junctions, are not disclosed in the article; the page
shows the ONT adapter as an inferred role token and does not draw the PacBio branch.

## Open questions

- Protocol Exchange full protocol and the ref. 46 linker-ligation method: linker annealing
  conditions and any stated end geometry.
- ONT SQK-LSK108/109 manual revisions used (end-prep/dA-tailing, adapter structure).
- Whether the UMI is used in the published analysis (not checked in this pass).
