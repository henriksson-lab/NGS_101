# BLISS — research notes

Defining source: Yan *et al.*, [10.1038/ncomms15058](https://doi.org/10.1038/ncomms15058).

## Evidence used

- 🟢 Fixed cells or tissue are immobilized, DSB ends are blunted, and a double-stranded adapter is ligated directly to breaks in situ.
- 🟢 The adapter contains a T7 promoter, RA5 sequence, sample barcode and an 8–12-nt UMI.
- 🟢 T7 in-vitro transcription linearly amplifies break-adjacent sequence before completion of the Illumina RNA-style library.
- 🟢 The authors' [Protocol Exchange v2](https://assets-eu.researchsquare.com/files/nprot-5597/v2/17ba35a1-4de4-4f0a-8953-6004ecad1900.pdf?c=1631828212), pp. 3 and 10, specifies custom-synthesized TruSeq Small RNA RA3 ligation to IVT RNA, RTP reverse transcription, then RP1 and a chosen RPI index primer for PCR. The [Illumina sequence table](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm) prints RA3 `TGGAATTCTCGGGTGCCAAGG`, RTP `GCCTTGGCACCCGAGAATTCCA`, and RP1 `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`, plus the RPI family. A sample-specific RPI member remains unassigned within the inspected source scope.
- 🟢 Supplementary Data 1 prints both strands of eleven classic BLISS adapters
  and ten longer-UMI adapter variants. The page uses classic adapter A1 as its
  representative exact duplex.
- 🟢 Supplementary Fig. 1a marks the sample-barcode end as the blunt ligation
  site and shows a 5′ phosphate on its bottom strand. The opposite bottom-strand end has
  a two-base 3′ overhang that prevents head-to-tail adapter concatemers.
- 🟡 RA3/RTP/RP1/RPI are a distinct library-finishing family; Supplementary Data 1
  lists the break adapters. The model shows RA3 ligation before RT. The method says
  these oligos were custom synthesized from TruSeq Small RNA sequences, but custom RA3 end modifications and the RPI index chosen for a sample remain
  unassigned within the inspected scope.

## Bounded source pass 2026-10-17 (round breadth-first-07)

Read: article PMC5437291 JATS lines 86 (workflow), 132 (adapter preparation), 136, 148
(prefix filtering), 154/158 (R1 prefix); Supplementary Information Fig. 1a caption
(lines 1–18); Supplementary Data 1 rows A1–A11; Protocol Exchange PDF lines 160–290
(Day 1–4 procedure). The Protocol Exchange file on disk is the Research Square
revision v2 (doi:10.21203/rs.2.1448/v2, posted 13 Nov 2019). The article cites
protex.2017.018; that original revision was not inspected, so step details below are
attributed to v2.

- 🟢 Model correction: the A1 oligos previously carried one extra G at the T7/RA5
  junction (`…CACTATAGGGGG TTCAG…`). Supplementary Data 1 prints
  `…CACTATAGGGGTTCAGAG…` (75-nt top, 77-nt bottom). This was a repository
  transcription mistake, not a source error. The self-test now compares A1 with the
  workbook.
- 🟢 Methods line 132: oligos were standard desalted IDT oligos, and "sense oligos" were
  phosphorylated with T4 PNK before annealing (95 °C 5 min, cooled to 25 °C over 45 min).
  Fig. 1a puts the 5′ phosphate on the bottom (ligating) strand and says it is added
  either during synthesis or with T4 PNK. The classic Supplementary Data 1 rows print no
  /5Phos/, so the drawn `/5Phos/` is the post-kinase state. Which strand "sense" means is
  not stated.
- 🟢 Protocol Exchange v2 steps now in the graph: LB1/LB2 lysis in situ, Quick Blunting
  (1 h, RT), T4 DNA ligase 16–18 h at 16 °C, 3–5 high-salt washes that remove unligated
  adapters, proteinase K extraction (55 °C overnight), column or phenol–chloroform
  cleanup, sonication to 300–500 bp, AMPure 0.8× concentration (SpeedVac alternative),
  MEGAscript T7 IVT 14 h at 37 °C, DNase I, RNAClean XP 1.8×, RA3 ligation (70 °C 2 min,
  then T4 Rnl2tr 1 h at 28 °C), RTP/SuperScript III RT (1 h at 50 °C), RP1 + chosen RPI
  PCR with NEBNext HF 2X (10–16 cycles), then AMPure 0.8×. The CRISPR-BLISS A-tailing
  branch is not drawn.
- 🟢 Published discrepancy, shown on the page: Supplementary Data 1 adapters A3 and A6
  print the same barcode (`TGATGATC` top, `GATCATCA` bottom). The model draws A1 only.
- 🟢 Read layout: R1 begins with the 8-nt UMI plus the 8-nt barcode prefix (Methods
  lines 148, 158). Most libraries were sequenced single-end (line 88, Supplementary Data 2).

Open: the cited protex.2017.018 original revision; Protocol Exchange Supplementary
Table 2; RA3 custom end modifications; the RPI member; exact RA5-side read-primer bases
and cycle counts; an exact final-library sequence (left as roles).

## Protocol boundary

The page covers original BLISS. Its suspension implementation, sBLISS, is a separately named protocol and is not silently merged here.

## Source-owned validation migration (2026-10-18)

The final construct and graph are now declared locally rather than selected by
`protocol_specs.get`. The explicit representative is classic adapter A1: its UMI
and sample barcode each occupy eight positions and carry structured identifier
features. The former generic ten-position payload defaults have been removed.
RP1/RPI-installed arms, selected RPI, exact RA5-side sequencing-primer geometry and
acquired cycles remain unassigned; no generic TruSeq endpoint is substituted.
The graph uses actual first-strand cDNA as the PCR precursor. IVT RNA and retained
genomic DNA occupy independent molecular rows; DNase changes the DNA population,
and RNA cleanup then yields the RNA-only material. cDNA contains the complement
of the A1 barcode rather than its forward orientation.

This bounded pass re-read three originals: published PMC5437291 JATS workflow,
adapter preparation and read-prefix Methods (lines 86–88, 130–136, 148–158);
published Supplementary Data 1 classic A1–A11 rows; and the explicitly 2019 v2
Protocol Exchange PDF Day 1–4 procedure (lines 160–290). Bundle primer-scraper hits
were checked against the workbook. The original article-cited protocol revision,
its separately cached attachments, custom RA3 end modifications, selected RPI and
exact acquired run remain deferred. The 2019 v2 operations are version-qualified;
complete equivalence with the article-cited original is not asserted. Earlier
Fig.1a-based adapter-end declarations are retained from the prior review, without
a new figure reading. This is not a full source audit or a publication-absence claim.
