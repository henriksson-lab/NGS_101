# scSNV-seq — source notes

Defining paper: Cooper et al., *Genome Biology* (2024), DOI
[10.1186/s13059-024-03169-y](https://doi.org/10.1186/s13059-024-03169-y).

🟢 The vector carries a 12-nt transcribed PuroR barcode and a 6-nt iBAR in the guide.
The workflow combines targeted single-cell DNA genotyping with 10x 5′ transcriptomes;
guide/iBAR and PuroR libraries use spiked-in RT primers and separate nested enrichments.

## Historical oligo audit (earlier scope, not a fresh batch06 source pass)

🟢 Sources checked: the defining PMC article and JATS, associated bioRxiv manuscript,
all five supplementary files, and especially Supplementary Table S2. That table publishes
the vector-construction oligos, the complete JAK1 Tapestri amplicon panel, both custom 10x
RT spike-ins, and the shared indexed nested primer.

🟢 The guide spike-in is oJR160,
`AAGCAGTGGTATCAACGCAGAGTACCAAGTTGATAACGGACTAGCC`; the PuroR spike-in is
`AAGCAGTGGTATCAACGCAGAGTACggcttgtactcggtcatggt`. Their 5′ SMART handle is
followed by a scaffold- or PuroR-specific binding region.

🟢 `Kosuke_Guide_i5longFwd` and `Puro_i5longFwd` are the same published primer:
`AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT`,
where the ten `N` bases are the library index.

🟡 The paper does not print the paired commercial 10x 5′ primer/bead construct or the
complete Mission Bio Tapestri sequencing shell. The exact custom oligos and target-panel
primers are shown separately. Historical vendor guides and the authors' linked code
still need a second source check before those commercial endpoints can be called absent.

## 2026-10-17 bounded source check

🟢 Read: published PMC10789043 JATS Methods "gRNA library cloning…" (line 115) and
"Base editing screens" (line 119), and Additional file 3: Table S2 (MOESM3 sheet
`Table_S2`). The chemistry paragraph is line 119: 10x 5′HT kit on Chromium X, a
scaffold-specific RT primer spiked in before loading, a nested PCR from the small SPRI
fraction after cDNA amplification for the guide + iBAR library, the analogous PuroR
spike-in and nested PCR, Tapestri (Mission Bio) single-cell genotyping with a custom
amplicon panel spanning JAK1 exons and promoter plus the gRNA + iBAR and PuroR
barcode amplicons, and sequencing on a NovaSeq 6000. Read lengths are not published.

🟢 Table S2 also publishes the integrated pKLV2 region, which fixes the transcribed
layout: a vector G, a 19-nt guide spacer, the scaffold loop `GTTTAAGAGCTATGCTGGA`, the
6-nt iBAR, and the scaffold continuing to the oJR160 priming site; and, in the PuroR
transcript, the NeoUTR3 5′ UTR, the 12-nt PuroR barcode and the PuroR start. The
gene-specific arms of both spike-in RT primers are now checked to be the reverse
complements of those vector stretches, and Table S2's amplicon rows name the two
barcode amplicons (`sanger_AMP40` iBAR, `sanger_AMP41` PuroR) beside 39 JAK1 amplicons.

🟡 Model corrections: the guide library no longer carries a single 26-nt "guide"
placeholder but the source-ordered spacer/scaffold/iBAR layout; the PuroR library
carries the published UTR, barcode and CDS start; both nested libraries are built with
the amplicon constructor (no dA ligation junction) and a **10-base** i5 index, because
the printed nested forward primer is exactly P5 + 10N + TruSeq Read 1; the genotype
library formerly used an inferred16-base barcode. Batch06 removes that borrowed physical width entirely; MissionBio identity remains distinct from Chromium identity. These were repository modelling assumptions, not source
errors.

🔴 Still not obtained within the cap: the 10x 5′ HT kit user guide (revision not named
in the Methods), the Mission Bio Tapestri library-preparation manual, and the authors'
code at MarioniLab/scSNV-seq. The i7-side nested primer and the complete bead construct
remain undisclosed in the inspected sources; that is a deferral, not a confirmed absence.


## Native graph migration, 2026-10-18 batch06

Fresh source scope: published PMC10789043 Methods, “gRNA library cloning…” and
“Base editing screens” (text113–121), plus Additionalfile3/TableS2 MOESM3,
Table_S2 customRT/nested primers, JAK1/vector target panel and integratedpKLV2
sequence. Two originals; full bundle scraper run and relevant hits verified.
Superseded preprint chemistry is not current authority. No absence conclusion.

🟢 The source samples separately expanded edited cells for Chromium5′HT and
Tapestri. Genotyping is not DNA split from the same Chromium droplet. The graph
now separates those physical populations before their assay-specific workflows.
MissionBio and Chromium cell identifiers are distinct; linkage through iBAR/PuroR
identifies clonal descendants computationally, not molecular identity of one cell.

Native customRT copies are ss and built from the printed primer plus the copied
upstream vector region; the priming site is not copied twice. Independent guide,
PuroR, otherRNA, gDNA and bead populations are separate Scene strands. After
source cDNA amplification the smallSPRI fraction is retained for both analogous
nested enrichments, with mixed aliquots. The graph exposes first enrichment and
printed-forward nested PCR stages; first-round primer identities/cycles and the
complete native reverse assembly remain unassigned. Actual incoming native
mixtures provide PCR geometry rather than final library references.

The printed P5–N10–TruSeqRead1 nested primer is retained, together with source
19N guide/6N iBAR/12N PuroR target features. Previously borrowed16/10 Chromium
and16base MissionBio widths, fixed i7 widths and complete TruSeq shells are
removed. Native10x transcriptome completion and Tapestri internal preparation
are explicit unread-source boundaries. NovaSeq6000 is positive; native run cycles,
cell/UMI widths and index acquisition remain unassigned. Off-target/RNA/DNA
recovery is not inferred from nested-PCR specificity or computational cell calling.

## Independent implementation correction (2026-10-18)

Confirmed separate expanded-cell sampling for Chromium and Tapestri, correct distinct cell identities/computational clone linkage, exact source RT-target computation and actual mixed precursors for native first/nested PCR. Found/fixed premature Chromium bead row in expanded cells and separately sampled Tapestri states: beads now appear only after GEM partitioning and in downstream copied/amplified mixtures. Native manufacturer and reverse-primer/index-run boundaries remain explicit.
