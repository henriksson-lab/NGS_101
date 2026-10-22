# TT-seq

## 1. Protocol boundary

TT-seq is the Schwalb et al. workflow: a five-minute 4-thiouridine pulse, spike-in addition,
RNA fragmentation **before** labelled-RNA purification, DNase treatment and a strand-specific
commercial library preparation. The ordering of fragmentation and selection distinguishes it
from the unfragmented 4sU-seq control in the same experiment.

## 2. Sources read

🟢 Schwalb et al., *TT-seq maps the human transient transcriptome*, *Science* (2016),
DOI `10.1126/science.aad9841`.
The original publisher PDF was read from the author-hosted copy listed in `ref/MANIFEST.md`.

🟢 The defining paper's supplementary methods (`aad9841-Schwalb-SM.pdf`), retrieved from
the supplementary-material mirror listed in `ref/MANIFEST.md`. The RNA-seq methods state:
500 µM 4sU for 5 min; spike-in addition during TRIzol extraction; BioRuptor fragmentation at
240 ng/µl for one 30 s ON / 30 s OFF cycle; labelled-RNA purification; 2 U Turbo DNase; NuGEN
Ovation Human Blood RNA-seq kit; HiSeq 1500. The processing section reports paired-end 50-base
reads and an additional six-base barcode read.

🟢 Dölken et al., *High-resolution gene expression profiling for simultaneous kinetic
parameter analysis of RNA synthesis and decay*, *RNA* (2008), DOI
`10.1261/rna.1136108`, the labelled-RNA purification method cited as reference 26 in
the TT-seq supplement. Its original publisher PDF (Methods, p. 12) was read from the
university-hosted copy listed in `ref/MANIFEST.md`.

## 3. Molecular state changes

🟢 Newly synthesized RNA incorporates 4sU during the pulse. The six synthetic spike-ins include
three made with a 1:10 4sUTP:UTP ratio and three made with UTP only.

🟢 Purified RNA is split. The TT-seq arm is fragmented before the labelled fraction is purified;
therefore a fragment is selected only when it contains a pulse-labelled nucleotide. This removes
the positional bias produced when an entire nascent transcript is captured through one labelled
region.

🟢 The supplement delegates labelled-RNA purification to its reference 26. That method
names EZ-Link Biotin-HPDP for 4sU thiol labelling, removes excess biotin reagent,
captures labelled RNA on µMACS streptavidin columns, washes, then releases it with DTT.
The model follows this inherited capture chemistry without assigning its 2008 reagent
volumes to the 2016 experiment.

🟢 Captured RNA is DNase-treated (2 units Turbo DNase) and passed to the strand-specific
NuGEN kit. The reaction graph carries this as its own digestion step between DTT release
and the kit's reverse transcription; no separate cleanup is printed before the kit.

🟡 The supplement does not print the NuGEN adapter oligos. The public page therefore marks its
standard Illumina-compatible final construct and primer landing sites as inferred. The six-base
single-index and paired-end run structure itself is directly reported.

## 4. Read layout

- Read 1: 50 bases of strand-specific cDNA.
- Index 1: six barcode bases.
- Read 2: 50 bases from the opposite cDNA end.
- No Index 2 read is reported for the defining experiment.

## 5. Open evidence boundary

The exact discontinued Ovation Human Blood RNA-seq adapter and index-primer sequences
remain unresolved because its relevant original kit guide was unavailable in this audit.
The inferred Illumina geometry is not an assertion about those historical kit oligos;
in particular its single dA junction base is a TruSeq-style placeholder needed for the
standard TruSeq Read 2 / Index 1 landing, not an A-tailing step printed by any inspected
source. The strand read by Read 1 and which index read (i7) carries the six-base barcode
are likewise not stated in the supplement's processing section.

The 2026-10-17 capped pass (round 21) re-read supplement Materials and Methods lines
21–49 and the original article p.1 lines 88–103, and checked Dölken 2008 Methods for
Biotin-HPDP, chloroform/isoamyl-alcohol removal, µMACS capture and DTT release; the
scraper's supplement hits (lines 202–460) are reporter-cloning primers, not library oligos.

## Supplied originals checked, 2026-10-22

The uploaded TT-seq article PDF was checked against its printed title and own DOI;
its first-page TT-seq description confirms fragmentation before labelled-RNA isolation.
It also contains the preceding article's tail, which is not TT-seq evidence. The
supplied supplement is byte-identical (MD5 recorded in MANIFEST.tsv) to the existing
published `aad9841-Schwalb-SM.pdf`; no duplicate original is retained.

Read supplementary Methods RNA-seq, spike-in mix preparation and Sequencing data
processing (printed pp2–3; text lines21–49). Ran the full bundle primer scraper;
its reporter-assay/termination-primer hits do not identify the commercial sequencing
adapters. Existing source-scoped chemistry agrees with this supplied original.
The typed sequencing plan now carries the source50/6/50 run lengths as structured
data, instead of leaving those values only in prose. The inferred NuGEN-compatible
adapter shell is still explicitly inferred; the supplied documents do not establish
its exact historical kit sequences. This is a bounded check, not an accepted absence audit.
