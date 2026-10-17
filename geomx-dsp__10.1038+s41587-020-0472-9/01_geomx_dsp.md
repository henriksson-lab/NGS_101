# GeoMx DSP NGS readout — source notes

## Scope

This page covers the named **GeoMx Digital Spatial Profiler (DSP) NGS readout**: a
photocleavable synthetic tag is released from an RNA probe or antibody inside a selected
area of illumination (AOI), collected into one well, PCR-completed with well-specific
GeoMx Seq Code primers and counted by short paired-end Illumina sequencing. It does not
describe the alternative nCounter readout.

## Sources read

- 🟢 Merritt et al., *Nature Biotechnology* (2020),
  [doi:10.1038/s41587-020-0472-9](https://doi.org/10.1038/s41587-020-0472-9), defines DSP
  with antibody or RNA probes attached to oligonucleotide tags through a photocleavable
  linker and spatially patterned UV release. The publisher PDF (version of record,
  printed doi, 22 pp.) was obtained on the institutional network on 2026-10-17 and its
  Online Methods were read (see the deferred-followup-03 pass below). The bioRxiv
  559021v2 preprint is superseded and not used for page facts.
- 🟢 NanoString, *GeoMx DSP NGS Readout User Manual: Library Preparation & GeoMx NGS
  Pipeline*, MAN-10153-01 (April 2022), is the authoritative commercial library-prep
  source. Figure 2 gives every region and length; Tables 6–9 give PCR, cleanup and run
  settings.
- 🟢 The current Bruker GeoMx knowledge base independently recommends paired-end 2×27
  sequencing and reports an expected 162-bp Seq Code product.

## Published construct facts

🟢 Figure 2 (MAN-10153-01 p. 12, inspected as the rendered page) labels the collected
aspirates "Photocleaved DNA oligos" and prints SPR1 20 nt, UMI+RTS ID 26 bp, SPR2 20 nt
before PCR and P5 29, i5 8, SPR1 33, UMI 14, RTS ID 12, SPR2 34, i7 8, P7 24 nt after PCR.

🟢 Before PCR the collection plate is dried down (bench overnight, or 65 °C for 1 h with
the lid open) and each well is rehydrated in nuclease-free water (Table 5: 10 µL for RNA,
80 µL for protein; 10 min at room temperature). Each 10-µL PCR then combines 2 µL 5X PCR
Master Mix, 4 µL Seq Code primer mix from the matching primer-plate well and 4 µL aspirate
(Table 7). After PCR, 4 µL of every well (including the NTC, excluding any positive
control) is pooled; the NTC and optional positive control are also purified separately as
QC comparators.

🟢 Before PCR, the 5′→3′ photocleaved tag is `SPR1(20 nt) – UMI(14 nt) – RTS ID(12 nt) –
SPR2(20 nt)`. The 12-nt readout-tag sequence identifier names the biological target. The
14-nt random UMI is used to remove PCR duplicates.

🟢 GeoMx Seq Code PCR extends SPR1 from 20 to 33 nt and SPR2 from 20 to 34 nt, and adds
P5(29), i5(8), i7(8) and P7(24). Each well receives a unique dual-index pair, so the index
pair carries the AOI/collection-well identity. Eighteen amplification cycles follow a
30-min 37 °C UDG incubation and 10-min 50 °C deactivation. Products from the plate are
pooled and purified with two 1.2× AMPure XP cleanups.

🟢 The commercial run (MAN-10153-01) is paired-end Read 1=27, Index 1=8, Index 2=8, Read 2=27. Read 1 therefore
contains UMI cycles 1–14 and RTS-ID cycles 15–26. Read 2 sees the reverse end: RTS-ID
cycles 1–12 and UMI cycles 13–26. The 27th base on each side enters the opposite PCR arm,
which permits the two reads to be stitched.

## What is not public

🔴 MAN-10153-01 prints region lengths and chemistry, but not the bases of SPR1, SPR2,
UMIs, target RTS-ID whitelist, Seq Code primers or UDI sets. The assay-specific RTS map
is supplied in a Probe Kit Configuration file, while the run-specific index values are
supplied in `SeqCodeIndices.csv`.

🟡 The final arm lengths and standard four-read Illumina geometry match canonical TruSeq
P5, Read 1, Read 2, P7 and index-primer regions. The schematic uses those canonical bases
only to compute and draw sequencing-primer binding, and marks every such segment inferred;
it does not claim that NanoString published the proprietary oligo sequences.

## Bounded source pass 2026-10-17 (round breadth-first-10)

Inspected MAN-10153-01 lines 555–660 (introduction, Figure 2 rendered p. 12, Table 1 i5
orientation), 1100–1460 (dry-down/rehydration Table 5, PCR Table 6/7, pooling and both
AMPure cleanups, Table 8 elution) and 1599–1607 (run lengths); Nature HTML preview
(abstract only) and MOESM1/MOESM2 headings. The graph now carries the dry-down and
rehydration step, the UDG pre-incubation in the PCR label and master-mix input, and the
4-µL-per-well pooling. Deferred: published Methods (paywalled), the bioRxiv preprint
NGS Methods, the 5X PCR Master Mix composition and UDG substrate, and the proprietary
tag/Seq Code bases. No missing-information claim is accepted by this capped pass.


## Bounded source pass 2026-10-17 (round deferred-followup-02)

Routes for the published Merritt 2020 article: OpenAlex and Unpaywall report it closed
with no repository full text; Europe PMC has no PMCID (PMID 32393914 only), so no PMC or
PMC Cloud author manuscript exists; a web search found only the bioRxiv versions and a
vendor white paper. The Methods remain paywalled (manual URL above). Lead only, not used
for page facts (superseded preprint; published Methods unreadable, so omission cannot be
verified): bioRxiv 559021v2 "NGS readout for photocleaved oligo counting" (pdf.txt lines
854–865) describes the prototype run as 4 µL aspirate, 1 µM i5 and i7 primers, 1X NSTG
PCR Master Mix, 37 °C 30 min, 50 °C 10 min, 95 °C 3 min, 18 cycles, two 1.2× AMPure XP
cleanups, and paired-end 2×75 on a MiSeq. This agrees with MAN-10153-01 except for the
run length (2×27 in the commercial manual); the page keeps the commercial manual's
settings. No page change. Still deferred: published Methods, master-mix composition,
proprietary tag/Seq Code bases. No missing-information claim is accepted.


## Bounded source pass 2026-10-17 (round deferred-followup-03)

Source: Merritt et al. 2020 publisher PDF `s41587-020-0472-9.pdf` (hand-downloaded by the
maintainer; printed doi:10.1038/s41587-020-0472-9, title and author list checked). Read
in its `.txt`: main text NGS readout paragraph and Fig. 6 legend (lines 2178–2191,
2704–2718); Online Methods "RNA profiling probe design", "Sample preparation for RNA
profiling" and nCounter hybridization (lines 3402–3440); "Next-generation sequencing
readout for photocleaved oligonucleotide counting" (lines 3526–3545). The oligo scraper
returns 0 sequences for the bundle, including this PDF.

- 🟢 RNA probes are, 5′→3′, a 35–50-nt target-complementary sequence, an IDT iSpPC
  UV-photocleavable linker and the indexing oligonucleotide; up to 13 probes per target.
  The NGS indexing oligonucleotide carries PCR primer sites, a UMI and a target barcode;
  each tile (probe) has its own barcode. No bases are printed.
- 🟢 NGS PCR: each ROI uniquely i5×i7 dual-indexed; 4 µL aspirate, 1 µM i5, 1 µM i7, 1×
  NSTG PCR Master Mix; 37 °C 30 min, 50 °C 10 min, 95 °C 3 min, 18 × (95 °C 15 s, 65 °C
  60 s, 68 °C 30 s), 68 °C 5 min; two 1.2× AMPure XP cleanups. This agrees with the
  commercial manual's program and cleanups (the graph is unchanged).
- 🟢 Published runs were paired-end MiSeq 2×75 or NextSeq 2×37, whereas MAN-10153-01
  specifies 2×27 + 8 + 8. This difference is now a graph-owned source note on the page;
  the model keeps the commercial 2×27 run. UMI deduplication used UMI-tools at Hamming
  distance 3 (analysis, not chemistry).

Still open: the master-mix composition and UDG substrate (the article names only "NSTG
PCR Master Mix"), the proprietary SPR1/SPR2, RTS-ID, Seq Code primer and UDI bases, the
article's index-read lengths, and the NTC/positive-control QC branches. No
missing-information claim is accepted by this capped pass.
