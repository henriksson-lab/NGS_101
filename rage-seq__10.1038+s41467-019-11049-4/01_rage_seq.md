# RAGE-Seq — research notes

Defining source: Singh *et al.*, “High-throughput targeted long-read single cell sequencing reveals the clonal and transcriptional landscape of lymphocytes,” *Nature Communications* (2019), [10.1038/s41467-019-11049-4](https://doi.org/10.1038/s41467-019-11049-4) (PMC6635368). This is the published version; the bioRxiv preprint 10.1101/424945 is superseded and not used as evidence.

## Sources read

- 🟢 Published article, PMC JATS text (`s41467-019-11049-4_PMC6635368.xml.txt`): Results “RAGE-seq workflow” (line 72); Methods “Droplet-based scRNA-Seq (10x Genomics)” (166), “Antigen–receptor capture probe design” (170), “Targeted capture” (174), “Nanopore sequencing” (178), “Demultiplexing nanopore sequencing data” (192).
- 🟢 Supplementary Information (MOESM1): Supplementary Fig. 1a legend (lines 75–83) and Supplementary Table 8 (lines 2244–2310; kits and flow cells per sample).

## Evidence used

- 🟢 Chromium Single Cell 3′ v2 libraries were prepared with two extra cDNA amplification cycles; 30 µL cDNA (35 µL for the tumour) went into the Illumina library, sequenced on a NextSeq 500 (150-bp paired-end). The remaining cDNA went to targeted capture.
- 🟢 Before capture, the cDNA was amplified with KAPA HotStart HiFi ReadyMix, 3 µM TSO primer `AAGCAGTGGTATCAACGCAGAGT` and 3 µM R1 primer `CTACACGACGCTCTTCCGATCT` (5 cycles for cell lines, 20 for primary cells), then purified with AMPure XP. 500 ng–1 µg went into capture.
- 🟢 Capture followed the Roche–NimbleGen SeqCap EZ double-capture protocol: two overnight hybridizations at 47 °C, with no universal or index HE oligos. After each capture, a PCR with 1 µM TSO and R1 primers (instead of the TS-PCR oligos) ran 5 cycles after the first capture and 20 after the second. The post-capture library was 0.6–2 kb.
- 🟢 The SeqCap RNA Choice design targeted 644 exons (~128 kb) of functional human IG/TCR V, J and constant genes; Supplementary Table 7 lists the genes. Roche–NimbleGen designed and synthesized the probes. Exact individual probe sequences are unassigned in this model; the linked design/vendor source chain remains outside this capped pass.
- 🟢 Long-read libraries used the ONT 1D ligation kit SQK-LSK108, or SQK-LSK308 for one sample (base-called as 1D). Flow cells were R9.4.1 FLO-MIN106, or R9.5.1 FLO-MIN107 for three cell-line samples. The native adapter/end-preparation geometry remains unassigned; the historical kit protocols were not read.
- 🟢 Demultiplexing: 16-nt cell barcodes from the matched CellRanger output are searched, on either strand, in the first and last 200 nt of reads longer than 250 nt. The 13 nt downstream of a barcode match are removed to drop the 10-nt UMI.
- 🟡 Model: the amplicon is drawn as R1 primer site – 16-nt cell barcode – 10-nt UMI – oligo(dT) – receptor cDNA – template-switch-derived bases – TSO primer site. The flanks come from the two printed primers, and the barcode/UMI order and lengths come from the demultiplexing Methods. The oligo(dT) and the remaining 10x template-switch bases are role-level.

## Protocol boundary

The page covers the targeted receptor branch: re-amplification, double capture, post-capture PCRs and the ONT library. The matched short-read branch follows the standard Chromium 3′ v2 protocol. It is drawn as an explicit, warned boundary rather than as a redrawn 10x library.

## Open questions

- The ONT SQK-LSK108/SQK-LSK308 kit protocols (end-prep/dA-tailing and adapter-end chemistry) were not read.
- The NimbleGen SeqCap EZ User's Guide v3.0 and Double Capture Technical Note, cited in Methods, were not read. They describe hybridization, washes and the replaced TS-PCR oligos.
- The 10x Chromium 3′ v2 user guide (bead-primer oligo(dT), TSO, Illumina library) was not read for this page.

## Native graph migration — 2026-10-18 batch04

One original was re-read: the published2019 Methods and demultiplexing sections. The graph begins with the explicitly source-amplified Chromium3′v2 full-length cDNA pool, retaining the upstream commercial RT/amplification boundary. The expression and targeted-capture aliquots carry the same transcriptome mixture. R1/TSO primers remain in the local source-owned module; their ends and source16nt barcode/10nt UMI lineage remain explicit. PCR geometries receive the actual incoming amplified or capture-bound material. Capture-bound DNA, probe and background rows are separate; wash transfer is physical but does not assert pure receptor recovery. BaselineLSK108 native joining stays unassigned; the one LSK3081D2 sample analysed as1D is retained as a separate historical alternative. Source NextSeq500150bps paired-end wording is positive but does not assign an exact per-read/index mask here: the read table leaves cycles unassigned rather than inventing150/150. Linked historical SeqCap/ONT/vendor recipes were not followed within this capped pass, and no publication-wide absence is concluded.

The local curator-written review target tuple is independent of graph contents. This capped pass does not complete the missing-information audit and accepts no source-absence claim.
