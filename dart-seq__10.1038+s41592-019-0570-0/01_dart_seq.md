# DART-seq — research notes

Defining source: Meyer KD, "DART-seq: an antibody-free method for global m6A detection",
Nat. Methods 16, 1275–1280 (2019), [10.1038/s41592-019-0570-0](https://doi.org/10.1038/s41592-019-0570-0)
([PMC6884681](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884681/), author manuscript).

## Sources read (bounded pass 2026-10-17)

- PMC6884681 JATS text: Results (DART-Seq concept, low input, long-read), Methods "Plasmids",
  "In vitro deamination assays", "cDNA synthesis and Sanger sequencing", "Next-generation
  sequencing", "Data analysis", Data availability (GEO GSE125780).
- Supplementary Protocol "Transcriptome-wide m6A detection with DART-Seq" (Springer MOESM10;
  identical file as PMC NIHMS1537259-supplement-2), sections General notes, A–D.
- Reporting summary (MOESM2) for platform naming only.
- Oligo scraper over the whole bundle: hits are cloning primers (YTH Fwd/YTH-HA Rev,
  YTHmut Fwd), the METTL3 sgRNA, the biotin-RNA pulldown bait and Sanger RT/PCR validation
  primers. No sequencing-library oligo is printed.

## Evidence used

- 🟢 DART-seq means deamination adjacent to RNA modification targets.
- 🟢 Rat APOBEC1 (pCMV-APOBEC1, Addgene #73019) is fused through a 15-amino-acid linker to
  YTHDF2 amino acids 385–579 with a C-terminal HA tag; YTHmut lacks aa 385–409 (Methods, Plasmids).
- 🟢 HEK293T cells are transfected with APOBEC1-YTH, APOBEC1-YTHmut or APOBEC1 alone using
  FuGENE HD and grown 24 h (Methods NGS; Supplementary Protocol A).
- 🟢 Total RNA: TRIzol then RNase-free DNase I, 20 min at 37 °C (Methods); the Supplementary
  Protocol B alternatively uses RNeasy Plus Mini with on-column or post-purification DNase.
- 🟢 1 µg total RNA into the NEBNext Ultra II Directional RNA Library Prep Kit for Illumina;
  low input (10 or 100 ng) uses the NEBNext Single Cell/Low Input RNA Library Prep Kit; rRNA
  removal "not necessary"; samples barcoded with NEBNext Multiplex Oligos for Illumina;
  sequenced on HiSeq 4000 (Methods NGS).
- 🟢 Supplementary Protocol C: paired-end or single-read possible; authors "typically use SR
  50bp sequencing on an Illumina HiSeq4000", 30–50 M reads per sample.
- 🟢 Strand-specific reads are reverse-complemented before alignment; C-to-U sites are called
  with CTK/CIMS and filtered against dbSNP, untransfected-cell editing and APOBEC1-alone sites
  (Methods, Data analysis; Supplementary Protocol D).
- 🟡 With a dUTP-type directional kit, single Read 1 is antisense to the RNA, so a C-to-U mark
  appears as G-to-A in the raw read; this is why the reads are reverse-complemented. The kit
  chemistry itself was not read in this pass.
- 🔴 (not printed in the inspected sources) adapter, index-primer and universal-primer
  sequences, which NEBNext Multiplex Oligos index set was used, index-read cycles.

## Model and page

The informative chemistry happens in living cells before library construction. Editing is
not a selection step: RNA without YTH-bound m6A remains in the total RNA and is sequenced;
background is removed computationally using the APOBEC1-alone and YTHmut controls. The page
models transfection, in-cell deamination, extraction and DNase digestion, then one warned
kit-internal step and an index PCR whose primers are named but not sequenced. The final
TruSeq-compatible arms, i5/i7 windows and their lengths are inferred placeholders, marked as
such; only Read 1 is shown as acquired (SR50).

Not drawn: the low-input kit branch, PacBio Iso-Seq long-read DART-seq (1 µg RNA, Sequel),
the in vitro deamination assay (TNT-expressed fusion, 4–6 h at 37 °C, then Ultra II
Directional library) and the RT-PCR/Sanger validation.

## Open questions (deferred, not confirmed absent)

- NEBNext Ultra II Directional (E7760) and Single Cell/Low Input (E6420) manuals, revision
  current in 2019: fragmentation, adaptor (U-loop/USER) geometry and PCR primers.
- NEBNext Multiplex Oligos set used and index-read lengths (possibly in GEO GSE125780 run
  metadata, not inspected).
- Iso-Seq library chemistry for the long-read branch.

## Protocol boundary

This is Meyer's m6A-editing DART-seq, distinct from the unrelated droplet-assisted
RNA-targeting method also called DART-seq.
