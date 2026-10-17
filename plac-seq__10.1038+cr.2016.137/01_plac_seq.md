# PLAC-seq — research notes

Defining source: Fang *et al.*, [10.1038/cr.2016.137](https://doi.org/10.1038/cr.2016.137).

## Evidence used

- 🟢 PLAC-seq means Proximity Ligation-Assisted ChIP-seq.
- 🟢 Crosslinked nuclei undergo restriction digestion, biotin fill-in and proximity ligation before chromatin shearing.
- 🟢 Immunoprecipitation after ligation enriches contacts anchored at the chosen protein or histone mark; biotin selects contact junctions for paired-end library construction.
- 🟢 An [original PLAC-seq GEO sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2295907) describes end repair, dA addition and Illumina adapter ligation after fragmentation, without naming a custom contact-specific oligonucleotide. The contact junction is created by restriction-end fill-in and proximity ligation, rather than an oligo bridge.
- 🔴 The record does not disclose the lot-specific adapter or index oligo sequences. The final library's standard Illumina shell is a model of the sequencing platform, not a transcription of an exact ordered adapter set from this GEO record.

## Published Data S1 (round 2026-10-17-breadth-first-15)

Sources read: published PMC article text (PMC5143423, lines 1161–1174, Fig. 1A legend
line 1169) and the published Supplementary information obtained from the
nature.com ESM links (the PMC copies sit behind the download gate): Data S1
"Methods & Materials" (`41422_2016_BFcr2016137_MOESM56_ESM.pdf`, text lines 1–73) and
Table S1 (`..._MOESM55_ESM.xlsx`). The bioRxiv preprint 10.1101/074294 is superseded and
was not used for this pass.

- 🟢 Fixation: 1% methanol-free formaldehyde, 15 min RT, quenched with 0.2 M glycine.
- 🟢 In situ proximity ligation "similar to" in situ Hi-C (Rao et al. 2014): 0.5–5 M
  cells lysed (0.2% IGEPAL), 0.5% SDS at 62 °C 10 min, quenched with Triton X-100;
  100 U MboI in NEBuffer 2, 2 h at 37 °C; MboI inactivated at 62 °C for 20 min;
  fill-in with 15 nmol each dCTP, dGTP, dTTP and biotin-14-dATP and 40 U Klenow,
  1.5 h at 37 °C; ligation at RT in 1.2 mL with 4000 U T4 ligase, 1% Triton X-100 and
  0.1 mg/mL BSA (duration not printed).
- 🟢 ChIP: nuclei pelleted, lysed in RIPA, sonicated on a Covaris M220 (75 W, 10%
  duty, 200 cycles/burst, 10 min, 7 °C), cleared at 14,000 rpm; pre-cleared 3 h with
  Protein G Sepharose; ~5% saved as input; 2.5 µg H3K27ac (ab4729) or H3K4me3
  (04-745) or 5 µg Pol II (ab817) antibody overnight; BSA-blocked Protein G beads 3 h;
  washes RIPA ×3, high-salt RIPA ×2, LiCl ×1, TE ×2; RNase A, proteinase K reversal at
  65 °C; phenol/chloroform/isoamyl alcohol and ethanol precipitation.
- 🟢 Biotin pull-down and library construction follow the in situ Hi-C protocol except
  20 µL MyOne Streptavidin T1 beads (instead of 150 µL) and a qPCR-chosen minimal
  PCR cycle number. No adapter or index oligo is printed in Data S1; the end-processing
  and adapter steps live in the cited Rao et al. 2014 protocol, which this pass did not
  inspect. The page shows that step as a warned boundary.
- 🟢 Read processing: paired-end reads, each end mapped separately with BWA-MEM (mm9),
  uniquely mapped (MAPQ > 10) cis pairs kept, ends within 500 bp of an MboI site,
  Picard duplicate removal. Read lengths are not stated in Data S1.
- 🟢 Table S1 sheet 2 lists 4C-seq validation primers (TruSeq-tailed); they belong to
  the 4C comparison, not to the PLAC-seq library.

## Open questions

- Library construction after pull-down (end repair, dA-tailing, adapters, index PCR):
  in the cited Rao et al. 2014 in situ Hi-C protocol, not inspected in this pass.
- Read lengths and adapter/index identities: not in Data S1; GEO GSE86150 run
  metadata not re-read in this pass.

## Protocol boundary

The order is defining: PLAC-seq ligates in intact nuclei before shearing and ChIP. HiChIP remains a separate named protocol.
