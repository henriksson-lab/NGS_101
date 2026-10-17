# PHAGE-ATAC

**Defining source:** Fiskin E. et al., "Single-cell profiling of proteins and chromatin
accessibility using PHAGE-ATAC", *Nature Biotechnology* 40, 374 (2022),
doi:10.1038/s41587-021-01065-5. Read as the author manuscript PMC10549956 (JATS XML from
the PMC Cloud bucket) and the published Supplementary Information PDF (Springer
`41587_2021_1065_MOESM1_ESM.pdf`). The bioRxiv preprint (doi:10.1101/2020.10.01.322420)
is superseded and no longer supplies active claims; earlier review entries that read it
are history.

## Sources read (bounded pass 2026-10-17)

- PMC10549956 Methods, text lines 106–108 (PAC-tag/RD1 phagemid cloning; EF79 introduces
  the PCR handle for PDT amplification), 136–148 (staining, fixation, lysis, 10x loading,
  linear PCR, split, PDT library, sequencing), 154–156 (PDT read processing).
- Supplementary Table 1 (MOESM1 text lines 10–193, EF147/EF91/EF77–79/EF170) and
  Supplementary Table 4 (lines 1174–1200, CDR3 barcode readouts with R1/R3 assignment).
- `scrape_primers.py` run on the source bundle; EF91 found once, exact, at MOESM1 line 81.

## Molecular boundary

🟢 Nanobodies are displayed as Nb–PAC–p3 fusions; the 12-amino-acid PAC tag is the
reading-frame-modified RD1, introduced so that the 10x scATAC bead oligo captures the
CDR3 barcode alongside accessible chromatin during droplet indexing (Methods line 108).

🟢 Cells are stained with phage (20 min, 4 °C), fixed in 1% formaldehyde (10 min RT,
glycine-quenched) and permeabilized in a custom lysis buffer (10 mM Tris pH 7.5, 10 mM
NaCl, 3 mM MgCl2, 0.1% NP-40, 1% BSA; 3 min on ice) because the standard 10x scATAC
lysis buffer lost pNb staining. Intracellular staining instead fixes and lyses first,
then stains. Tagmentation, chip loading and encapsulation follow the Chromium Single Cell
ATAC Solution protocol.

🟢 After linear PCR the emulsion is broken, products are cleaned with MyOne silane beads
and eluted in 40 µl; 5 µl goes to the PDT library and 35 µl to the scATAC library.

## PDT library (published Methods)

| PDT PCR oligo | 5′→3′ sequence (Supplementary Table 1) |
|---|---|
| EF147 | AATGATACGGCGACCACCGAGA |
| EF91 | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGATACCGCGGTGTATTATTGC |

🟢 PDT PCR: 100 µl, 2× KAPA polymerase, EF147 + EF91; 95 °C 3 min, 20 × (95 °C 20 s,
60 °C 30 s, 72 °C 20 s), 72 °C 5 min. Double-sided SPRIselect: 65 µl beads, keep 160 µl
supernatant, add 192 µl beads, two 80% ethanol washes, elute 40 µl EB.

🟢 Indexing PCR: 15 ng, 100 µl with 50 µl 10x Amp-Mix, 7.5 µl SI-PCR Primer B and 2.5 µl
10x i7 sample-index primer; 98 °C 45 s, 6 × (98 °C 20 s, 67 °C 30 s, 72 °C 20 s),
72 °C 1 min; 120 µl SPRIselect, elute 40 µl EB; Qubit HS and Bioanalyzer HS QC.

🟡 EF91's 3′ 21 nt (GATACCGCGGTGTATTATTGC) occur inside EF79 (handle mutagenesis primer) and begin every
R3-read CDR3 barcode in Supplementary Table 4, so EF91 primes at the FR3 handle and
adds a Nextera ME/Read 2 tail; EF147 is a P5 primer. The resulting region order drawn on
the page is derived, not printed.

🟢 Supplementary Table 4: CDR3 barcodes are read in R3 for most clones and in R1 for
hashtags PH-C and PH-D.

## Sequencing

🟢 PDT and ATAC libraries pooled; paired-end 2 × 34 cycles, NextSeq 550 High Output;
demultiplexed with CellRanger-ATAC mkfastq into R1/R2/R3; PDT counts by kallisto|bustools
with a dummy poly(A) UMI (phages have no UMI). Index-read cycles are not stated.

## Remaining questions

- 🔴 in the paper: sequences of the 10x ATAC gel-bead oligo, Amp-Mix, SI-PCR Primer B and
  i7 sample-index primers; the version-matched Chromium Single Cell ATAC guide is a named
  follow-up, not read here.
- Exact PDT template junctions after bead capture of the PAC-tag RD1 (depends on the
  bead oligo).
- The PMC "Integrated supplementary information" PDF is download-gated; the publisher
  MOESM1 was read instead. MOESM2/MOESM3 (image-only/near-empty text) not inspected.
- The ATAC library steps after the split follow the vendor protocol and are drawn at role
  level only.
