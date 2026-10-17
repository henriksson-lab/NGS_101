# NanoSeq — research notes

Defining source: Abascal *et al.*, "Somatic mutation landscapes at single-molecule
resolution", *Nature* 593, 405–410 (2021),
[10.1038/s41586-021-03477-4](https://doi.org/10.1038/s41586-021-03477-4).

## Sources read

- The published article text is behind the publisher page (manifest row `(manual)`); not read.
- Accepted author manuscript (Cambridge repository 1810/319175; local
  `author_manuscript.pdf.txt`): Results lines 264–279 (fragmentation and ddBTP A-tailing),
  Methods "BotSeqS and NanoSeq library preparation protocols" lines 1268–1369, "DNA
  quantification, dilution and PCR amplification" lines 1370–1386 and 1461–1464, and the
  read-processing Methods lines 1590–1622.
- Published Supplementary Information (Springer MOESM1): Note 2 (restriction-enzyme choice,
  lines 63–74), Note 3 (Mung Bean alternative, lines 76–110), Note 4 (modified A-tailing,
  lines 111–128). Supplementary Table 3 (MOESM3 sheet "SupTable3 – Restriction enzymes")
  lists the candidate enzymes only.
- The oligo scraper's hits in MOESM3 are indel alleles from Supplementary Table 6, not oligos.
- The author manuscript cites "Supplementary Note 3" for the dATP/ddBTP mix. In the
  published Supplementary Information that topic is Note 4. This is an ordinary numbering
  change between versions, not a published discrepancy.

## Evidence used

- 🟢 Input: 10 ng genomic DNA (or LCM sections in 20 µl) bound with 100 µl of a 50:50
  water/AMPure XP mix, then on-bead fragmentation with 0.5 µl 5 U/µl HpyCH4V in CutSmart
  (25 µl, 37 °C, 15 min), followed by a 2.5x AMPure cleanup.
- 🟢 HpyCH4V was chosen from 14 candidate 4-bp cutters for its genome and coding coverage.
  AluI and HpyCH4V leave blunt ends (SI Note 2). HpyCH4V cuts at TGCA, so read ends begin
  with `TCA` after the adapter (Methods).
- 🟢 A-tailing: 15 µl reactions with Klenow fragment (3'→5' exo-) in NEBuffer 4 and
  either 1 mM dATP or 1 mM equimolar dATP/ddBTPs, at 37 °C for 30 min. With ddBTPs present,
  extension from an internal nick stops at an incorporated ddBTP, making that strand
  unamplifiable (SI Note 4).
- 🟡 This repository previously called this step "restrictive end repair". The sources
  describe restriction fragmentation without end repair, followed by modified A-tailing.
  The model now uses the source's step.
- 🟢 Ligation: the A-tailing product is added to 22.4 µl of mix with 0.33 µl 15 µM xGen
  Duplex Seq Adapters (IDT 1080799) and 0.56 µl 400 U/µl T4 DNA ligase, at 20 °C for 20 min.
  This is followed by a 1x AMPure cleanup.
- 🟢 Quantification: KAPA library quantification (KK4835) with two added HPLC-purified
  primers, NanoqPCR1 and NanoqPCR2 (verbatim below; no modifications printed). These are
  qPCR oligos, not the PCR primers used to amplify the library.

```
{{oligo: nanoseq_chem.NANOQPCR1}}
{{oligo: nanoseq_chem.NANOQPCR2}}
```
- 🟢 The library is diluted to the target amount (typically 0.3 fmol) and amplified in 50 µl
  with NEBNext Ultra II Q5 Master Mix and dried UDI-containing PCR primers. Cycling is
  98 °C, then 98 °C 10 s / 65 °C 75 s, with 10–16 cycles depending on input. Two consecutive
  0.7x AMPure cleanups follow, and reads are 150 bp paired-end (HiSeq2500, HiSeqX, NovaSeq).
- 🟢 Read structure: reads begin `NNNTCA` or `NNNXTCA`. `NNN` is a 3-nt random barcode, `X`
  an optional spacer, `T` the adapter overhang, and `CA` the HpyCH4V end. Four adapter
  bases are clipped.
- 🟡 Final library: both qPCR primers are the first 20 bases of the TruSeq Read 1 and Read 2
  arms (`illumina.TRUSEQ_READ1`, `illumina.TRUSEQ_READ2`). The model therefore draws an
  *inferred* TruSeq-type shell (P5/i5/Read 1 … Read 2/i7/P7, with 8-nt index windows) around
  `NNN`‑`T`‑`CA…TG`‑`A`‑`nnn`. The optional spacer base is not drawn. The 5′ phosphate on the
  restriction ends is general enzymology, not a sentence in the source.
- 🔴 Not printed in the inspected sources: the complete xGen adapter strands and their
  modifications, the UDI primer sequences, and the index lengths.

## Variants not drawn

- BotSeqS comparison libraries: sonication to 450 bp and direct ligation (Methods).
- Mung Bean NanoSeq: sonication, Mung Bean nuclease, T4 PNK, then dATP/ddBTP A-tailing
  (Methods; SI Note 3).
- The no-ddBTP NanoSeq condition (plain dATP).

## Open questions

- The author Protocol Exchange method (Lensing *et al.*,
  [10.21203/rs.3.pex-1298/v1](https://doi.org/10.21203/rs.3.pex-1298/v1)) has not been
  read. It may print the adapter and UDI details.
- IDT xGen Duplex Seq Adapter documentation, for the strand sequences and the spacer
  design.
- The published Methods text has not been compared with the author manuscript.

## Protocol boundary

NanoSeq is not Nanopore sequencing. It is a distinct Illumina rare-variant protocol and
does not replace the general Duplex Sequencing page.
