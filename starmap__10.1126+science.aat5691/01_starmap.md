# STARmap — source notes

**Defining source:** Wang et al., *Science* (2018),
[doi:10.1126/science.aat5691](https://doi.org/10.1126/science.aat5691).

## Evidence used

- 🟢 The defining paper's supplementary methods describe each SNAIL pair as a
  19–25-nt primer probe and 19–25-nt padlock probe binding adjacent RNA regions,
  with 6-nt complementary probe ends and one ligation site.
- 🟢 The primer probe carries the 12-nt template that receives both 6-nt padlock
  ends. It therefore splints the padlock's single closure junction and remains as
  the primer for rolling-circle amplification; it is not incorporated into the circle.
- 🟢 The same methods specify four SNAIL pairs per gene, a 5-base gene code,
  T4 DNA ligase, phi29 rolling-circle amplification with aminoallyl-dUTP,
  polyacrylamide embedding, proteinase-K clearing, and six SEDAL cycles.
- 🟢 In SEDAL an 11-nt partly degenerate reading probe and an 8-nt fluorescent
  decoding probe form a 19-bp ligated duplex only at a correct interrogation base.
- 🟡 The schematic uses symbolic target arms and code bases while the original
  Tables S1 and S2 are being recovered. The public manuscript supplement's contents
  page and methods explicitly cite those sequence tables, but the available PDF copy
  ends before the tables. The sequences are therefore unverified here, not established
  as unpublished. A later STARmap-family table is a different panel and revision.
- 🟢 A later STARmap-family protocol, doi:10.1038/s41596-025-01248-3, is retained as
  supporting evidence for the named SEDAL reagents, but its oligos are not presented
  as if they were the original 2018 panel.

## Bounded pass 2026-10-17 (round breadth-first-19)

Read the NSF public manuscript bundle (par.nsf.gov 10063217, `10063217.html.txt`):
main text lines 104–150, "STARmap probe design" 1113–1127, thin-section procedure
1128–1180, thick-tissue procedure 1240–1270, Fig. S1 legend 1386–1391 and Fig. S3
legend 1437–1449. The scraper finds no sequences in the defining article or the NSF
bundle; its hits are the 2025 successor panel only.

- 🟢 Thin-section route (lines 1147–1168): probes 100 nM each in 2X SSC/10% formamide
  hybridization buffer, 40 °C overnight, PBSTR and 4X SSC washes; T4 DNA ligase
  (Thermo EL0011) 1:50 RT 2 h; Phi29 (Thermo EP0094) 1:50 with 250 µM dNTP and 20 µM
  aminoallyl-dUTP, 30 °C 2 h; 20 mM acrylic acid NHS ester 2 h; 4% acrylamide/0.2%
  bis with APS/TEMED 1 h; Proteinase K 0.2 mg/ml 37 °C 1 h. The graph now separates
  hybridization (with unbound-probe waste), ligation, RCA, AA-NHS acrylation and
  embedding, and names these reagents.
- 🟢 Thick-tissue route (lines 1245–1264) embeds and digests (overnight) before
  ligation (12 h) and RCA (12–24 h), uses 5′-amino/AA-NHS primer probes and
  BSPEG9 crosslinking, and reads with OR1–OR7 plus four 1-base decoders. Not drawn.
- 🟢 SEDAL (Fig. S3A legend; main text 139–148): reading probe 11 nt, 5′-phosphate and
  3′ inverted dT; 8-nt fluorescent decoding probe with 3′-OH; T4 ligase joins them
  only on a perfect match, leaving a 19-nt product; 60% formamide strips; each cycle
  adds one degenerate base to shift the frame. Main 160/1020-gene runs use R1–R6 and
  sixteen 2-base decoders. Model corrections: reading probe now carries 3′ InvT,
  decoders are described as 2-base (not "one of four"), the rolony is drawn as
  primer-initiated circle complements (it was labelled as circle repeats), the
  input RNA no longer shows probes, and the 5-base gene code carries a
  `feature_barcode` feature.
- 🟢/discrepancy: probe design (line 1118) gives 20–25-nt halves; Fig. S1A gives
  19–25 nt per probe arm. Shown on the page; arms remain symbolic.
- Open: original Tables S1/S2 (bases, gene codes, 2-base decoder assignments); the
  1-nt flanking base of the 17-nt SEDAL unit (Fig. S16 text line 1625) is not drawn;
  which dinucleotide each of the six cycles reads is not modelled.

## Boundary

STARmap produces spatially anchored DNA rolonies read by microscopy and iterative
ligation chemistry. It does not produce an Illumina R1/I1/I2/R2 library.
