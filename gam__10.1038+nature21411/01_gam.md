# GAM — research notes

Beagrie et al., *Nature* (2017),
[doi:10.1038/nature21411](https://doi.org/10.1038/nature21411).

## Sources read

- 🟢 Defining article, Methods and Supplementary Guide (PMCID: PMC5366070).
- 🟢 Primary sequencing deposit GSE64881.
- 🟢 The later Multiplex-GAM paper
  ([doi:10.1038/s41592-023-01903-1](https://doi.org/10.1038/s41592-023-01903-1))
  was read only to separate the original one-profile-per-library protocol from the later
  multi-profile variant.

## What constitutes one GAM observation

GAM does not ligate spatial neighbours. A random approximately 220-nm plane through one
fixed nucleus selects many physically separate pieces of chromosomes. One column of the
eventual segregation matrix is the set of genomic windows detected in that nuclear profile
(NP); contacts are inferred because nearby loci co-segregate across many independent NPs.

1. 🟢 Mouse ESC nuclear architecture is strongly fixed, embedded in 2.1 M sucrose, frozen
   and cut into approximately 220-nm ultrathin cryosections.
2. 🟢 Sections are placed on PEN membrane slides. One visually identified NP is cut and
   laser-catapulted into an adhesive PCR cap. One cap in each eight-cap strip is empty and
   follows the entire WGA procedure as a contamination control.
3. 🟢 DNA is extracted and whole-genome amplified with the Sigma WGA4 single-cell kit.
   The defining methods disclose the temperatures, times and reagent roles but not the
   proprietary adapter bases.
4. 🟢 WGA material is purified. A TruSeq DNA HT library uses up to 1.1 µg and is Pippin
   selected to 300–500 nt; the alternative TruSeq Nano DNA HT workflow uses up to 200 ng
   and a bead selection centred near 350 nt.
5. 🟢 Ninety-six indexed NP libraries are pooled. The original mESC data comprise 408
   passing NPs after excluding samples with less than 15% uniquely mapped mouse reads.

## Final molecule and read entry

- 🟢 WGA leaves a 30-base amplification adapter on both ends of each genomic piece.
- 🟢 The final sequencing shell is the named TruSeq DNA HT/Nano chemistry. Its sample
  indexes identify an NP library, not a molecular contact or cell barcode.
- 🟢 Pools of 96 were sequenced "in single end 100 bp rapid-run mode on two lanes of an
  Illumina HiSeq"; because each library has 30 bp WGA adaptors at both ends, "the flow
  cell was not imaged for the first 30bp of each run" (dark cycles; custom recipe with
  Illumina). There is no Read 2 observation in this defining run.
- 🟡 Whether the 30 dark cycles are part of the 100 or precede 100 imaged cycles is not
  stated. An earlier version of this note and the model asserted 30 + 100 = 130 cycles;
  the model now uses the printed 100-cycle Read 1 and leaves the imaged genomic length
  open.
- 🔴 WGA4 adapter base sequences are proprietary. The schematic preserves their published
  30-nt length and marks the bases inferred/unknown rather than inventing an oligo.

## Protocol/version boundary

The later Multiplex-GAM protocol can put several NPs into one library and uses Nextera;
later GAM applications also introduced a public GAT-7N/GAM-COM DeepVent WGA. Those are
useful developments, but neither is silently substituted for the original 2017 GAM page.

## Bounded source pass 2026-10-17 (round breadth-first-10)

The source directory was missing; `tools/get_sources.py --slugs gam__10.1038+nature21411`
fetched the PMC Cloud author-manuscript JATS/PDF, NIHMS supplements and Springer ESM.
Read: PMC5366070 JATS Methods lines 171–201 (cell culture, cryosections, NP isolation,
WGA, library preparation, sequencing, mapping); Supplementary Guide (whole, 31 lines);
Supplementary Table 1 header and kit/index columns. GEO GSE64881 returned a reCAPTCHA
page this pass (manual row) and was not re-read.

- 🟢 Fixation 4% then 8% EM-grade PFA in 250 mM HEPES (10 min, 2 h), 2.1 M sucrose,
  ~220 nm Tokuyasu cryosections on PEN membrane slides, sucrose washed off (PBS, water).
- 🟢 WGA4 with "minor modifications to the previously described protocol" (ref. 38):
  proteinase K in single cell lysis and fragmentation buffer, 50 °C 4 h; 99 °C 4 min;
  library preparation buffer and stabilisation solution, 95 °C 4 min; library
  preparation enzyme 16/24/37 °C 20 min each, 75 °C 5 min; then amplification master
  mix and WGA DNA polymerase with the supplier's program. The model now separates
  lysis from the proprietary library-preparation step (a warned magic boundary).
- 🟢 Supplementary Table 1: 85 NP libraries TruSeq HT and 386 TruSeq Nano, plus 5 Nano
  negative controls, all with 8-nt i7 and i5 indexes; the positive control was
  TruSeq DNA LT with a 6-nt index. Acquired index cycles are not stated.
- Published discrepancy (positive control): the NP-isolation paragraph names
  WGA-amplified E14 mESC genomic DNA, while the WGA paragraph adds 1 µl human genomic
  DNA to an empty lid. Shown on the page; not resolved.
- Open: TruSeq DNA HT/Nano normally fragment input; the Methods still state that each
  library carries 30 bp WGA adaptors at both ends. The model follows that statement.
- Deferred: WGA4 kit technical bulletin; ref. 38 predecessor protocol (Baslan et al., Nat Protoc 2012, doi:10.1038/nprot.2012.039); GEO GSE64881.
