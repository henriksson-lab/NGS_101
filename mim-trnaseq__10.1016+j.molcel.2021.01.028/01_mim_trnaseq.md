# mim-tRNAseq

**Defining source:** Behrens *et al.*, *Molecular Cell* (2021),
doi:10.1016/j.molcel.2021.01.028. The PMC full text and Supplementary Table S4 were read.

## Evidence used

- 🟢 Full text lines 82–92 establish barcoded 3′-adapter ligation, prolonged low-salt
  TGIRT reverse transcription and preservation of modification-induced mismatches.
- 🟢 STAR Methods lines 398–410 establish Rnl2tr KQ ligation, the modified RT oligo,
  CircLigase closure, five-to-six-cycle PCR and 150-cycle NextSeq sequencing.
- 🟢 Supplementary Table S4 gives all four 3′ adapters and the RT/PCR oligos verbatim;
  the public page shows a representative adapter and every construct-changing primer.
- 🟢 Table S4 writes `*` between the last two 3′ bases of both library PCR primers
  (`...GCT*C` and `...GT*G`). Supplementary Methods S1 says these phosphorothioate
  bonds are essential because KAPA HiFi has strong 3′→5′ exonuclease activity. The
  page's oligo list preserves the `*`; the DNA-only construct segments represent the
  resulting base sequence and cannot encode a backbone bond.
- 🟢 Table S4 also gives template-switching RNA/DNA starter oligos; the defining
  library branch described in STAR Methods uses the primer-dependent RT oligo with
  5′ `pRN` notation and an internal `iSp18` spacer. These are distinct initiation
  modes and should not be conflated.

## Boundary

Unlike DM-tRNA-seq, mim-tRNAseq does not use AlkB: its information includes the
misincorporation signatures of native modifications.
