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
- 🟢 STAR Methods line 404 prints both the manufacturer's high-salt TGIRT buffer
  (450 mM NaCl) and the low-salt buffer (75 mM KCl) used for the improved
  primer-dependent branch. The 2017 InGex small-RNA/CircLigase protocol describes
  a high-salt template-switching branch and cannot establish this low-salt recipe.
- 🟢 The distributor-hosted InGex [TGIRT50.pdf](https://fnkprddata.blob.core.windows.net/domestic/data/datasheet/IGX/TGIRT50.pdf)
  and [TGIRT10.pdf](https://fnkprddata.blob.core.windows.net/domestic/data/datasheet/IGX/TGIRT10.pdf)
  are byte-identical five-page copies of that 2017 small-RNA/CircLigase protocol
  (MD5 `4ca4abd4a8e6ef898670ea8d2cac9051`). The TGIRT50 filename does not
  identify a distinct enzyme product sheet or the authors' lot.
- 🟢 Supplemental mmc2–4 are result matrices (Tables S1–S3), and mmc6 is a
  one-page cover; the chemistry methods and oligos are in mmc1 and mmc5.
- 🔴 The historical InGex TGIRT50 vendor revision used by the authors has not
  been obtained. Enzyme-dependent details beyond those printed in the original
  article and supplement remain unresolved.

## Boundary

Unlike DM-tRNA-seq, mim-tRNAseq does not use AlkB: its information includes the
misincorporation signatures of native modifications.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC8062790 JATS376–410, tRNA purification, barcoded adapters, RT and cDNA circularization/PCR.

Corrected the displayed RT primer from invented RNA rNrN to source-written degenerate RN. Changed NaOH RNA hydrolysis from strand separation to chemical conversion; specified CircLigase60°C3h and KAPAHiFi5–6 PCR cycles.

Deferred: MethodsS1 extended protocol, adapter preadenylation with Mth ligase, preparative deacylation/PNK/gel stages, RNA-fragment fates and run mask; no universal charging-assay branch substituted. This pass does not establish a complete source review or publication-wide absence.

## Native graph migration (2026-10-18, batch02)

Inspected: Published PMC8062790 JATS Results84/92, STAR Methods382–410 and preprocessing424; published TableS4 mmc5 XLSX oligo rows5–13. Mandatory whole-bundle scraper reviewed and modified oligos checked directly.

Selected primer-dependent low-salt 42°C/16h I1 branch is explicit, with I2–I4 and alternative TS/high-salt/short-RT/charging branches separate. Native graph includes adapter preadenylation, physical terminal ddC, source ligation/PAGE, RT and persistent RNA until NaOH hydrolysis, actual spacer-containing ss-cDNA and ss-circle, one-fifth circle aliquot and actual-circle PCR input. Non-nucleotide iSp18 remains physical; the blocked adapter suffix beyond the RT-binding site is not automatically copied into cDNA. Full circle-PCR molecular proof is unresolved; final candidate roles do not become the initial template. Narrow adapter source-guard exemption documents scraper omission of terminal C in /ddC/, independently covered by the complete modified-oligo regression.

Deferred: Full ss-circle amplification geometry, complete internal final construct, acquired read/index cycle mask beyond reported total150 cycles; MethodsS1 not newly re-read and no vendor recipe borrowed. This capped pass does not establish publication-wide absence or complete molecular proof.
