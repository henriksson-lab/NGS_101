# SELECT-seq — preprint source intake

Defining preprint: *SELECT-seq allows Pre-Sequencing Enrichment of SNP Edits in One-Pot
Single-Cell Whole-Transcriptome Sequencing*, DOI
[10.64898/2026.07.07.736965](https://doi.org/10.64898/2026.07.07.736965).

## Sources and molecular path

🟢 The July 2026 bioRxiv v1 [original](https://www.biorxiv.org/content/10.64898/2026.07.07.736965v1)
is available. Its Methods (“SELECT-seq workflow,” extracted lines 619–679) say that
single-cell lysis and RT follow [Smart-seq3 Protocol V.3](https://doi.org/10.17504/protocols.io.bcq4ivyw).
Whole-transcriptome PCR retains that chemistry while adding SNP-selective primers.
These have 3′ phosphorothioate linkages and a terminal base matching the target SNP.
Cas12a and a matched crRNA detect the resulting amplicon with a quenched fluorescent
ssDNA reporter before the whole-transcriptome product is purified, tagmented and
indexed for paired-end 150-cycle sequencing.

🟢 The original's [Supplementary Table 1](https://www.biorxiv.org/content/10.64898/2026.07.07.736965v1.supplementary-material)
(`media-2.xlsx`, sheets `primers`, `crRNA`, `pegRNA`, `probes`) prints ordering
sequences. For the NRF2 T80K demonstration it lists forward primers
`ttcgctcagttacaactagatgaaga*g*a*c` (WT) and
`ttcgctcagttacaactagatgaaga*g*a*a` (T80K), reverse primer
`TCAACAGGGGCTACCTGAGC*A*A*C`, and an NRF2 Cas12a crRNA with direct repeat
`TAATTTCTACTAAGTGTAGAT` and spacer `AGGTGACTGAGCCTGATTAGT`.
The Methods' fluorescent reporter is
`/56-FAM/TTATTTTAT/ZEN/TTT/3IABkFQ/`. Asterisks mark published
phosphorothioate linkages; uppercase/lowercase is preserved as printed.

🟢 The cited Smart-seq3 V.3 protocol prints a 5′-biotin oligo-dT30VN,
5′-biotin N8/rGrGrG template-switch oligo, and amplification primers in its
“List of oligos” (PDF text lines 76–120). Those are inherited rather than
SELECT-specific. The exact index primer member used for each SELECT-seq well is
not given in the inspected Methods and table.

🟡 This first pass checked the original PDF, supplementary figures/table and linked
Smart-seq3 V.3 protocol with the oligo scraper. An independent second source pass
is needed before declaring any sequence absent or the protocol fully reviewed.

## Reaction detail (bounded pass 2026-10-17, round breadth-first-18)

Read in the defining bioRxiv v1 PDF text, section “SELECT-seq workflow”
(lines 619–679) and “Purification and library preparation”, “Genotyping of
Cas12a-positive single cells”, “Library purification and sequencing”
(lines 680–813):

- 🟢 Lysis and RT follow Smart-seq3 v3 unchanged. The multiplex PCR is the
  Smart-seq3 v3 whole-transcriptome PCR with the forward and reverse WTA
  primers *reduced* to 0.05 µM final, SNP-specific forward and reverse primers
  added at 0.1 µM final each, the same thermal profile, and the cycle number
  raised to 35.
- 🟢 Cas12a detection is run on the unpurified PCR product: 5 µL of mix giving
  50 nM Cas12a (NEB M0653), 100 nM crRNA, 1 µM double-quenched FAM/ZEN/IBFQ
  ssDNA reporter and 1× NEBuffer r2.1; 37 °C for 30 min with fluorescence
  monitoring, then 65 °C for 10 min to inactivate Cas12a, then −20 °C storage.
  Wells are called SNP-positive or negative on a ΔF threshold. The reporter is
  a fluorescence readout, not a library component.
- 🟢 Purification: 10 µL of reaction + 30 µL water + 24 µL AMPure XP, per the
  Smart-seq3 v3 protocol; Qubit dsDNA HS quantification; dilution to
  100 pg µL⁻¹; Smart-seq3 v3 tagmentation terminated with 0.1 % SDS; index PCR
  for 14 cycles. Indexed libraries are AMPure-purified and sequenced on the
  Element Aviti platform, paired-end 150 cycles.
- 🟢 Inherited Smart-seq3 v3 oligos, verbatim with modifications
  (“List of oligos”, PDF text lines 76–120):
  `Smartseq3_OligodT30VN` = `/5Biosg/ACGAGCATCAGCAGCATACGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN`;
  `Smartseq3_N8_TSO` = `/5Biosg/AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG`;
  `Fwd_PCR_primer` = `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G`;
  `Rev_PCR_primer` = `ACGAGCATCAGCAGCATAC*G*A`. The protocol states custom
  Nextera index primers are used but assigns no member to a SELECT-seq well.
- 🟢 Supplementary Table 1 `primers` lists four allele variants
  (`_F_a/_t/_g/_c`) per SNP target plus one reverse primer, each with three
  terminal phosphorothioate linkages (e.g. `GRIK1_F_a`
  `gcaacaaagacaagtccagca*a*t*a`, `GRIK1_R` `ggtgggagatcttcagaaca*a*c*g`), an
  ACTB internal control pair, and the NRF2 WT/T80K pair already quoted above.
  The `crRNA` sheet gives four PIK3CA crRNAs and one NRF2 crRNA, all with the
  same 21-base direct repeat `TAATTTCTACTAAGTGTAGAT`. The `probes` sheet gives
  six reporters (FAM/HEX/SUN × `TTATT` or `TTATTTTAT/ZEN/TTT`).
- 🟢 The three-step genotyping PCR (first PCR with
  `cgacggaaagagtatgagctgg` / `gcaacctgggagtagttggc`, second PCR adding the
  Nextera tails `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` and
  `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, then 10-cycle i5/i7 index PCR) is a
  *separate amplicon validation* library on purified SELECT-seq material, and a
  further gDNA amplicon library genotypes pooled prime-edited cells. Neither is
  the SELECT-seq single-cell library; they must not be merged with it.

🟡 Still open: the index member per well; whether the SNP-specific amplicons
are deliberately retained in the tagmented library or only in the targeted
readout; Smart-seq3 v3 revision details beyond the inspected oligo list.
Deferred originals: Smart-seq3 v3 protocols.io body steps beyond the oligo
list, the Hagemann-Jensen Smart-seq3 paper, NEB M0653 and IDT reporter product
documentation, and the code repository named in Code Availability.


## Publication follow-up (2026-10-30)

The defining original inspected here is explicitly **bioRxiv v1, posted8July2026, not certified by peer review**. Its listing in a catalogue section is not evidence of a peer-reviewed successor. The actual source endpoint is **Element Aviti paired-end150cycles**; it must not be relabelled Illumina merely because some adapters match familiar sequences.

The cached exact cited Smart-seq3 Protocol V.3 prints19February2020 and DOI10.17504/protocols.io.bcq4ivyw. The source positively establishes inherited RT/WTA, SELECT-specific SNP amplification, Cas12a fluorescence and downstream tagmented/indexed DNA library. There is still no source-owned molecular graph or page builder. Constructing that inherited front-end with actual optical reporter and library allocations is substantive model work; this remains a **model-scope hold**, not evidence that the original or its primers are unavailable. Ordinary unread operative V.3 sections and unassigned per-well index values remain separate from unavailable documents. The separate genotyping libraries must not be merged with the transcriptome-library trace.
