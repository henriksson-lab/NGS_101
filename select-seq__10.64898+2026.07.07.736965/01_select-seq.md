# SELECT-seq — source intake

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
