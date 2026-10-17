# INSIGHT — research notes

Defining paper: Wu *et al.* (2021), DOI [10.1126/sciadv.abe5054](https://doi.org/10.1126/sciadv.abe5054). The published PMC7880595 full text was read: Materials and Methods lines 155 (beacon annealing), Table 1 reagents lines 190-229, Table 2 oligonucleotides lines 235-266, Steps 1-2 lines 268-365 (lysis, NASBA reaction mixtures, 41 °C incubation), Step 3 lines 369-406 (pooling, column purification, indexing PCR, Table 5 conditions, size selection, MiSeq run) and Step 4 line 410 (read trimming and the 102-nt reference).

## Reaction detail (round 2026-10-17-breadth-first-11)

- 🟢 Stage 1 NASBA runs at 41 °C for 90–120 min after a 65 °C/2 min or 95 °C/5 min pre-incubation and 10 min at 41 °C; primers are at 25 nM each in the reaction (Tables 3 and 4).
- 🟢 Table 1 lists only ProtoScript II reverse transcriptase, RNase H and T7 RNA polymerase as NASBA enzymes. 🟡 Model correction: the graph previously declared a separate "NASBA DNA polymerase" for second-strand synthesis; that input was removed and the reverse transcriptase now carries the step.
- 🟢 Stage 2: 2 µl of each end product are pooled, column-purified on a QIAquick PCR purification column, then amplified with the P5/P7 end primers in KAPA HiFi HotStart ReadyMix for 15 cycles (95 °C 3 min; 98 °C 20 s / 60 °C 15 s / 72 °C 30 s; 72 °C 4 min), with long elongation and minimal cycles to reduce barcode hopping (Step 3, Table 5).
- 🟢 AMPure XP double size selection (0.55× then 0.75×) follows the PCR; 30% PhiX is added before sequencing. 🟡 Model correction: the graph now shows that selection with its rejected off-size fraction.
- 🟢 Read layout: 150-bp paired-end MiSeq, MiSeq Reagent Kit v2 (300 cycles). Analysis trims both reads to the first 80 nt and merges them (Step 4).
- 🟢 The merged reference is 102 nt: `NNNNN` + the 92-nt S-gene region + `NNNNN`, the first 5 nt being the right barcode and the last 5 the reverse complement of the left barcode. 🟡 Model correction: the drawn amplicon was a 46-nt placeholder and now carries those 92 published bases; the selftest pins the transcription and checks that the amplicon ends at the forward-primer site and starts at the reverse-primer site.

## Published discrepancy

🟢 Table 2's barcoded forward primer handle is written `tgactggagttcagacgtgtgctcttccgatct` (33 nt), one base shorter at the 5′ end than the 34-nt TruSeq Read 2 arm inside the same table's P7 end primer (`...GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`). The drawn library carries the full-length arm, which the P7 end primer supplies during the pooled PCR; the NASBA primer is listed verbatim. The page shows this.

## Evidence

- 🟢 Table 2 gives the P8 target primers, 5-nt sample barcodes, T7 promoter, Illumina handle, outer indexed PCR primers and custom T7-containing sequencing primer.
- 🟢 NASBA uses reverse transcriptase, RNase H and T7 RNA polymerase at 41 °C.
- 🟢 Completed reactions are pooled only after the sample barcode pair is installed; residual primers are removed before a minimal-cycle universal PCR.
- 🟢 The same Table 2 prints the alternative readout oligos, which do not enter the sequenced library: a 2′-O-methyl toehold beacon `FAM-AUUGACAGUCUACUAAUUUGGUUAAAAACAAAUGUGUCAA-BHQ1dT-UUCAACUUCAAUG-propyl`, PCRD capture probe A `FAM-AAAAGUCUACUAAUUUGGUUAAAA`, and probe B `ACAAAUGUGUCAAUUUCAACUUCA-Biotin`. These are RNA sequences as written; the beacon's internal quencher and terminal propyl are retained verbatim.

The lateral-flow and molecular-beacon readouts are useful parallel endpoints, but the schematic concentrates on the molecules that enter NGS.
