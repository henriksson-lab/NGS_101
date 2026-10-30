# Apt-seq

Defining original: **Combined aptamer and transcriptome sequencing of single cells**, DOI[10.1038/s41598-018-21153-y](https://doi.org/10.1038/s41598-018-21153-y).

Capped source/model pass,2026-10-29. Private and incomplete;structural graph validity is separate from source completeness. No missing-information claim has passed the two-reader gate.

## Originals inspected

- `s41598-018-21153-y_PMC5811598.xml`: Published article Methods85–99 and Results/Fig1: aptamer folding/staining, five washes, source Drop-seq capture, Maxima DNA/RNA copying classes, 13cycle joint PCR, aptamer-specific Nextera overhangs and MiSeqv3 paired-end acquisition.
- `s41598-018-21153-y_supp_41598_2018_21153_MOESM1_ESM.pdf`: Published SI Tables1–2, originalpp3–4: TD05 sequence and flanks, printed B followed by21A, aptamer-specific Nextera-tailed primers, smartPCR primer and mixed DNA/rGrGrG TSO. Rendered originalp3 confirms printedB and bold identity region.

## Chemistry and scope

🟢 The original Apt-seq uses DNA affinity aptamers whose own sequences identify the binders. The local sequence example is TD05; it does not stand for the other seven panel members or guarantee specific binding. Aptamers are folded in PBS/Mg/glucose at95°C5min then ice10min, pooled, incubated with mixed Ramos/3T3 cells30min on ice, and washed through five cold centrifugation/resuspension cycles. These are physical cell/wash separations, not computational enrichment.

🟢 The source uses Drop-seq cell/bead droplets and polyT capture. MaximaHMinus produces two distinct classes: cDNA–RNA hybrids from RNA and barcoded double-stranded aptamer DNA from DNA templates. Calling both ordinary RNA reverse transcription would erase the aptamer's physical identity. Full historical bead polymer, cell-barcode/UMI widths and workup details remain unassembled; no current Drop-seq or inDrops recipe is substituted.

🟢 Source joint PCR is13cycles, with aptamer-specific primers at1:20 relative to the smartPCR primer. The aptamer primers append `GTCTCGTGGGCTCGG`; short aptamers are expected by the source to experience little Nextera tagmentation, but this is not enforced as perfect exclusion. The mixed sample receives tagmentation; no unreported physical RNA/aptamer separation is invented. Further processing follows the cited Drop-seq original/protocol. MiSeqReagentKitv3 paired-end acquisition is reported; cycles, sample-index identities and exact read-primer geometry remain unassembled. In-silico aptamer/transcript separation is not physical library fractionation.

🟢 Published TD05: `AAGGAGCAGCGTGGAGGATA AACACCGTGGAGGATAGTTCGGTGGCTGTTCAGGGTCTCCTCCCGGTG TTAGGGTGTGTCGTCGTGGT BAAAAAAAAAAAAAAAAAAAAA`. The original table renders B immediately before the21A tail. It is retained as the printed constrained nucleotide class (C/G/T), not interpreted as a biotin modification or silently replaced by a fixed base. The recognition sequence carries feature identity; neither that identity nor the tail is a per-molecule UMI.

🟢 Published TC/TD/TE Nextera primer: `GTCTCGTGGGCTCGGAAGGAGCAGCGTGGAGGATA`;smartPCR `AAGCAGTGGTATCAACGCAGAGT`;TSO `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`. The final three TSO guanosines have RNA sugars; they are not silently rendered as an all-deoxy oligo. Aptamer folded-state assertions end after copying/amplification, not carried as a predicted PCR-product structure.

Component ancestry: Methods99 explicitly follows Drop-seq (ref13) and its online protocol(ref70) with listed modifications. That supports the bead/capture/library operational relationship only; whole later recipes remain unread. The source mentions inDrops software adaptation for analysis; this does not establish inDrops molecular preparation inheritance. The2025 JACS Apt-seq homonym is not used as authority or a second counted entry.

## Deferred originals and questions

- Direct cited Drop-seq defining original and online protocol ref70; full historical workflow and cleanup/index steps
- Complete bead, cDNA/aptamer product, end geometry and sequencing primer/run declarations
- Remaining source-linked chemistry originals and independent review

## Validation

Local source transcription and physical branch/state checks are recorded in the source-review JSON. Full chemical products,library completion and complete source coverage are not established.
