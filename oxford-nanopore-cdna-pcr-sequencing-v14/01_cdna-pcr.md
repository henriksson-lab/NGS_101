# Oxford Nanopore cDNA-PCR Sequencing V14 (SQK-PCS114)

## 1. What it is

SQK-PCS114 converts RNA to full-length strand-switched cDNA, incorporates a UMI, PCR-amplifies molecules with rapid-attachment tags, and attaches a motor-loaded Rapid Adapter for nanopore sequencing.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *cDNA-PCR Sequencing V14 (SQK-PCS114)*, current vendor protocol: https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114
- 🟢 Oxford Nanopore, *Chemistry Technical Document*, CHTD_500 v1 rev AV (24 July 2026), section 15, publishes CRTA, RTP, SSPII, both cDNA Primer sequences, and the Rapid Adapter/RAP T strand.
- 🟢 The protocol explicitly divides library preparation into reverse transcription and strand switching, full-length selection by PCR, and Rapid Adapter addition.
- 🟢 The current protocol states that strand switching incorporates a UMI and PCR primers add terminal tags used for rapid-adapter attachment.
- 🟢 SSPII is published with a partially degenerate `V/T` UMI region and terminal `mGmGmG`; the modification spelling is retained verbatim in the ordering list.
- 🔴 The complete motor-bearing Rapid Adapter assembly, motor attachment, and attachment chemistry are not disclosed. The published RA strand must not be mistaken for that entire assembly.

## 3. Molecular path

This protocol is kept separate from Direct cDNA because PCR is both the enrichment/amplification step and the means of installing the Rapid Adapter attachment boundary.
