# Oxford Nanopore cDNA-PCR Sequencing V14 (SQK-PCS114)

## 1. What it is

SQK-PCS114 converts RNA to full-length strand-switched cDNA, incorporates a UMI, PCR-amplifies molecules with rapid-attachment tags, and attaches a motor-loaded Rapid Adapter for nanopore sequencing.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *cDNA-PCR Sequencing V14 (SQK-PCS114)*, current vendor protocol: https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114
- 🟢 The protocol explicitly divides library preparation into reverse transcription and strand switching, full-length selection by PCR, and Rapid Adapter addition.
- 🟢 The current protocol states that strand switching incorporates a UMI and PCR primers add terminal tags used for rapid-adapter attachment.
- 🔴 CRTA, RTP, SSPII, cPRM and Rapid Adapter molecular sequences are proprietary.

## 3. Molecular path

This protocol is kept separate from Direct cDNA because PCR is both the enrichment/amplification step and the means of installing the Rapid Adapter attachment boundary.
