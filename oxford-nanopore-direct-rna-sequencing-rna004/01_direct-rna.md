# Oxford Nanopore Direct RNA Sequencing (SQK-RNA004)

## 1. What it is

SQK-RNA004 prepares native polyadenylated RNA for nanopore sequencing without converting the sequenced molecule into amplified DNA.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *Direct RNA Sequencing Kit SQK-RNA004* protocol and current chemistry technical document.
- 🟢 Garalde DR et al. *Highly parallel direct RNA sequencing on an array of nanopores.* Nature Methods 15:201–206 (2018). [doi:10.1038/nmeth.4577](https://doi.org/10.1038/nmeth.4577).
- 🟢 The current workflow ligates a reverse-transcription adapter to poly(A) RNA, synthesizes a complementary stabilizing strand, and ligates a motor-loaded sequencing adapter.
- 🟢 Oxford Nanopore states that only RNA, not the reverse-transcribed strand, is sequenced.
- 🔴 Current adapter bases, motor identity and attachment sites are proprietary.

## 3. Molecular path

The poly(A)-proximal RNA 3′ end is adapted first. It is therefore the end presented to the pore; the native RNA travels through the pore from 3′ toward 5′ while the cDNA strand serves as structural support.
