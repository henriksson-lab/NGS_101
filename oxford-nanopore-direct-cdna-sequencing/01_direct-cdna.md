# Oxford Nanopore Direct cDNA Sequencing (SQK-DCS109)

## 1. What it is

The historical SQK-DCS109 workflow makes strand-switched duplex cDNA without PCR, repairs and dA-tails it, then ligates a motor-loaded dT-overhang sequencing adapter.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *Chemistry Technical Document*, Direct cDNA Sequencing Kit section: https://nanoporetech.com/document/chemistry-technical-document
- 🟢 The vendor describes RT primer annealing, strand-switch RT, RNA removal, second-strand synthesis, end repair/dA tailing and ligation of a dT-tailed sequencing adapter.
- 🟢 The workflow is PCR-free and starts from substantially more poly(A) RNA than cDNA-PCR.
- 🔴 Kit oligo, sequencing-adapter and motor sequences are proprietary.
- 🟢 SQK-DCS109 is historical rather than a current V14 kit; its status remains visible instead of conflating it with PCS114.

## 3. Molecular path

The cDNA duplex is treated like a genomic-DNA ligation library after second-strand synthesis. The A:T overhang junction, rather than a PCR-added rapid tag, installs the sequencing adapter.
