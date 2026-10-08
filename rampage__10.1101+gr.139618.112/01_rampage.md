# RAMPAGE

## 1. What it is

RAMPAGE selects 5′-complete capped transcripts and uses paired-end sequencing to connect a base-resolution promoter start to downstream transcript structure.

## 2. Why this CAGE-class representative

- 🟢 Batut P et al. *High-fidelity promoter profiling reveals widespread alternative promoter usage and transposon-driven developmental gene expression.* Genome Research (2013). [doi:10.1101/gr.139618.112](https://doi.org/10.1101/gr.139618.112).
- 🟢 Batut P & Gingeras TR. *RAMPAGE: promoter activity profiling by paired-end sequencing of 5′-complete cDNAs.* Current Protocols (2013). [doi:10.1002/0471142727.mb25b11s104](https://doi.org/10.1002/0471142727.mb25b11s104).
- 🟢 RAMPAGE combines template-switching, cap trapping and exonuclease depletion with paired-end Illumina sequencing.
- 🟡 RAMPAGE is represented rather than a “CAGE family” page because it has a protocol-specific defining paper, a complete open molecular protocol, and a distinct paired-end construct.
- 🟢 The open protocol gives the exact TSO, RT primer, final PCR primers and custom Read 1/Read 2 sequencing primers. The TSO carries a six-base inline sample barcode; its identity varies across the published barcode set.

## 3. Molecular path

Template switching creates a 5′-completion handle, while cap trapping independently selects capped RNA–cDNA hybrids. The insert end next to that handle marks the transcription start site.
