# RAMPAGE

## 1. What it is

RAMPAGE selects 5′-complete capped transcripts and uses paired-end sequencing to connect a base-resolution promoter start to downstream transcript structure.

## 2. Why this CAGE-class representative

- 🟢 Batut P et al. *High-fidelity promoter profiling reveals widespread alternative promoter usage and transposon-driven developmental gene expression.* Genome Research (2013). [doi:10.1101/gr.139618.112](https://doi.org/10.1101/gr.139618.112).
- 🟢 Batut P & Gingeras TR. *RAMPAGE: promoter activity profiling by paired-end sequencing of 5′-complete cDNAs.* Current Protocols (2013). [doi:10.1002/0471142727.mb25b11s104](https://doi.org/10.1002/0471142727.mb25b11s104).
- 🟢 RAMPAGE combines template-switching, cap trapping and exonuclease depletion with paired-end Illumina sequencing.
- 🟡 RAMPAGE is represented rather than a “CAGE family” page because it has a protocol-specific defining paper, a complete open molecular protocol, and a distinct paired-end construct.
- 🟢 The open protocol gives the exact TSO, RT primer, final PCR primers and custom Read 1/Read 2 sequencing primers. The TSO carries a six-base inline sample barcode; its identity varies across the published barcode set.
- 🟢 The detailed protocol says there are 40 six-base barcode variants with pairwise Hamming distance at least three, printed in Supplementary File 1. The accessible PMC author-manuscript XML gives the TSO scaffold but not those 40 assignments; the archived supplementary PDF is behind PMC's download gate in the source manifest.
- 🟢 A later [Gingeras Lab ENCODE library-production document](https://www.encodeproject.org/documents/db4fd5dc-145b-44a7-b3ea-0cddb5b458cb/@@download/attachment/296007.pdf), PDF pp. 9–10, prints four sets of five six-base RAMPAGE indexes (20 strings) for a particular tissue library. It also repeats the general 40-barcode description, but this partial list cannot establish all 40 original Supplementary File 1 members or which were used in the defining study. Its printed table labels the sets 1, 2, 3 and 5, while the reaction text calls for set 4; that ENCODE sample's barcode assignment is therefore unresolved. 🔴 The original supplement remains inaccessible; a direct PMC source-bucket probe returned HTTP 404.
- 🔴 [Philippe Batut's author thesis](https://repository.cshl.edu/31200/1/Batut%20Thesis.pdf) is search-indexed as reproducing RAMPAGE Table 25B.11.1 with barcode numbers 01–40. Direct thesis retrieval timed out, so its sequences and correspondence to the protocol supplement remain unchecked. The exact 40-member family and study-specific assignments remain blocked pending an original-source read.

## 3. Molecular path

Template switching creates a 5′-completion handle, while cap trapping independently selects capped RNA–cDNA hybrids. The insert end next to that handle marks the transcription start site.
