# TT-seq

## 1. Protocol boundary

TT-seq is the Schwalb et al. workflow: a five-minute 4-thiouridine pulse, spike-in addition,
RNA fragmentation **before** labelled-RNA purification, DNase treatment and a strand-specific
commercial library preparation. The ordering of fragmentation and selection distinguishes it
from the unfragmented 4sU-seq control in the same experiment.

## 2. Sources read

🟢 Schwalb et al., *TT-seq maps the human transient transcriptome*, *Science* (2016),
DOI `10.1126/science.aad9841`.

🟢 The defining paper's supplementary methods (`aad9841-Schwalb-SM.pdf`), retrieved from
the supplementary-material mirror listed in `ref/MANIFEST.md`. The RNA-seq methods state:
500 µM 4sU for 5 min; spike-in addition during TRIzol extraction; BioRuptor fragmentation at
240 ng/µl for one 30 s ON / 30 s OFF cycle; labelled-RNA purification; 2 U Turbo DNase; NuGEN
Ovation Human Blood RNA-seq kit; HiSeq 1500. The processing section reports paired-end 50-base
reads and an additional six-base barcode read.

## 3. Molecular state changes

🟢 Newly synthesized RNA incorporates 4sU during the pulse. The six synthetic spike-ins include
three made with a 1:10 4sUTP:UTP ratio and three made with UTP only.

🟢 Purified RNA is split. The TT-seq arm is fragmented before the labelled fraction is purified;
therefore a fragment is selected only when it contains a pulse-labelled nucleotide. This removes
the positional bias produced when an entire nascent transcript is captured through one labelled
region.

🟢 The supplement delegates labelled-RNA purification to its reference 26. The established
4sU purification chemistry couples biotin to the thiol, binds streptavidin and elutes labelled RNA
under reducing conditions; these operations are shown as selection chemistry, not as an adapter.

🟢 Captured RNA is DNase-treated and passed to the strand-specific NuGEN kit.

🟡 The supplement does not print the NuGEN adapter oligos. The public page therefore marks its
standard Illumina-compatible final construct and primer landing sites as inferred. The six-base
single-index and paired-end run structure itself is directly reported.

## 4. Read layout

- Read 1: 50 bases of strand-specific cDNA.
- Index 1: six barcode bases.
- Read 2: 50 bases from the opposite cDNA end.
- No Index 2 read is reported for the defining experiment.

## 5. Open evidence boundary

The exact discontinued Ovation Human Blood RNA-seq adapter sequences remain undisclosed here.
They must not be replaced with an unmarked modern TruSeq kit merely because both run on Illumina.
