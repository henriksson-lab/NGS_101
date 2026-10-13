# SurePlex with VeriSeq PGS

## Defining sources

🟢 Illumina, *SurePlex Summary Protocol Reference Guide* 15053626 v02, and *VeriSeq PGS Library Prep Reference Guide* 15052877 v05.

## Molecular workflow

🟢 SurePlex takes one cell in 2.5 µl, performs cell extraction, preamplification and whole-genome amplification. The commercial guide names the buffers and enzymes but does not publish the preamplification-primer bases.

🟢 The associated VeriSeq PGS workflow uses 1 ng unpurified SurePlex dsDNA. An engineered transposome fragments and tags it, neutralization stops tagmentation, and limited-cycle PCR uses S503/S504 and N701–N712 primers to make single- or dual-index libraries.

🟢 The intended MiSeq acquisition is single-read: Read 1 is analysed, index reads demultiplex samples, and Read 2 is not acquired.

🔴 SurePlex and PicoPLEX share a commercial lineage, but they are kept as separate named protocols because PicoPLEX DNA-seq integrates dual indexing whereas SurePlex hands its WGA product to the separate VeriSeq PGS library kit.
