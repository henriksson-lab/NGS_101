# CEL-Seq2 — source note

The detailed source audit, oligo transcriptions, protocol steps, and comparison with the
original CEL-Seq method are retained in the companion
[CEL-Seq / CEL-Seq2 research note](../cel-seq__10.1016+j.celrep.2012.08.003/01_cel-seq.md).

**Defining source:** Hashimshony T *et al.* “CEL-Seq2: sensitive highly-multiplexed
single-cell RNA-Seq.” *Genome Biology* 17, 77 (2016).
doi:[10.1186/s13059-016-0938-8](https://doi.org/10.1186/s13059-016-0938-8).

The schematic uses the CEL-Seq2 RT-primer and randomhexRT sequences printed in
Supplementary Table S2 (Additional file 5, MOESM5) and the 96-member primer list in
the detailed protocol (Additional file 4, MOESM4). Additional file 6 (MOESM6) is the
separate Fluidigm C1 workflow, not an oligo table; this schematic follows the manual
well workflow. The detailed protocol calls for Illumina
RP1 and indexed RPIX, directing readers to Illumina for their sequences. The
[official Illumina TruSeq Small RNA adapter table](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm)
(document 1000000002694 v22) prints both oligos and the RPI index table. The model uses
RP1 and representative RPI1; the actual indexed RPIX member used in each experiment is
not specified. Their 3′ ends land on 16 and 21 known bases, respectively,
of the selected CEL-Seq2 cDNA template. The transcript insert and cell UMI/barcode are
variable; sequencing-primer placement remains unverified.
