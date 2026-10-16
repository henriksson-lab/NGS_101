# microSPLiT — source note

The complete source audit and oligo transcription are retained in the companion
[SPLiT-seq / microSPLiT research note](../split-seq__10.1126+science.aam8999/01_split-seq.md).

**Defining source:** Kuchina A *et al.* “Microbial single-cell RNA sequencing by
split-pool barcoding.” *Science* 371, eaba5257 (2021).
doi:[10.1126/science.aba5257](https://doi.org/10.1126/science.aba5257).

The schematic uses a mirrored copy of the authors' 2021 *Science* Supplementary
Table S3 (`scg_aba5257_table_s3.xlsx`, sheet `General Oligonucleotides`,
lines 1–19). The workbook identifies itself as Table S3; its embedded metadata
names Anna Kuchina as the last editor. The original publisher attachment is still
gated, so byte identity of this mirror has not been verified. The table supplies
the split-pool oligos, BC_0243/BC_0244 adapter duplex, and BC_0027 plus indexed
BC_0076–BC_0083 terminal PCR pairs. The page instantiates BC_0076 (TSBC07)
as a representative final-library index. BC_0062 and BC_0108 cDNA-amplification
primers are printed in the later 2024 detailed protocol, not in the mirrored 2021
Table S3.

**Adapter source-version difference.** The mirrored 2021 Table S3, General
Oligonucleotides row BC_0243 (line 15), prints
`ACACTCTTTCCCTACACGACGCTCTTCCGATC*T` with the adjacent note
`*phosphorothioated base` 🟢. Its BC_0244 opposite strand (line 16) is
`GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` 🟢; the first 32 bases pair, leaving
BC_0243's terminal T as the 3′ overhang 🟡. The 2019 bioRxiv v2 Auxiliary
Table (`869248_media-2.xlsx.txt`, lines 14–15) and 2024 *Nature Protocols*
Table 3 (`...MOESM2_ESM.xlsx.txt`, lines 121–122) print the same BC_0243 bases
without the asterisk 🟢. The molecular model uses the 2021 Table S3 terminal
C*T phosphorothioate annotation. The reason for the notation difference is
unknown; the official 2021 attachment remains inaccessible for direct comparison.

**Later primary-source cross-check (2024).** The authors' detailed Nature Protocols
procedure, doi:[10.1038/s41596-024-01007-w](https://doi.org/10.1038/s41596-024-01007-w),
explicitly describes fragmentation followed by ligation of the BC_0243/BC_0244
adapter duplex and PCR with BC_0027 plus an indexed reverse primer (Figure 2 and
Part 2, steps 134–157). Its Supplementary Table 3 (`...MOESM2_ESM.xlsx`, sheet
`Table 3 General Oligonucleotide`) prints the same representative BC_0027,
BC_0062, BC_0108, BC_0076, BC_0127 (`rGrG+G`), BC_0243 and BC_0244 sequences
used in the model 🟢. It additionally gives **24** indexed final-PCR primers:
BC_0076–BC_0083 and BC_0246–BC_0261 🟢. Thus the 24 sublibraries mentioned
in the earlier source can be assigned distinct published primer sequences;
the model's BC_0076 remains a representative rather than a complete index list.
The 2024 procedure confirms this later implementation. The mirrored 2021 table
allows a targeted sequence comparison, while direct verification against the
publisher-hosted 2021 attachment remains blocked 🔴.
