# microSPLiT — source note

The complete source audit and oligo transcription are retained in the companion
[SPLiT-seq / microSPLiT research note](../split-seq__10.1126+science.aam8999/01_split-seq.md).

**Defining source:** Kuchina A *et al.* “Microbial single-cell RNA sequencing by
split-pool barcoding.” *Science* 371, eaba5257 (2021).
doi:[10.1126/science.aba5257](https://doi.org/10.1126/science.aba5257).

The schematic uses the authors' Supplementary Table S3, available both with the preprint
and through the later detailed protocol. It supplies the complete split-pool oligos,
BC_0062 and BC_0108 cDNA-amplification primers, BC_0243/BC_0244 adapter duplex, and
BC_0027 plus indexed BC_0076–BC_0083 terminal PCR pairs. The page instantiates BC_0076
(TSBC07) as a representative exact final library.

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
The 2024 procedure confirms this later implementation, while the unavailable
2021 published Supplementary Table S3 still prevents an exact version-by-version
comparison 🔴.
