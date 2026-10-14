# TAIL-seq — research notes

Defining source: Chang et al., 2014, [doi:10.1016/j.molcel.2014.02.007](https://doi.org/10.1016/j.molcel.2014.02.007).

## Sources read

- 🟢 Defining article's publisher abstract and workflow figure, which establish the assay and main read roles; full Methods were not available in the source fetch.
- 🟢 Later implementations were read only to confirm the role of biotin capture and partial RNase T1 digestion; their altered splints and oligos were not imported into the original protocol.
- 🔴 `get_sources.py` records the original *Molecular Cell* full text and supplement as manual retrieval; no source-verified original oligo list is present locally. Later mTAIL-seq and aTAIL-seq sources cannot establish the 2014 adapter termini or modifications.

## Molecular path

🟢 A biotinylated adapter is ligated to native RNA 3-prime ends. Partial RNase T1 digestion shortens the upstream RNA, and streptavidin capture selects fragments that still carry the native terminus. A 5-prime adapter, RT and PCR form the library.

🟢 Read 1 identifies the transcript; the long paired read traverses terminal additions and poly(A). Tail length is determined from raw fluorescence because standard base calling is unreliable through long homopolymers.

🟡 Adapter bases not verified in the defining source remain inferred placeholders.

🔴 The exact biotinylated 3′ RNA adapter, 5′ RNA adapter, adapter-complementary RT primer and PCR primer sequences remain open. The available original article record describes their roles but does not print their bases. Raw FASTQ could potentially reveal ligated adapter sequence, but it cannot prove a biotin or ligation-end chemical modification; no such inference has been made here.
