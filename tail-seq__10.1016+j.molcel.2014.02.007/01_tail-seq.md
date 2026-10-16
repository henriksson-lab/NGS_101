# TAIL-seq — research notes

Defining source: Chang et al., 2014, [doi:10.1016/j.molcel.2014.02.007](https://doi.org/10.1016/j.molcel.2014.02.007).

## Sources read

- 🟢 Defining article's full Methods in an [article mirror](https://phmpt.org/wp-content/uploads/2024/06/125752_S2_M3_33_chang-2014.pdf), especially library preparation and sequencing. The article points to Supplemental Experimental Procedures for reaction details; those remain unavailable in the local source bundle.
- 🟢 Later implementations were read only to confirm the role of biotin capture and partial RNase T1 digestion; their altered splints and oligos were not imported into the original protocol.
- 🔴 `get_sources.py` records the original *Molecular Cell* full text and supplement as manual retrieval, while the external mirror provides the main article. No source-verified original oligo list is present locally. Later mTAIL-seq and aTAIL-seq sources cannot establish the 2014 adapter termini or modifications.

## Molecular path

🟢 A biotinylated adapter is ligated to native RNA 3-prime ends. Partial RNase T1 digestion shortens the upstream RNA, and streptavidin capture selects fragments that still carry the native terminus. A 5-prime adapter, RT and PCR form the library.

🟢 Read 1 identifies the transcript; the long paired read traverses terminal additions and poly(A). Tail length is determined from raw fluorescence because standard base calling is unreliable through long homopolymers.

🟡 Adapter bases not verified in the defining source remain inferred placeholders.

🔴 The exact biotinylated 3′ RNA adapter, 5′ RNA adapter, adapter-complementary RT primer and PCR primer sequences remain open. The full article Methods describe their roles but do not print their bases; the linked Supplemental Experimental Procedures have not yet been checked. Raw FASTQ could potentially reveal ligated adapter sequence, but it cannot prove a biotin or ligation-end chemical modification; no such inference has been made here.
