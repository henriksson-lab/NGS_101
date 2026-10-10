# TAIL-seq — research notes

Defining source: Chang et al., 2014, [doi:10.1016/j.molcel.2014.02.007](https://doi.org/10.1016/j.molcel.2014.02.007).

## Sources read

- 🟢 Defining article, workflow figure and Methods available from the publisher record.
- 🟢 Later implementations were read only to confirm the role of biotin capture and partial RNase T1 digestion; their altered splints and oligos were not imported into the original protocol.

## Molecular path

🟢 A biotinylated adapter is ligated to native RNA 3-prime ends. Partial RNase T1 digestion shortens the upstream RNA, and streptavidin capture selects fragments that still carry the native terminus. A 5-prime adapter, RT and PCR form the library.

🟢 Read 1 identifies the transcript; the long paired read traverses terminal additions and poly(A). Tail length is determined from raw fluorescence because standard base calling is unreliable through long homopolymers.

🟡 Adapter bases not verified in the defining source remain inferred placeholders.
