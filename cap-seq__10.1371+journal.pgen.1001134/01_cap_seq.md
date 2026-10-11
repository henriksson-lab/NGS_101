# CAP-seq — CXXC affinity purification sequencing

## Name boundary

This page covers **CXXC Affinity Purification plus sequencing**, the named DNA method introduced by Illingworth and colleagues. It is not the similarly named capped-RNA “Cap-seq”/CIP–TAP workflow; that is a distinct protocol and should receive its own page if added.

## Sources read

- 🟢 Illingworth *et al.*, “Orphan CpG islands identify numerous conserved promoters in the mammalian genome,” *PLOS Genetics* (2010), DOI `10.1371/journal.pgen.1001134`, PMCID `PMC2944787` (defining paper and methods).
- 🟢 Earlier MAP/CAP preparation methods cited by the defining paper were used only to interpret the named chromatography steps; no unprinted oligo was imported from them.

## Molecular interpretation

- 🟢 CAP enriches genomic fragments containing clusters of **unmethylated CpGs** by binding them to recombinant CXXC immobilized on nickel-charged Sepharose.
- 🟢 The sequencing adapters are ligated to input DNA **before** affinity purification. This matters: selection acts on an already sequenceable library rather than on naked DNA followed by a separate library build.
- 🟢 Binding occurs in column buffer containing 0.1 M NaCl. Human libraries are washed at 600 mM NaCl; the optimized mouse workflow uses 560 mM because mouse CpG islands are, on average, less CpG-dense. Bound material is eluted at 1 M NaCl.
- 🟢 The enriched libraries were sequenced single-end on the Solexa platform.
- 🔴 The defining paper refers to earlier Solexa library methods but does not print the adapter or sequencing-primer bases, so the diagram labels these ends by role.

## Final readout

A mapped read marks an adapter-ligated genomic fragment retained because its clustered unmethylated CpGs bound CXXC under the chosen salt conditions. CAP-seq reports affinity enrichment, not a base-by-base conversion call.
