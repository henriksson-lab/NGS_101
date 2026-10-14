# m6A-seq / MeRIP-seq — research notes

Defining source: Dominissini et al., 2012, [doi:10.1038/nature11112](https://doi.org/10.1038/nature11112).

## Sources read

- 🟢 Author-hosted PDF of the defining *Nature* paper, Methods (“RNA preparation” and “m6A-seq”): poly(A) enrichment, RNA fragmentation, anti-m6A immunoprecipitation, competitive elution and input/IP sequencing. <https://www.weizmann.ac.il/molgen/Sorek/files/publications/Dominissini_Nature_2012.pdf>
- 🟢 Publisher Supplementary Figures PDF; its figures do not supply a library oligo table.

## Boundary

🟢 m6A-seq/MeRIP-seq enriches a population of RNA fragments; the modification does not become a barcode or a deterministic base change. Peak position is consequently broader than single-nucleotide resolution.

🟡 The defining article does not make one sequence-level adapter design intrinsic to the assay, so the final conventional Illumina library arms remain inferred.

🟢 The original Methods specify library generation with an **Illumina mRNA sequencing kit** and sequencing on GAIIx with genomic sequencing kit V2. They do not provide a kit catalog number or the RT primer, adapter or PCR oligo sequences. The later author protocol, [Dominissini et al. 2013](https://doi.org/10.1038/nprot.2012.148), describes a detailed implementation but should not silently determine the 2012 library construct.

🔴 The exact 2012 kit generation and its adapter/primer synthesis specification remain unresolved. This is a kit provenance gap rather than evidence for a new m6A-specific adapter.
