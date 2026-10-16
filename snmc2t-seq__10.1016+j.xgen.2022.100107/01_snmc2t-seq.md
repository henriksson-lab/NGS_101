# snmC2T-seq — source notes

Published defining study: Luo et al., *Cell Genomics* (2022), [doi:10.1016/j.xgen.2022.100107](https://doi.org/10.1016/j.xgen.2022.100107). The published paper names the extended assay snmCAT-seq; snmC2T is its transcriptome-plus-methylome core.

- 🟢 The published PMC full text, JATS, STAR Methods, Key Resources table, supplementary files and associated bioRxiv package were read. The two cited snmC-seq2 source supplements were also checked for the inherited library finish.
- 🟢 Smart-seq or Smart-seq2 cDNA synthesis substitutes 5-methyl-dCTP for dCTP, producing fully methylated cDNA.
- 🟢 After bisulfite treatment, RNA-derived reads remain highly methylated while most genomic cytosines convert; read-level non-CG methylation partitions the two modalities.
- 🟢 Both molecule types are carried through one snmC-seq2 library, not physically separated libraries.
- 🟢 snmCAT-seq adds M.CviPI GpC marking before sorting; that optional extension is distinguished from the snmC2T core on the page.
- 🟢 The Key Resources table prints dT30VN_4, N6_2, TSO_3 and ISPCR23_2 verbatim, including `/5SpC3/`, riboguanosines and the locked terminal guanosine. These are now listed on the schematic instead of being treated as unknown Smart-seq oligos.
- 🟢 The inherited snmC-seq2 finish publishes all eight inline-barcode/N9 random primers and 768 paired P5/i5 and P7/i7 primers. The page draws representative pair 1 and computes its PCR landing sites from the shared snmC-seq2 endpoint model.
- 🔴 Adaptase-added bases are proprietary/unprinted and remain the only inferred bases.

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC9004682 JATS540–554; newly retrieved through existing source fetcher.

Added nuclear-branch cDNA amplification (15 cycles cultured nuclei, 14 brain), UDG treatment with source-specific EndoVIII distinction, and paired150-base acquisition. Named SuperScriptII and distinguished whole-cell snmC-seq finish from nuclear snmC-seq2. Source primers verified in published Methods; superseded preprint not used.

Deferred: Full nuclear methylome finish and recovery, selected index acquisition, methylation propagation through amplification and exact Adaptase bases; endonuclease branch remains summarized. This capped pass does not establish a complete source review or publication-wide absence.
