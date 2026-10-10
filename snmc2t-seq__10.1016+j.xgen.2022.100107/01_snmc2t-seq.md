# snmC2T-seq — source notes

Published defining study: Luo et al., *Cell Genomics* (2022), [doi:10.1016/j.xgen.2022.100107](https://doi.org/10.1016/j.xgen.2022.100107). The published paper names the extended assay snmCAT-seq; snmC2T is its transcriptome-plus-methylome core.

- 🟢 PMC full text and STAR Methods were read.
- 🟢 Smart-seq or Smart-seq2 cDNA synthesis substitutes 5-methyl-dCTP for dCTP, producing fully methylated cDNA.
- 🟢 After bisulfite treatment, RNA-derived reads remain highly methylated while most genomic cytosines convert; read-level non-CG methylation partitions the two modalities.
- 🟢 Both molecule types are carried through one snmC-seq2 library, not physically separated libraries.
- 🟢 snmCAT-seq adds M.CviPI GpC marking before sorting; that optional extension is distinguished from the snmC2T core on the page.
- 🔴 Adaptase-added bases are proprietary/unprinted and remain inferred.
