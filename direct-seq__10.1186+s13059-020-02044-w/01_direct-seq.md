# Direct-seq — source notes

Defining paper: Song et al., *Genome Biology* (2020), [doi:10.1186/s13059-020-02044-w](https://doi.org/10.1186/s13059-020-02044-w).

- 🟢 PMC full text, figures and supplementary oligo tables were read.
- 🟢 A 30-nt mixed A/G “8A8G” tract is inserted into the sgRNA scaffold so ordinary poly(dT) primers capture the Pol III guide transcript.
- 🟢 In the 10x 3′ v3 workflow, guide and endogenous RNA acquire the same cell barcode and UMI during RT.
- 🟢 Pre-amplified cDNA is split: one branch makes the normal mRNA library and nested guide-specific PCR enriches a guide-index library before P5/P7 completion.
- 🟢 The exact capture tract and terminal PCR primers are printed in the supplementary tables.
