# SHERRY — source notes

Defining paper: Di et al., *PNAS* 2020, [doi:10.1073/pnas.1919800117](https://doi.org/10.1073/pnas.1919800117).

## Sources read

- 🟢 PMC7022195 full-text XML, especially Results and Methods “SHERRY Library Preparation and Sequencing”.
- 🟢 PNAS supplementary methods, pp. 2–4, including the transposome oligos and hybrid-tagmentation conditions.

## Construct-changing steps

- 🟢 Poly(A) RNA is copied with poly(T)30VN and SuperScript II; the ISPCR tail and TSO used by Smart-seq2 are omitted.
- 🟢 Tn5 is loaded with `CTGTCTCTTATACACATCT` plus either `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` or `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` and directly tagments the RNA/cDNA hybrid.
- 🟢 Extension uses SuperScript II for 10/200-ng input or Bst 2.0 WarmStart for 100-pg input, followed by indexed common-primer PCR.
- 🟢 The paper's strand tests show that tagmented cDNA contributes most of the final library, although adapter transfer to both RNA and DNA was demonstrated.
- 🟡 The final duplex is drawn as the heterologous S5/S7 amplifiable Nextera product; homotypic Tn5 products are not sequencing-library endpoints.

The paper sequenced on NextSeq 500 or HiSeq 4000. It does not report a UMI or cellular barcode in standard SHERRY.
