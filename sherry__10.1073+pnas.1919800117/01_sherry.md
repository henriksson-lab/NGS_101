# SHERRY — source notes

Defining paper: Di et al., *PNAS* 2020, [doi:10.1073/pnas.1919800117](https://doi.org/10.1073/pnas.1919800117).

## Sources read

- 🟢 PMC7022195 full-text XML, especially Results and Methods “SHERRY Library Preparation and Sequencing”.
- 🟢 PNAS supplementary methods, pp. 2–4, including the transposome oligos and hybrid-tagmentation conditions.

## Construct-changing steps

- 🟢 Poly(A) RNA is copied with poly(T)30VN and SuperScript II; the ISPCR tail and TSO used by Smart-seq2 are omitted.
- 🟢 Tn5 is loaded with `CTGTCTCTTATACACATCT` plus either `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` or `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` and directly tagments the RNA/cDNA hybrid.
- 🟢 These are the unmodified standard SHERRY adaptor sequences printed in the
  supplementary methods. The supplementary Fig. S4 experiments separately make
  `rCrArG` or terminal `rG` substitutions at the 3′ adaptor ends; those experimental
  variants are not silently applied to the standard library construct.
- 🟢 Extension uses SuperScript II for 10/200-ng input or Bst 2.0 WarmStart for 100-pg input, followed by indexed common-primer PCR.
- 🟢 The paper's strand tests show that tagmented cDNA contributes most of the final library, although adapter transfer to both RNA and DNA was demonstrated.
- 🟡 The final duplex is drawn as the heterologous S5/S7 amplifiable Nextera product; homotypic Tn5 products are not sequencing-library endpoints.

The paper sequenced on NextSeq 500 or HiSeq 4000. It does not report a UMI or cellular barcode in standard SHERRY.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC7022195 JATS175–180, Library preparation and sequencing.

Specified input-dependent gap-extension mixture (SuperScript II/Q5 versus Bst2.0 WarmStart) and 12/15/25 PCR cycles for 200ng/10ng/100pg RNA. Candidate non-amplifiable Tn5 products remain unresolved rather than presumed physically discarded.

Deferred: Published supplementary oligos, transfer strand geometry and actual gap-fill products; final 1:1 bead cleanup and read acquisition. This pass does not establish a complete source review or publication-wide absence.
