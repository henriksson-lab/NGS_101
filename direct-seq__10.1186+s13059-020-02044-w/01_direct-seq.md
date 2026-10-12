# Direct-seq — source notes

Defining paper: Song et al., *Genome Biology* (2020), [doi:10.1186/s13059-020-02044-w](https://doi.org/10.1186/s13059-020-02044-w).

- 🟢 PMC full text, figures and supplementary oligo tables were read.
- 🟢 A 30-nt mixed A/G “8A8G” tract is inserted into the sgRNA scaffold so ordinary poly(dT) primers capture the Pol III guide transcript.
- 🟢 In the 10x 3′ v3 workflow, guide and endogenous RNA acquire the same cell barcode and UMI during RT.
- 🟢 Pre-amplified cDNA is split: one branch makes the normal mRNA library and nested guide-specific PCR enriches a guide-index library before P5/P7 completion.
- 🟢 The exact capture tract and terminal PCR primers are printed in the supplementary tables.

## Programmed transcript and capture

🟢 Supplementary Tables S2–S3 give the complete Tail-8A8G scaffold. The 30-nt `AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA` tract sits between the programmed scaffold and Pol III terminator. In the 10x 3′ implementation a published 72-nt human Gln tRNA precedes the guide spacer and supplies the nested-PCR priming site.

🟡 The 10x poly(dT) bead primer can anneal across the capture tract, but the three encoded G residues are intentional mismatches. The schematic therefore places the exact oligos without claiming a perfect 30-bp duplex. Reverse transcription still gives the guide transcript the normal 16-nt cell barcode and 12-nt UMI.

## Nested guide-library PCR

🟢 The first enrichment PCR uses `P5_Read1` and `tRNA_Read2`. The 3′ 23 nt of `tRNA_Read2` exactly match bases 25–47 of the human Gln tRNA; its 5′ 34 nt are the TruSeq Read 2 primer site. A third PCR retains `P5_Read1` and replaces the reverse primer with one of two printed P7/index/Read-2 primers.

🟡 The final construct therefore makes its read structure computable: Read 1 is the ordinary 16-nt cell barcode plus 12-nt UMI. With a 90-cycle Read 2, cycles 1–23 are the tRNA primer site, 24–48 the remaining tRNA sequence, 49–68 the guide spacer and 69 onward the programmed scaffold. The guide spacer—not the 8A8G tract—is the perturbation identifier.
