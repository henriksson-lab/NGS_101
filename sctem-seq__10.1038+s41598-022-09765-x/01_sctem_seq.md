# scTEM-seq

**Defining source:** doi:10.1038/s41598-022-09765-x.

## Molecular boundary

🟢 DNA and mRNA from one cell are separated. Bisulfite-converted DNA is amplified with indexed SINE Alu or LINE-1 primer panels while the RNA branch is processed separately.

## Evidence boundary

🟢 The source publishes 28 SINE Alu primer pairs, internal dual indexes and 0–5 N staggers; they are a panel, not one universal oligo.

## Primer evidence

🟢 Supplementary Table S1 in the published article's MOESM2 workbook is an
ordering-level table of the SINE Alu and LINE-1
PCR primer panels. Its first-generation rows include:

| Primer | 5′→3′ sequence as published |
|---|---|
| Sine.Alu_F | CTACACGACGCTCTTCCGATCTAAAAATATAAAAAATTAGTTGGG |
| Sine.Alu_R | CAGACGTGTGCTCTTCCGATCTNNNNNNNNCCCAAACTAAAATACAATAAC |
| LINE.L1_F | CTACACGACGCTCTTCCGATCTTTTTGAGTTAGGTGTGGGATATA |
| LINE.L1_R | TGCTGAACCGCTCTTCCGATCTNNNNNNNNAAAATCAAAAAATTCCCTTTC |

🟢 The earlier preprint's `media-1.xlsx` is a distinct source version. Its
`Sine.Alu_F` first-generation row prints
`AAAAATATAAAAAATTAGTTGGG` without the 5′ partial-adapter arm
`CTACACGACGCTCTTCCGATCT` printed in the published MOESM2 row above. The
published version is used by the model/page because the defining 2022 article
and its supplement are the final protocol source. The preprint's shorter row
cannot be silently substituted into the published library geometry; whether
the difference reflects a changed ordering design or table correction is not
established by the accessible text.

🟢 The second-generation indexed SINE Alu rows in that workbook give each
individual eight-base index and 0–5 N stagger; use the row corresponding to a
well rather than substituting the four pilot primers above. The article methods
(PMC8986802, “scTEM-seq library preparation”) say the 5′ partial-adapter arms
are completed with NEBNext dual-index oligos, catalogue E7780S. The matched
RNA branch follows Smart-seq2, so its RT and amplification oligos belong to that
source and the particular implementation.

🟡 The current shared model does not transcribe this published panel or a
particular NEBNext index choice. Its unresolved read-primer placement is thus
partly a model gap; the source does not identify one universal sample-specific
indexed primer. The large MOESM1 DOCX was over the source tool's size cap,
but the primer table itself was retrieved as MOESM2.
