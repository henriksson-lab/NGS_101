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

The published article and its supplements define the protocol version shown here. Superseded preprint designs are not used.

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

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC8986802 JATS124–135, Library preparation and Sequencing.

Separated35-cycle TE PCR from1.2xAMPure purification/pooling and5-cycle NEBNext indexing PCR. Named uracil-tolerant KAPA polymerase. Removed superseded-preprint primer comparison from public caveat; published primer row remains authoritative.

Deferred: Published SupplementaryTablesS1/S2 full primer-family validation and binding geometry, final0.9xAMPure cleanup, transcriptome predecessor and exact run lengths. This capped pass does not establish a complete source review or publication-wide absence.
