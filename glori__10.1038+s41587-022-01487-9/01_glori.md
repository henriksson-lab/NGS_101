# GLORI

**Defining source:** Liu *et al.*, *Nature Biotechnology* (2023),
doi:10.1038/s41587-022-01487-9. **Detailed protocol:** Shen *et al.*, *Nature
Protocols* (2024), doi:10.1038/s41596-023-00937-1. The methods, public data record
and protocol description were read.

## Evidence used

- 🟢 Glyoxal protects G and catalyses nitrite conversion of unmodified A to inosine;
  m6A resists conversion and remains A.
- 🟢 The detailed procedure specifies fragmentation, glyoxal protection, 750 mM nitrite
  conversion at 16 °C for 8 h plus 4 °C overnight, heat deprotection and end repair.
- 🟢 Defining libraries used modified eCLIP with a 10-nt inline UMI and paired-end HiSeq
  X Ten sequencing; NEBNext Small RNA and SMARTer alternatives were also demonstrated.
- 🟡 Exact complete eCLIP amplification oligos are not reproduced in the accessible text;
  their roles and structured UMI are shown without filling sequence gaps.

## Construct logic

GLORI is a negative conversion assay: an unmodified genomic A should be read as G, while
the retained-A fraction estimates m6A stoichiometry. The model therefore keeps chemical
state distinct from the downstream library branch and carries UMI identity explicitly.
