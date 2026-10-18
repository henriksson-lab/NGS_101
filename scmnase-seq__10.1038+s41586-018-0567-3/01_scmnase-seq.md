# scMNase-seq research notes

The source review, evidence record, and inferred oligo reconstruction are preserved in
the combined [scDNase-seq / scMNase-seq note](../scdnase-seq__10.1038+nature15740/01_scdnase-scmnase.md#scmnase-seq-supplementary-methods-moesm1).

The defining source for this protocol is Lai *et al.*, *Nature* (2018),
doi:[10.1038/s41586-018-0567-3](https://doi.org/10.1038/s41586-018-0567-3).

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC8353605 JATS54–59; original MOESM1 PDF text27–48.

Added 23-cycle index PCR, 160–300 bp E-gel selection and source paired-end HiSeq2500 read roles. Removed upstream schematic as adapter evidence and assumed i7 acquisition; canonical geometry and six-base width remain illustrative. Public page records summarized recovery stages.

Deferred: Cited adapter method33 (Barski); exact adapter identity, index width and run lengths; graph recovery stages. This is a capped check, not a full source audit or confirmation of missing information.

## Source-owned graph migration, 2026-10-18

Replaced shared nuclease_accessibility recipe with a local scMNase graph. Added the source-described stop/proteinaseK and DNA recovery steps. No carrier assumption inherited from scDNase. Physical universal-adapter end geometry remains unavailable; inferred TruSeq end-prep/adapter/primer drawings are conditional references only. Removed the six-base illustrative index length from molecular declarations and use an unknown-length structured cell-index role. PCR source limitation references the actual ligated precursor, while the plan preserves source paired-end read roles and leaves run/index acquisition unknown.

Source scope: Lai etal.2018 published Supplementary Methods MOESM1 PDF pp2–3, extracted lines27–46: cell lysis/CaCl2, MNase incubation, SDS/EGTA stop/PK, phenol–chloroform/ethanol/glycogen recovery, cited universal-adapter ligation,23cycle indexPCR,160–300bp E-gel and paired-end HiSeq2500.

Deferred: Cited universal-adapter method33 Barski etal. remains deferred; exact study adapter/indexed primer identities and index length remain unassigned; Exact input-output molecular reaction relations and read/index cycle mask remain unresolved. This bounded check accepts no missing-information claim.
