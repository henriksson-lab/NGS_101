# Final chemistry-definition migrations

Status: complete as of 2026-10-19. The original 182-protocol cohort and the six
additional protocols listed below are migrated. PEM-seq had two queued dependencies;
an indirect-import review also found and migrated UMI-4C. No ownership queue remains open.

Source: [chemistry ownership queue](catalogue/chemistry_decoupling.tsv).

## Protocols

- [x] [10x Chromium Single Cell ATAC v1](10x-chromium-single-cell-atac/): `nextera`.
- [x] [Duplex-Repair](duplex-repair__10.1093+nar+gkab855/): `batch_ngs`.
- [x] [LamPORE](lampore__10.1101+2020.08.07.20161737/): `batch_ngs`.
- [x] [Olink Explore](olink-explore__10.1016+j.mcpro.2021.100168/): `batch_ngs`.
- [x] [PEM-seq](pem-seq__10.1038+s41421-019-0088-8/): `batch_ngs`, `protocol_model`.
- [x] [TIME-seq](time-seq__10.1038+s43587-023-00555-2/): `batch_ngs`.

## Completed migration scope

- `batch_ngs` (five protocols): replace shared library defaults with source-owned adapter, index and read declarations. Keep canonical sequences and neutral assembly/primer-site verification.
- `nextera` (Chromium ATAC): review `TAGMENTATION_OUTCOMES` / `amplifiable`; distinguish physical amplification, sequencer compatibility and computational filtering using protocol-specific evidence.
- `protocol_model.da_tailed` (PEM-seq): replace schematic end labels with explicit substrate strands and source-qualified end geometry; preserve unknowns.

For future additions, use protocol-owned native graph declarations with actual incoming molecular states, explicit identifier/modification lineage and physical fates. Complete source-qualified declaration review, independent critical implementation review and affected checks before resolving queue entries. Ownership completion does not certify all chemistry or resolve source-absence reviews.

## Status (2026-10-19)

All six are migrated, plus UMI-4C, which executed `protocol_model.da_tailed` through the
`conformation_replication_specs` re-export. Each has an independent implementation review.
See [round 2026-10-19-14](catalogue/migration_rounds/2026-10-19-14.tsv) and AUDIT.md
"Batch14 migration checkpoint". The `batch_ngs`, `nextera` and `protocol_model.da_tailed`
queue rows are `ownership_resolved` with no direct protocol callers.

## Code cleanup (2026-10-20)

Obsolete family modules, library constructors and graph/PCR shortcuts were retired
after checking explicit imports, Markdown expressions and remaining callers. Live
drawing primitives, page metadata and the Chromium partition helper remain. See the
cleanup checkpoint in [AUDIT.md](AUDIT.md) for scope and verification.
