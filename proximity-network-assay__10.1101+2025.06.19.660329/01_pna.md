# PNA / Pixelgen Proxiome — research notes

## Defining sources read

- 🟢 Karlsson F, Simonetti M, Galonska C, *et al.* **Cells as Single-Molecule Protein Proximity Networks.** bioRxiv (2025; version 3 posted 2026). DOI: [10.1101/2025.06.19.660329](https://doi.org/10.1101/2025.06.19.660329).
- 🟢 Pixelgen Technologies. **UM002: User Manual Pixelgen Proxiome Kit v2, Immuno 155 v2**, v1.00 (2026). The current download is linked from the [vendor manual page](https://www.pixelgen.com/user-manual-v2/).

## What the assay is

🟢 PNA is the method underlying the Pixelgen Proxiome Kit. It is related to PLA, but it does not use a separately designed reporter pair for every protein pair. Each targeted protein is represented by many copies of an antibody-derived PID–UMI address in a localized rolling-circle product; generic linker molecules then capture edges between nearby type-1 and type-2 nodes.

## Molecular sequence of operations

1. 🟢 PFA-fix cells, bind sample hashtags, pool eight samples, then bind the barcoded antibody panel and a second antibody reagent.
2. 🟢 Every antibody target has type-1 and type-2 oligo-conjugate forms. Each carries a protein ID (PID) and UMI.
3. 🟢 Hybridize the corresponding padlock probe. Polymerase copies PID and UMI across the gap and ligase closes the padlock.
4. 🟢 Short in-situ RCA produces a localized repeated PID–UMI node at each bound antibody.
5. 🟢 A proximity oligo hybridizes to one type-1 and one type-2 RCP. The second gap-fill/ligation copies both PID–UMI addresses into one linker amplicon. The type split prevents one RCP from supplying both sides of a link.
6. 🟢 Digest the cell/RCP scaffold, pre-amplify for four cycles, clean with SPRI beads, and add unique dual indexes in a second PCR.

## Library and reads

🟢 The published linker-amplicon order is:

`TruSeq R1 handle — UMI-1 — PID-1 — linker-binding sequence — linker — PID-2 — UMI-2 — Nextera R2 handle`

🟢 Proxiome v2 uses paired-end sequencing with minimum `Read 1 44`, `i7 8`, `i5 8`, `Read 2 78`, and 15% PhiX. Both indexes are acquired. The manual recommends demultiplexing on i7 alone and masking i5 by default because i5 quality can be lower.

🔴 The commercial manual names the padlock, proximity oligo, pre-amplification and indexing reagents but does not disclose their complete sequences. Exact primer landings and field boundaries therefore cannot be reconstructed from public product documentation.

## Boundary with Molecular Pixelation

🟢 PNA is not the same library chemistry as the earlier Molecular Pixelation (MPX) method (Karlsson *et al.*, *Nature Methods* 2024, DOI [10.1038/s41592-024-02268-9](https://doi.org/10.1038/s41592-024-02268-9)). MPX is useful lineage/context, while the PNA paper and Proxiome v2 manual define this page.
