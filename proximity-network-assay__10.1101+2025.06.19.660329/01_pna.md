# PNA / Pixelgen Proxiome — research notes

## Defining sources read

- 🟢 Karlsson F, Simonetti M, Galonska C, *et al.* **Cells as Single-Molecule Protein Proximity Networks.** bioRxiv (2025; version 3 posted 2026). DOI: [10.1101/2025.06.19.660329](https://doi.org/10.1101/2025.06.19.660329).
- 🟢 Pixelgen Technologies. **UM002: User Manual Pixelgen Proxiome Kit v2, Immuno 155 v2**, v1.00 (2026). The current download is linked from the [vendor manual page](https://www.pixelgen.com/user-manual-v2/).
- 🟢 Pixelgen Technologies. **Method for spatial mapping of cell targets.** [WO2025158200A1](https://patents.google.com/patent/WO2025158200A1/en). The patent documents the same two-conjugate, two-padlock, RCA and proximity-linking topology and identifies an incorporated sequence-listing file.
- 🟢 Pixelgen's open-source [Pixelator](https://github.com/PixelgenTechnologies/pixelator) software contains machine-readable `proxiome-v1` and `proxiome-v2` assay definitions, public PID panels, and real paired FASTQ test data. The values below use the Proxiome v2 definition at commit `20b5e4d2f16a7bc3afb7d133fc94b575d754b6b5`.

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

🔴 A literal-sequence audit of the defining paper and commercial manual found no complete oligo molecules: both expose roles, ordering and read lengths, but not complete primer sequences.

🟢 The v2 manual **does** print the seven supplied sample-index pairs in Table A6. In primer order 1–7, the i7 sequences are `TTCTCCAA`, `CACCTACT`, `ACAACCAC`, `GAGAACCT`, `ACTTACCG`, `CTAGCAGA`, `AGCCTGTC`; the forward-orientation i5 sequences paired with them are `TTGTGGAA`, `GGAGGTGT`, `AATTCGTG`, `GGTTGCTC`, `TAGGATGA`, `TAGACACG`, `ACACTGAT`. These are index windows, not complete indexing PCR primer molecules. Table A6 also prints their reverse complements for sample-sheet use.

🟢 The software fills in much of the missing internal construct. Its Proxiome v2 definition is `UMI-1(28) — PID-1(10) — CAAGTGACGCTGGTCAAACACTCATGTC — UEI(15) — CGCTTAGATGTCGG — PID-2(10) — UMI-2(28)`, totaling 133 nt. The public v2 panel supplies the 10-nt PID whitelist; for example, HLA-ABC is `ACTTCCTAGG` at both ends.

🟡 Combining that 133-nt definition with the documented 44/78 run gives a base-resolved read interpretation without guessing oligo synthesis sequences: R1 covers UMI-1 + PID-1 + 6 nt of LBS-1; R2 approaches from the other end and covers UMI-2 + PID-2 + all of LBS-2 + UEI + 11 nt of LBS-1.

🟢 The [manuscript dataset page](https://software.pixelgen.com/datasets/pna-methods-paper/) releases the original paired FASTQs and identifies the experiments as UM001 / Proxiome v1. Pixelator's v1 schema defines the paper-era 141-nt insert with 32-nt LBS-1 `CAAGTGACGCTGGGCATTCAAACACTCATGTC` and 18-nt LBS-2 `GCTTCGCTTAGATGTCGG`. Small ranged samples of the public FASTQs confirm the expected orientation and boundaries: R1 contains UMI-1(28) + PID-1(10) + 6 nt of LBS-1; R2 contains UMI-2(28) + PID-2(10) + reverse-complemented LBS-2(18) + UEI(15) + 7 nt of LBS-1. The current v2 design shortens each LBS by four bases.

🔴 The software begins at the internal amplicon. It does not expose the full pre-amplification/indexing primers, padlocks, antibody-conjugate oligos or physical flow-cell ends, so those cannot be reconstructed from FASTQs alone.

🟡 The related Pixelgen patent says that exemplary sequence elements and a separate ST.26 sequence-listing XML were filed. Those are evidence for a patent implementation, not automatically for the current Proxiome v2 kit; they should only be drawn after the individual patent sequences can be mapped unambiguously to COa, COb, PL1, PL2, linker and gap-fill oligos.

## Boundary with Molecular Pixelation

🟢 PNA is not the same library chemistry as the earlier Molecular Pixelation (MPX) method (Karlsson *et al.*, *Nature Methods* 2024, DOI [10.1038/s41592-024-02268-9](https://doi.org/10.1038/s41592-024-02268-9)). MPX is useful lineage/context, while the PNA paper and Proxiome v2 manual define this page.
