# PNA / Pixelgen Proxiome — research notes

## Defining sources read

- 🟢 Karlsson F, Simonetti M, Galonska C, *et al.* **Cells as Single-Molecule Protein Proximity Networks.** bioRxiv (2025; version 3 posted 2026). DOI: [10.1101/2025.06.19.660329](https://doi.org/10.1101/2025.06.19.660329).
- 🟢 Pixelgen Technologies. **UM002: User Manual Pixelgen Proxiome Kit v2, Immuno 155 v2**, v1.00 (2026). The current download is linked from the [vendor manual page](https://www.pixelgen.com/user-manual-v2/).
- 🟢 Pixelgen Technologies. **UM001: User Manual Proxiome Kit Cell Surface**, v1.01, for the [historical Immuno 155 v1.1 kit](https://www.pixelgen.com/wp-content/uploads/2025/12/um001-user-manual-proxiome-kit-cell-surface-v.1.01.pdf). Table 12 gives R1 44, R2 78 and two 8-cycle index reads; Appendix Table A6 prints its sample-index windows. Its reagent list names padlock, proximity and PCR oligos without printing their full molecules. This revision is separate from the v2 assay drawn here.
- 🟢 Pixelgen Technologies. **Method for spatial mapping of cell targets.** [WO2025158200A1](https://patents.google.com/patent/WO2025158200A1/en). The patent documents the same two-conjugate, two-padlock, RCA and proximity-linking topology and identifies an incorporated sequence-listing file.
- 🟢 Pixelgen's open-source [Pixelator](https://github.com/PixelgenTechnologies/pixelator) software contains machine-readable `proxiome-v1` and `proxiome-v2` assay definitions, public PID panels, and real paired FASTQ test data. The values below use the Proxiome v2 definition at commit `20b5e4d2f16a7bc3afb7d133fc94b575d754b6b5`.

## Bounded pass 2026-10-17 (round breadth-first-15)

Read: preprint v3 PDF text lines 105–125 (Results, assay design), 178–187 (Figure 1 legend) and 598–620 (Methods, Proximity Network Assay and sequencing); UM002 v1.00 text Table of contents, Introduction lines 255–275, STEP 3–5 (barcoded/2nd antibody binding through indexing PCR clean-up) and STEP 6.1–6.3 (271-bp size, Table 16); UM001 v1.01 RCA step (line 1609) only for the RCA-time comparison. The scraper found no printed kit oligo sequences in these texts. Manual Figure 2a/2b (molecular description) and the patent sequence listing were not inspected in this pass.

## What the assay is

🟢 PNA is the method underlying the Pixelgen Proxiome Kit. It is related to PLA, but it does not use a separately designed reporter pair for every protein pair. Each targeted protein is represented by many copies of an antibody-derived PID–UMI address in a localized rolling-circle product; generic linker molecules then capture edges between nearby type-1 and type-2 nodes.

## Molecular sequence of operations

1. 🟢 PFA-fix cells, bind sample hashtags, pool eight samples, then bind the barcoded antibody panel overnight (16–20 h, 4 °C), wash four times, and stabilize it with the kit 2nd Antibody (10 min, 4 °C) followed by washes (UM002 STEP 3.1–3.2).
2. 🟢 Every antibody target has type-1 and type-2 oligo-conjugate forms. Each carries a protein ID (PID) and UMI.
3. 🟢 Hybridize the corresponding padlock probe (Padlock Oligos, 30 min at 45 °C, then washes). In Gap-fill 1 (30 min at 25 °C) polymerase copies PID and UMI across the gap and ligase closes the padlock (UM002 STEP 4.1–4.2).
4. 🟢 Short in-situ RCA produces a localized repeated PID–UMI node at each bound antibody: UM002 STEP 4.3 E (and UM001 v1.01) gives 15 min at 30 °C.
   ⚠️ Discrepancy, shown on the page: the preprint v3 Methods ("Proximity Network Assay") say RCA ran for 10 minutes while stating that the manufacturer's instructions were followed. The model uses the UM002 value for the drawn v2 kit; which manual revision the paper experiments followed is not stated.
5. 🟢 Proximity Oligos are added directly to the RCA reaction and hybridized for 30 min at 45 °C, then washed; each linker binds one type-1 and one type-2 RCP. Gap-fill 2 (30 min at 25 °C) copies both PID–UMI addresses into one linker amplicon (UM002 STEP 4.3 G–P, 4.4). The type split prevents one RCP from supplying both sides of a link.
6. 🟢 Count cells, digest 8,000 cells per pool with the kit Digestion program, pre-amplify for four cycles, clean with 1:1 SPRI/AMPure XP beads, add unique dual indexes in a second PCR (8, 6 or 5 cycles for 1,000, 4,000 or 8,000 cells), and clean again with 1:1 beads. The final library is approximately 271 bp (UM002 STEP 5.1–5.6, 6.1).
7. 🟡 Unclosed padlocks, unlinked RCPs and digested cell material are not described as physically removed; the page keeps them as unresolved fates rather than discards.

## Library and reads

🟢 The published linker-amplicon order is:

`TruSeq R1 handle — UMI-1 — PID-1 — linker-binding sequence — linker — PID-2 — UMI-2 — Nextera R2 handle`

🟢 Proxiome v2 uses paired-end sequencing with minimum `Read 1 44`, `i7 8`, `i5 8`, `Read 2 78`, and 15% PhiX. Both indexes are acquired. The manual recommends demultiplexing on i7 alone and masking i5 by default because i5 quality can be lower.

🟡 The reviewed preprint and commercial manual establish the roles, reaction order and
read lengths, and the manual prints index windows. We have not yet verified complete
kit oligo molecules from all linked sources. The patent refers to a separate ST.26
sequence listing, so their absence cannot be accepted without checking that listing
and mapping any relevant entries to the kit revision.

🟢 The v2 manual **does** print the seven supplied sample-index pairs in Table A6. In primer order 1–7, the i7 sequences are `TTCTCCAA`, `CACCTACT`, `ACAACCAC`, `GAGAACCT`, `ACTTACCG`, `CTAGCAGA`, `AGCCTGTC`; the forward-orientation i5 sequences paired with them are `TTGTGGAA`, `GGAGGTGT`, `AATTCGTG`, `GGTTGCTC`, `TAGGATGA`, `TAGACACG`, `ACACTGAT`. These are index windows, not complete indexing PCR primer molecules. Table A6 also prints their reverse complements for sample-sheet use.

🟢 The software fills in much of the missing internal construct. Its Proxiome v2 definition is `UMI-1(28) — PID-1(10) — CAAGTGACGCTGGTCAAACACTCATGTC — UEI(15) — CGCTTAGATGTCGG — PID-2(10) — UMI-2(28)`, totaling 133 nt. The public v2 panel supplies the 10-nt PID whitelist; for example, HLA-ABC is `ACTTCCTAGG` at both ends.

🟡 Combining that 133-nt definition with the documented 44/78 run gives a base-resolved read interpretation without guessing oligo synthesis sequences: R1 covers UMI-1 + PID-1 + 6 nt of LBS-1; R2 approaches from the other end and covers UMI-2 + PID-2 + all of LBS-2 + UEI + 11 nt of LBS-1.

🟢 The [manuscript dataset page](https://software.pixelgen.com/datasets/pna-methods-paper/) releases the original paired FASTQs and identifies the experiments as UM001 / Proxiome v1. Pixelator's v1 schema defines the paper-era 141-nt insert with 32-nt LBS-1 `CAAGTGACGCTGGGCATTCAAACACTCATGTC` and 18-nt LBS-2 `GCTTCGCTTAGATGTCGG`. Small ranged samples of the public FASTQs confirm the expected orientation and boundaries: R1 contains UMI-1(28) + PID-1(10) + 6 nt of LBS-1; R2 contains UMI-2(28) + PID-2(10) + reverse-complemented LBS-2(18) + UEI(15) + 7 nt of LBS-1. The current v2 design shortens each LBS by four bases.

🟡 The examined assay-definition file begins at the internal amplicon and does not
encode the full pre-amplification/indexing primers, padlocks, antibody-conjugate oligos
or physical flow-cell ends. The released FASTQs alone also cannot establish the
unsequenced ends; the patent sequence listing and kit sources are still open checks.

🟡 The related Pixelgen patent says that exemplary sequence elements and a separate ST.26 sequence-listing XML were filed. Those are evidence for a patent implementation, not automatically for the current Proxiome v2 kit; they should only be drawn after the individual patent sequences can be mapped unambiguously to COa, COb, PL1, PL2, linker and gap-fill oligos.

## Boundary with Molecular Pixelation

🟢 PNA is not the same library chemistry as the earlier Molecular Pixelation (MPX) method (Karlsson *et al.*, *Nature Methods* 2024, DOI [10.1038/s41592-024-02268-9](https://doi.org/10.1038/s41592-024-02268-9)). MPX is useful lineage/context, while the PNA paper and Proxiome v2 manual define this page.
