# INVADEseq — research notes

## Citation and sources read

- Galeano Niño JL, Wu H, LaCourse KD, et al. “INVADEseq to identify cell-adherent or invasive bacteria and the associated host transcriptome at single-cell-level resolution.” *Nature Protocols* 18, 3355–3389 (2023). [doi:10.1038/s41596-023-00888-7](https://doi.org/10.1038/s41596-023-00888-7). Full author manuscript and supplements read from PMC10790651.
- Galeano Niño JL, et al. “Effect of the intratumoral microbiota on spatial and cellular heterogeneity in cancer.” *Nature* 611, 810–817 (2022). [doi:10.1038/s41586-022-05435-0](https://doi.org/10.1038/s41586-022-05435-0). Original application linked by the protocol paper.
- The INVADEseq protocol cites 10x CG000086 Rev J for its v1 workflow. The later official CG000086 Rev K guide is locally readable and prints the 5′ gel-bead architecture, but its equivalence to the exact cited revision has not been established.

## What changes relative to ordinary 10x 5′ gene expression

🟢 The protocol adds a reverse primer against conserved bacterial 16S rRNA to the GEM reverse-transcription mix while retaining the standard poly(dT) RT primer. Both bacterial and host first strands template-switch onto the barcoded gel-bead oligo, so molecules from one GEM share a cell barcode but carry independent UMIs.

🟢 Amplified cDNA is split. One aliquot follows standard 10x 5′ GEX preparation. A second aliquot undergoes two targeted PCRs with the same enrichment-primer pair, with a 955–1,215 bp BluePippin selection between them, followed by Dual Index TT Set A PCR and another size selection.

🟢 This is targeted taxonomic capture, not bacterial whole-transcriptome sequencing: the 16S primer enriches a conserved ribosomal locus whose variable regions support taxonomic assignment.

## Oligos transcribed from the protocol

| Oligo | Sequence as published, 5′→3′ | Function |
|---|---|---|
| 16S 1100R | `GGGTTGCGCTCGTTG` | Spiked into GEM RT; reverse-primes bacterial 16S rRNA. |
| 16S Enrich Forward | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC` | P5-bearing primer that binds the 10x/barcode end and retains cell barcode plus UMI. |
| 1061R Enrich Nested Reverse | `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTCACGRCACGAGCTGACGAC` | Read-2-arm-bearing nested reverse primer within the original 1100R target. |

The IUPAC `R` in 1061R is retained verbatim. The forward primer is canonical P5 followed by a 20-nt 10x-end binding region. The reverse primer is the canonical TruSeq Read 2 sequence followed by the 20-nt nested 16S-binding region.

## Construct logic

🟢 The paper explicitly says that the forward enrichment primer aligns to the Illumina/10x end and preserves the cell barcode and UMI. It explicitly requires a nested reverse primer carrying the published Read 2 adapter sequence.

🟡 The completed dual-index construct is derived from those primers plus the specified 10x Dual Index TT PCR. The bacterial branch therefore retains the normal 10x orientation: Read 1 starts at the barcode end; Read 2 starts at the nested 16S end. Standard sequencing-primer identities are located from the final construct by the shared model.

🟡 The schematic uses a variable-length 16S cDNA region. The paper selects complete enrichment products between 955 and 1,215 bp, but that selection does not establish a single insert length or a fixed number of 16S bases in either read.

## Read layout

🟢 The bacterial enrichment library is sequenced on MiSeq v3 with paired 300-base reads. The GEX library follows the 10x 5′ workflow and is sequenced separately on an Illumina high-throughput instrument.

🟡 For the bacterial branch, raw Read 1 contains the 16-base cell barcode and 10-base UMI first, followed by the template-switch sequence and 16S-derived sequence. Read 2 enters the 16S amplicon from the nested 1061R end. Both biological reads are retained; the two index reads demultiplex the sample.

The protocol’s downstream text refers to a BAM-derived “Read 1” used for taxonomic analysis after Cell Ranger processing. That processed-file label must not be used to reverse the physical raw-read geometry of the completed library.

## Open boundary

🟢 The original 2022 application supplement (MOESM1, Methods lines 645–693) also prints 1100R and both full enrichment primers. It reports Next GEM 5′ v2 for the GEX library while citing CG000086 Rev J for the bacterial enrichment. The 2023 detailed protocol explicitly describes its main procedure as v1 and gives v2 alternatives. These are source-version boundaries, not grounds to assume that the two kits have identical internal oligos.

🟢 The later official CG000086 Rev K Sequences appendix (extracted PDF lines 6845–6904 and 7084–7233) prints a 5′ gel-bead TSO primer with a 16-base cell barcode, 10-base UMI and modified `rGrGrG` terminus, plus finished v1 library architecture. 🔴 The exact cited Rev J guide, the Next GEM v2 guide for the alternative branch, and the particular sample index well remain unverified. The page therefore uses only the source-established roles and explicitly leaves revision-specific geometry unresolved.

🟢 The 2023 protocol Materials list (JATS lines 154–156) identifies the v2 reagent kit as PN-1000263 and the separate 5′ library kit as PN-1000002, while its Step 85 note (line 758) calls PN-1000002 a "5′ v2 kit." This appears to conflate the two kit roles. The page follows the explicit Materials identities; which kit combination the Step 85 note intended remains uncertain.
