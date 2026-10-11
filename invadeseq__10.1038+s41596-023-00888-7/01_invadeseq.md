# INVADEseq — research notes

## Citation and sources read

- Galeano Niño JL, Wu H, LaCourse KD, et al. “INVADEseq to identify cell-adherent or invasive bacteria and the associated host transcriptome at single-cell-level resolution.” *Nature Protocols* 18, 3355–3389 (2023). [doi:10.1038/s41596-023-00888-7](https://doi.org/10.1038/s41596-023-00888-7). Full author manuscript and supplements read from PMC10790651.
- Galeano Niño JL, et al. “Effect of the intratumoral microbiota on spatial and cellular heterogeneity in cancer.” *Nature* 611, 810–817 (2022). [doi:10.1038/s41586-022-05435-0](https://doi.org/10.1038/s41586-022-05435-0). Original application linked by the protocol paper.
- 10x Chromium Next GEM Single Cell 5′ v1/v2 kit documentation is the authority for the commercial backbone; the INVADEseq paper supplies the protocol-specific spike-in and enrichment primers.

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

🔴 The complete commercial gel-bead oligo and every version-specific reagent sequence are not printed in the protocol paper. The page therefore shows known barcode/UMI lengths and canonical Illumina sites, but does not invent proprietary internal bases.
