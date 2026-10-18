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

🟡 The completed dual-index construct is derived from those primers plus the specified 10x Dual Index TT PCR. The bacterial branch therefore retains the normal 10x orientation: Read 1 starts at the barcode end; Read 2 starts at the nested 16S end. The native complete sequencing-primer geometry remains unassigned; source roles and reported acquisition are kept separately.

🟡 The schematic uses a variable-length 16S cDNA region. The paper selects complete enrichment products between 955 and 1,215 bp, but that selection does not establish a single insert length or a fixed number of 16S bases in either read.

## Read layout

🟢 The bacterial enrichment library is sequenced on MiSeq v3 with paired 300-base reads (2023 protocol Step 121, JATS line 920; 2022 MOESM1 Methods lines 691–692). The GEX library follows the 10x 5′ workflow and is sequenced separately on an Illumina high-throughput instrument.

🟢 Host GEX acquired cycles are source-specific: the 2022 application supplement (MOESM1 lines 663–664) sequenced GEX on a NextSeq 2000 P3-100 flow cell at 26 × 10 × 10 × 90, which is a separate historical application run. The 2023 protocol's Equipment setup (JATS line 306) instead says HiSeq 4000 paired-end 75-bp reads, and Step 84 (line 732) leaves NextSeq 2000/NovaSeq read lengths to the user. 🟡 The current graph selects the2023 protocol paired75bp Equipment example and retains26/10/10/90 only as historical application context.

🟢 Library completion (2023 JATS lines 694–712 and 866–916; 2022 MOESM1 lines 683–690): GEX sample-index PCR uses 14–16 cycles depending on input; the bacterial branch takes 50 ng of second-enrichment product into a 15-cycle Dual Index TT Set A PCR, a single-sided SPRIselect cleanup, and then a repeat of the same 955–1,215-bp BluePippin program (Step 119). The model now shows that post-index selection as its own graph step. 🟡 The post-index length change against an unchanged collection window is not computed, because selected historicalTT primer/product assembly remains unassigned in this model.

🟡 For the bacterial branch, the laterRevK reference raw Read1 contains a16basebarcode and10baseUMI first, followed by the template-switch sequence and 16S-derived sequence. Read 2 enters the 16S amplicon from the nested 1061R end. That later-reference raw geometry is not proof of exacthistoricalRevJ bead/TT acquisition; nativewidths/indexcycles remain unassigned in the current graph.

The protocol’s downstream text refers to a BAM-derived “Read 1” used for taxonomic analysis after Cell Ranger processing. That processed-file label must not be used to reverse the physical raw-read geometry of the completed library.

## Open boundary

🟢 The original 2022 application supplement (MOESM1, Methods lines 645–693) also prints 1100R and both full enrichment primers. It reports Next GEM 5′ v2 for the GEX library while citing CG000086 Rev J for the bacterial enrichment. The 2023 detailed protocol explicitly describes its main procedure as v1 and gives v2 alternatives. These are source-version boundaries, not grounds to assume that the two kits have identical internal oligos.

🟢 The later official CG000086 Rev K Sequences appendix (extracted PDF lines 6845–6904 and 7084–7233) prints a 5′ gel-bead TSO primer with a 16-base cell barcode, 10-base UMI and modified `rGrGrG` terminus, plus finished v1 library architecture. 🔴 The exact cited Rev J guide, the Next GEM v2 guide for the alternative branch, and the particular sample index well remain unverified. The page therefore uses only the source-established roles and explicitly leaves revision-specific geometry unresolved.

🟢 The 2023 protocol Materials list (JATS lines 154–156) identifies the v2 reagent kit as PN-1000263 and the separate 5′ library kit as PN-1000002, while its Step 85 note (line 758) calls PN-1000002 a "5′ v2 kit." This appears to conflate the two kit roles. The page follows the explicit Materials identities; which kit combination the Step 85 note intended remains uncertain.

## Native graph migration — 2026-10-18 batch03

The source pass re-read one original: defining2023 protocol PMC10790651, sections listed in the matching bounded review. The graph selects that historical protocol, keeps the Materials/Step85 kit-role discrepancy explicit, and does not borrow internal bead/index geometry from laterRevK. Host and bacterial first RT products are separate single strands; template-switched copies are distinct from those first strands. Both physical aliquots retain the amplified mixed cDNA pool. Source introduction and Steps85–86 explicitly state that host cDNA remains the dominant background after bacterial enrichment. Bacterial target nonamplification is not a physical host-removal step. GEX uses sourcefragmentation/endprep, adapterligation and TT PCR inputs; bacterial35/20cycle enrichment has separate sourcecleanups, intermediateBluePippin and post-TT repeatedBluePippin. The2023 Equipment setup positively gives GEX paired75bp, whereas Step84 permits user-selected alternatives; the2022application26/10/10/90 run is a separate prior evidence scope, not this graph's run. BacterialPE300 is positive. TT acquisition/indexwidths and exactnativepriming remain unassigned. SourceindexPCR tables print72°C20min extension; this pass does not silently replace that value with a standard20s recipe.

Curator-written local review targets identify the source/implementation scope; targets are not derived from graph contents. This capped pass does not establish publication-wide absence or complete missing-information review.

### Implementation review corrections

Implementation correction draws source-supported opposed3′dA overhangs on the actual fragment duplex. Complementary role strands are derived with complement_segments to preserve matching widths. No phosphate or unknown adapter strand is supplied; ligation geometry stays unavailable. These corrections reference the original source pass above; no additional originals were read.
