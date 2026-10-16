# Slide-seqV2 — research notes

## 1. What it is

Slide-seqV2 captures tissue RNA on a monolayer of spatially decoded 10-µm beads. A bead-specific split-pool barcode and a molecule-specific random sequence travel with each captured transcript.

## 2. Sources read

- 🟢 Stickels et al., *Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2*, Nature Biotechnology 39, 313–319 (2021), DOI [10.1038/s41587-020-0739-1](https://doi.org/10.1038/s41587-020-0739-1); open full text [PMC8606189](https://pmc.ncbi.nlm.nih.gov/articles/PMC8606189/).
- 🟢 Online Methods lines for bead sequences, puck sequencing and library preparation.
- 🟢 [Supplementary Information, Table 6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-020-0739-1/MediaObjects/41587_2020_739_MOESM1_ESM.pdf), lists the template-switch, second-strand, cDNA PCR, TruSeq P5, and all puck-decoding ligation oligos with modifications. The Methods specify which ones enter each reaction.

## 3. Key oligo as printed

🟢 The Online Methods print three bead-capture oligos (lines 93–107):

| Bead preparation | 5′→3′ oligo as printed |
|---|---|
| Chemgenes | `TTTTTTTTCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNNT30` |
| In-house vs1 | `TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT25` |
| In-house vs2 | `TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNVVT30` |

The paper defines `PC` as a photocleavable linker, `J` as split-pool barcode bases, `N` as mixed UMI bases, and `V` as A/C/G. Fifteen J positions are synthesized; fourteen are used for spatial matching. The schematic selects the Chemgenes version, but displays only its barcode–UMI–oligo-dT core inside a normalized final Illumina duplex; the full bead oligo above is the source transcription.

🟢 Table 6 library-preparation oligos, 5′→3′ as printed (spaces caused by PDF line wrapping removed):

| Role | Ordered oligo |
|---|---|
| Template switch | `AAGCAGTGGTATCAACGCAGAGTGAATrG+GrG` |
| Second strand, dN-SMRT | `AAGCAGTGGTATCAACGCAGAGTGANNNGGNNNB` |
| TruSeq PCR handle | `CTACACGACGCTCTTCCGATCT` |
| SMART PCR primer | `AAGCAGTGGTATCAACGCAGAGT` |
| TruSeq P5 hybrid primer after Nextera | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` |

🟢 Puck decoding is itself an oligo-dependent sequencing reaction. Table 6 lists the following primer families; `N` is a degenerate base, `/5phos/` a 5′ phosphate, and `/3InvdT/` an inverted terminal dT:

| Role | Ordered oligo |
|---|---|
| TruSeq−1 / TruSeq / TruSeq+1 / TruSeq+2 | `/5phos/GATCGGAAGAGCGTCGTGTAG`; `/5phos/AGATCGGAAGAGCGTCGTGTAG`; `/5phos/NAGATCGGAAGAGCGTCGTGTAG`; `/5phos/NNAGATCGGAAGAGCGTCGTGTAG` |
| TruSeq(8b)+3 | `/5phos/NNNAGATCGGA/3InvdT/` |
| 3UP−1 / 3UP / 3UP+1 | `TCTCGGGAACGCTGAAG`; `TCTCGGGAACGCTGAAGA`; `TCTCGGGAACGCTGAAGAN` |
| UP−1 / UP / UP+1 / UP+2 | `/5phos/CTCGGGAACGCTGAAGA`; `/5phos/TCTCGGGAACGCTGAAGA`; `/5phos/NTCTCGGGAACGCTGAAGA`; `/5phos/NNTCTCGGGAACGCTGAAGA` |
| UP(7b)+3 / UP(7b)+4 | `/5phos/NNNTCTCGGG/3InvdT/`; `/5phos/NNNNTCTCGGG/3InvdT/` |
| Monobase 5A / 5G / 5C / 5T | `/6FAM/IIINNNAN`; `/aquaphluor593/IIINNNGN`; `/CY3/IIINNNCN`; `/Cy5/IIINNNTN` |
| Monobase 3A / 3G / 3C / 3T | `/5Phos/NANNNIII/6FAM/`; `/5Phos/NGNNNIII/aquaphluor593/`; `/5Phos/NCNNNIII/Cy3/`; `/5Phos/NTNNNIII/Cy5/` |

## 4. Construct-changing steps

1. Decode bead barcodes by in-situ ligation sequencing.
2. Hybridize tissue poly(A) RNA and reverse-transcribe with a template-switch oligo.
3. Digest tissue, recover beads, remove excess primers and synthesize the second strand.
4. PCR-amplify cDNA and use Nextera library preparation with the published TruSeq-P5 hybrid primer.

## 5. Evidence boundary

🟢 The three bead sequences, TSO, second-strand and preamplification primers, hybrid P5 primer, and puck-decoding oligos are printed in the paper and Table 6. 🟡 The schematic normalizes the final platform arms to canonical Illumina roles so sequencing-primer placement can be computed; it does not claim that normalization as another published oligo. The N700 kit index varies by sample and is not one fixed oligo. The publisher Table 6 is available, but the PMC-hosted Sup_info and Sup_Tab_3/4/5 attachments were gated at this review; source-review completeness remains open.
