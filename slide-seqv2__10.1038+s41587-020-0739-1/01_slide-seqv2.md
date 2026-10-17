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

The paper defines `PC` as a photocleavable linker, `J` as split-pool barcode bases, `N` as mixed UMI bases, and `V` as A/C/G. Fifteen J positions are synthesized (Bead synthesis); fourteen are decoded on the puck. The printed vs2 oligo shows only 8 + 6 J (see section 4). The schematic selects the Chemgenes version, but displays only its barcode–UMI–oligo-dT core inside a normalized final Illumina duplex; the full bead oligo above is the source transcription.

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

🟢 Online Methods "Puck Sequencing" and "Slide-seqV2 library preparation" (PMC8606189 JATS text lines 119–158 and 195–333):

1. Decode bead barcodes: per primer, 5 µM in 4× SSC for 40 min, wash, 40-min T4 DNA ligase ligation with 5(base) or 3(base) fluorescent oligos, image, strip with 80% formamide; 8 five-prime, 3 three-prime and 3 SEDAL (2 h, primer in the ligation mix) ligations read 14 of 15 J bases.
2. Hybridize tissue RNA (6× SSC, 2 U/µl RNase inhibitor, 30 min room temperature).
3. First strand with template switching: Maxima H− RT and 50 µM TSO (10 µl in 200 µl), 1.5 h at 52 °C (v1 used 1 h at 42 °C with Ficoll).
4. Tissue digestion: equal volume of 2× SDS/EDTA/proteinase K buffer, 37 °C 30 min; beads pipetted off, glass discarded, three wash-buffer spins.
5. Exonuclease I (NEB M0568), 37 °C 50 min; three washes. 🟡 The Methods do not state the purpose.
6. 0.1 N NaOH, 5 min room temperature; quench with wash buffer; three washes.
7. Second strand: dN-SMRT oligo (2 µl of 1 mM) with Klenow fragment (NEB M0210) in Maxima RT buffer with dNTPs, 37 °C 1 h; three washes.
8. cDNA PCR on beads: Terra PCR Direct (Takara 639270) with TruSeq PCR handle and SMART PCR primers; 95 °C 3 min, 4 cycles (98/65/72 °C), 9 cycles (98/67/72 °C), 72 °C 5 min.
9. Two successive 0.6× AMPure XP cleanups (30 µl beads per 50 µl); 600 pg tagmented with Nextera XT.
10. Library PCR with Truseq5 (Table 6 "Truseq P5") and N700-series index primers, 12 cycles; 0.6× AMPure XP.
11. NovaSeq S2 100-cycle kit, 42 bases Read 1, 8 bases i7, 50 bases Read 2.

🟡 Final library (derived): P5–TruSeq Read 1 arm (the bead's TruSeq PCR handle extended by Truseq5)–J8–`TCTTCAGCGTTCCCGAGA`–J7–UMI8–T30–cDNA–Nextera mosaic end–S7–i7–P7. TruSeq Read 1 covers the 41 nt of barcode, constant and UMI in its 42 cycles; Nextera Read 2 and i7 read the tagmented end.

🔴 Custom bead vs2 prints 8 + 6 J bases (line 105), whereas Bead synthesis states 15 split-pool cycles in blocks of 8 and 7 (line 111). The page shows this notice and models only the Chemgenes oligo.

## 5. Evidence boundary

🟢 The three bead sequences, TSO, second-strand and preamplification primers, hybrid P5 primer, and puck-decoding oligos are printed in the paper and Table 6. 🟡 Model correction (2026-10-17): the schematic previously drew a TruSeq adapter on both ends, put all 15 J bases before the bead constant and used Nextera i5/i7 indexing primers. It now draws the Truseq5 × Nextera N7xx endpoint printed in the Methods, the J8–constant–J7 bead layout, and the Exonuclease I, NaOH and Klenow second-strand steps. The Illumina P5/P7, S7 and mosaic-end bases come from the canonical sequence modules, not from this paper. The N700 kit index varies by sample and is not one fixed oligo. The publisher Table 6 is available, but the PMC-hosted Sup_info and Sup_Tab_3/4/5 attachments were gated at this review; source-review completeness remains open.
