# Slide-seq (v1)

**Defining paper:** Rodriques et al., *Science* (2019), doi:10.1126/science.aaw1219. The open author manuscript was read from PMCID PMC6927209, together with its published Supplementary Materials (PMC NIHMS1036108 supplement-1.docx, hand-placed by the maintainer; title and author list match the article). Earlier passes used the bioRxiv 563395 v1 supplement; that preprint is now superseded for every value on the page.

## Molecular boundary

🟢 Uniquely DNA-barcoded 10-µm beads are packed into a puck and decoded spatially by SOLiD sequencing-by-ligation. A fresh-frozen tissue section is transferred to that mapped surface.

🟢 Released mRNA is captured on bead oligos and converted into a Drop-seq-style 3′ RNA library. The bead barcode supplies position and the UMI identifies the captured molecule.

🟢 The published supplement (docx text lines 22–31) prints two bead-oligo designs. Sequence 2 contains the TruSeq handle, an 8-base split-pool barcode block, `TCTTCAGCGTTCCCGAGA`, a 7-base barcode block, an 8-base random UMI and oligo(dT)30. Sequence 1 has a 6-base first barcode block and additionally carries a photocleavable linker (`TTTT-PC-GCCGG…`) and a T7-promoter region. The preprint wrote this as a photocleavable thymidine (`PCTGCCGG…`); the published wording is used.

🟢 Published Table S1 (lines 310–326) prints the 58-nt TruSeq5 tagmentation primer, the 23-nt SMART PCR primer and the 22-nt TruSeq PCR-handle primer. Libraries were sequenced as 42 cycles Read 1, 8 cycles i7 and 50 cycles Read 2 on a NovaSeq S2 flow cell. The page therefore models the source-defined mixed endpoint: TruSeq Read 1 reads bead barcode and UMI, while Nextera Read 2 enters the transcript insert from the tagmented end.

## Library steps (published supplement, docx text lines 39–129)

Provenance: the published author manuscript (PMC6927209, line 137) says materials and methods are in the online supplement. Round deferred-followup-03 read the PMC author-manuscript supplement-1.docx; the Science-hosted `aaw1219_Rodriques_SM.pdf` (HTTP 403) was not obtained. A word-level comparison of the library-preparation section with the bioRxiv v1 supplement found only unit formatting and one reworded Bioanalyzer sentence, so every condition below is unchanged; the line numbers now point into the published docx text.

1. 🟢 Puck decoding (lines 39–51): per primer, hybridize 5 µM primer in 4× SSC for 20 min, wash in SOLiD buffer F, ligate 20 min with T4 DNA ligase (rapid) and 40×-diluted SOLiD SR-75 sequencing oligo, cleave with SOLiD buffers C then B (twice), ligate again, then strip with 80% formamide. Two ligations on each of ten primers (Table S1), six of them constant bases, give 14-base bead barcodes (line 44).
2. 🟢 Hybridization (lines 64–65): puck in 6× SSC with 2 U/µl RNase inhibitor for 15 min at room temperature.
3. 🟢 First strand (lines 67–77): 1 h at 42 °C with Maxima H− RT, Ficoll PM-400, dNTPs, RNase inhibitor and 50 µM template-switch oligo (10 µl in 200 µl).
4. 🟢 Tissue digestion (lines 79–87): equal volume of 2× buffer (SDS, EDTA, proteinase K), 37 °C for 40 min.
5. 🟢 Bead release and washes (lines 89–95): beads pipetted off, glass discarded, three Tris/EDTA/Tween-20 washes with 3000 RCF spins, then water.
6. 🟢 cDNA PCR (lines 97–115): KAPA HiFi HotStart ReadyMix with Truseq_PCR_handle and SMART PCR primers; 95 °C 3 min, 4 cycles (98/65/72 °C), 9 cycles (98/67/72 °C), 72 °C 5 min.
7. 🟢 Cleanup and tagmentation (lines 117–129): 0.6× AMPure XP; 600 pg tagmented with Nextera XT; 12-cycle PCR with Truseq5 and N700-series index primers; 0.6× AMPure XP.
8. 🟢 Sequencing (line 129): NovaSeq S2, 42 bases Read 1, 8 bases i7, 50 bases Read 2.

## Published discrepancies

🟢 Published Table S1 (docx text lines 314 and 316; the bioRxiv v1 Table S1, rendered p.22, prints the same values) prints the TSO as `AAGCTGGTATCAACGCAGAGTGAATrG+GrG` and the SMART PCR primer as `AAGCAGTGGTATCAACGCAGAGT`. The printed TSO lacks the primer's `AG` at positions 5–6 and contains only its 3′ 17 nt. 🟢 The same Qiagen order number (#339414YCO0076714, supplement line 111) is printed in the Slide-seqV2 Methods, whose Supplementary Table 6 gives `AAGCAGTGGTATCAACGCAGAGTGAATrG+GrG`. 🔴 Whether this is a transcription error is unresolved; the published supplement repeats the preprint value. The page shows this notice, displays the TSO as printed and keeps the cDNA PCR geometry unavailable.

🟢 Decoding primers (published Table S1 lines 317–326, all `/5Phos/`): Truseq `AGATCGGAAGAGCGTCGTGTAG`, Truseq -1 `GATCGGAAGAGCGTCGTCTAG`, Truseq -2 `ATCGGAAGAGCGTCGTGTAG`, TruSeq-3 `TCGGAAGAGCGTCGTGTAG`, TruSeq-4 `CGGAAGAGCGTCGTGTAG`, UP `TCTCGGGAACGCTGAAGA`, UP-1 `CTCGGGAACGCTGAAGA`, UP-2 `TCGGGAACGCTGAAGA`, UP-3 `CGGGAACGCTGAAGA`, UP-4 `GGGAACGCTGAAGA`. 🟡 Computed in `chem.DECODING_MISMATCHES`: Truseq is the reverse complement of the bead TruSeq handle and UP of the bead linker, and each -n primer is its parent recessed by n 5′ bases. Nine primers match. Truseq -1 is printed with C where the series gives G (`…GTCGTCTAG` against `…GTCGTGTAG`, fourth base from the 3′ end). 🔴 Whether the printed base is a typo is unresolved. The page shows a discrepancy notice and lists the primer as printed. The decoding step is modelled at role level, so neither value enters the graph.

## Open questions

- Science-hosted supplement PDF (`aaw1219_Rodriques_SM.pdf`): not obtained; identity with the PMC author-manuscript supplement unchecked.
- Drop-seq bead synthesis (Macosko 2015) and SOLiD SR-75 vendor chemistry: not followed in this pass.
- Decoding primers are now printed (Table S1). Still open: which primer and ligation pairs read which barcode positions, and the SOLiD probe chemistry (vendor). The decoding step stays at role level.

## Evidence boundary

🟢 This page describes the original Slide-seq workflow, not Slide-seqV2, whose template-switching and second-strand changes have their own page.


## Critical re-review, 2026-10-17

Candidate internal/SMART-end tagmentation fragments remain unresolved. The represented PCR and AMPure steps do not demonstrate class-specific physical removal.

Inspected scope: Local tagmentation, PCR and AMPure declarations and terminal-fate audit; no new reaction-source assertion. This bounded correction does not complete the linked-source or two-reader absence audit.
