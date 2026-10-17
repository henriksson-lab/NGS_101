# TAPS

## 1. What it is

TAPS is a bisulfite-free methylation-sequencing method that converts 5mC and 5hmC, rather than unmodified cytosine, into bases read as thymine.

## 2. Authoritative source

🟢 Liu et al., *Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution*, Nature Biotechnology (2019), DOI `10.1038/s41587-019-0041-2`.

Sources read (round 2026-10-17-breadth-first-20):

- 🟢 Published Supplementary Information (Springer MOESM1): Supplementary Note 7 "TET-Assisted Pyridine borane Sequencing (TAPS)", the step-by-step 100-ng genomic protocol (text lines 924–1059); Supplementary Note 2 model-DNA oligos (lines 474–590).
- 🟢 Authors' accepted manuscript ("TAPS nbt final.pdf", Oxford ORA record uuid:91d3f09f-9c55-497b-9dc2-103420db9343): Methods "mTet1 Oxidation", "Pyridine borane reduction", "Whole-genome TAPS", "Low-input whole-genome TAPS", "Cell-free DNA TAPS" (text lines 1845–1920, 2240–2350). The typeset published article remains paywalled and was not read.
- The 2018 bioRxiv preprint is superseded and is not used for current claims.

## 3. Construct-changing operations

🟢 Supplementary Note 7: sonicate genomic DNA with spike-ins (e.g. 0.5% methylated lambda) to ~200 bp on a Covaris M220 and select 200–400 bp with 0.55X–1X AMPure XP; KAPA Hyper Prep Kit (KK8504) end repair, A-tailing and ligation (ligation extended to 30 min; post-ligation clean-up, elute 20 µl).

🟢 mTet1CD oxidation (37 °C, 80 min) with Fe(II), α-KG, ascorbate, ATP and DTT, then proteinase K (50 °C, 1 h), Bio-Spin P-30 and 1.8X AMPure XP. Supplementary Note 7 makes a second round optional "for more complete oxidation"; the accepted-manuscript whole-genome Methods oxidize twice, while low-input and cell-free DNA libraries oxidize once. The graph draws the two-round whole-genome variant.

🟢 Pyridine borane reduction: 1 M pyridine borane, 600 mM sodium acetate pH 4.3, 37 °C 16 h at 850 rpm; Zymo-Spin IC clean-up. TET oxidation produces 5caC, reduction produces dihydrouracil (DHU), and PCR reads DHU as thymine while unmodified cytosines stay C.

🟢 Library amplification: KAPA HiFi HotStart Uracil+ ReadyMix, 5 cycles (100 ng input; Methods use 8 and 13 cycles for 10 ng and 1 ng), using NEBNext Multiplex Oligos for Illumina index primers instead of the kit primer mix; 1X AMPure XP.

## 4. Oligos, library and read layout

🟢 Supplementary Note 7 ligation adapter ("Illumina multiplexing adapters without the phosphorothioate bonds", IDT, HPLC), annealed at 15 µM:

- `5’-{{= taps.ADAPTER_TOP }}`
- `5’-/5Phos/{{= taps.ADAPTER_BOTTOM }}`

🟡 Computed: the bottom oligo's first {{= len(taps.R2_STEM) }} bases pair with `{{= taps.R1_STEM }}` of the top oligo, leaving the top oligo's 3′-T to pair with the insert's 3′-dA and two unpaired fork arms (`{{= taps.R1_FORK }}` and `{{= taps.R2_FORK }}`). The ligated molecule therefore carries no P5/P7; these are added by PCR.

🟡 The universal and index primer sequences, and the index chosen for each sample, are not printed. The final library's P5 end, the i7-side extension beyond the adapter, the i7 index and P7 are drawn from the standard TruSeq layout and marked inferred.

🟢 Supplementary Note 2 prints a separate synthetic model library (TruSeq Universal Adapter, TruSeq Adapter Index 6 with `GCCAAT` as carried in the library, and a USER-removed uracil linker). It validates the chemistry and is not the genomic library construct; its Index 6 is not carried into the schematic.

🟢 Accepted-manuscript Methods: whole-genome TAPS libraries were sequenced paired-end 80 bp on a NextSeq 500 with 1% PhiX. The index-read cycle count is not stated.

## 5. Discrepancy and open questions

- Kit name: accepted-manuscript Methods say KAPA HyperPlus; published Supplementary Note 7 lists and uses the KAPA Hyper Prep Kit KK8504, after Covaris sonication. The model uses Supplementary Note 7 and shows the discrepancy on the page. Whether the typeset published Methods also say HyperPlus has not been checked.
- NEBNext Multiplex Oligos (E7335S) vendor original not read: primer sequences and index identity remain open.
- Acquired index-read cycles are not stated.
