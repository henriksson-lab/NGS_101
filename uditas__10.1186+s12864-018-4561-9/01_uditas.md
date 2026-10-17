# UDiTaS — research notes

Giannoukos et al., *BMC Genomics* (2018), [doi:10.1186/s12864-018-4561-9](https://doi.org/10.1186/s12864-018-4561-9), PMCID: PMC5861650.

## Assay target

The sequenced material is targeted human genomic loci, quantifying editing outcomes including structural rearrangements. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (UMI-bearing asymmetric tagmentation ends and target-anchored PCR primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Full open-access article and Methods.
- 🟢 Additional file 2, the detailed UDiTaS protocol and oligonucleotide table.
- Round 2026-10-17-breadth-first-21 bounded pass: PMC5861650 JATS lines 70–113 (Background/Results/Methods), Table 2 lines 155–231, Bioinformatics pipeline lines 235–259; Additional file 2 DOCX text sections 3–13 and Appendix A (lines ~20–325).

## Construct logic

1. 🟢 A custom hyperactive Tn5 transposome carries P5, one of eight i5 pooling barcodes, a 10-base random UMI, Nextera S5 and the mosaic end. Its bottom oligo is 5′ phosphorylated and ends in ddC.
2. 🟢 Tagmentation simultaneously fragments genomic DNA and marks random endpoints. Unlike a fixed pair of locus primers, the random endpoint lets molecules with short indels, long deletions, inversions and translocation junctions enter the same assay.
3. 🟢 Round-one PCR pairs the common P5 primer with one target-specific anchor carrying the full TruSeq Read 2 sequence at its 5′ end.
4. 🟢 Round-two PCR uses P5 again and a P7–i7 primer. The final molecule is dual-indexed: i7 is a conventional sample index, while the opposite index read traverses the 10-nt UMI and then the eight-base i5 pooling barcode.
5. 🟢 The defining experiments were sequenced on MiSeq (article Results). Additional file 2 section 13 run sheet: Read 1 300 bases, Index read 1 8 bases, Read 2 300 bases, Index read 2 18 bases; the article's pipeline (Methods, Bioinformatics pipeline iv) aligns *paired* reads. An earlier version of this note said there was no Read 2; that was a repository transcription error, corrected in round 2026-10-17-breadth-first-21.
6. 🟢 Appendix A tabulates the "i5 Barcode READ (Tx10 for UMI)" as the i5 barcode exactly as written in the top oligo followed by ten UMI cycles (e.g. N501 `TAGATCGCTTTTTTTTTT`). 🟡 That order is the forward-strand (MiSeq) Index 2 read primed from the grafted P5 oligo; the model declares that primer explicitly and derives barcode cycles 1–8 and UMI cycles 9–18 from the final construct. The earlier page described a reverse-complement Index 2 read (UMI first), which contradicted this table.

## Reaction details (Additional file 2)

- 🟢 Adapter annealing: top and bottom oligo, each 50 µM in 100 µL, 95 °C 2 min then slow ramp to 25 °C (~12 h). Transposome: 85.7 µL Tn5 (1.85 mg/mL, pTXB1-Tn5 expressed after Picelli et al. 2014) + 14.3 µL annealed oligos, 60 min room temperature.
- 🟢 Tagmentation: 50 ng gDNA, 2 µL transposome, 1× TAPS-DMF (TAPS-NaOH pH 8.5, MgCl2, DMF) in 20 µL, 55 °C 7 min; Zymo DNA Clean & Concentrator-5 cleanup.
- 🟢 Round 1 PCR: Platinum SuperFi master mix, 30 mM TMAC, 500 nM gene-specific primer, 250 nM i5 (P5) primer; 98 °C 2 min, 12 × (98 °C 10 s, 65 °C 10 s, 72 °C 90 s), 72 °C 5 min; AMPure XP 0.9×.
- 🟢 Round 2 PCR: SuperFi, 500 nM i5 primer, i7-barcoded primer; 15 cycles same profile; AMPure XP 0.9×.
- 🟢 Pool equimolar, AMPure 0.9×, BluePippin size selection 400–850 bp (1.5% cassette), AMPure 0.9×, KAPA qPCR quantification; MiSeq 10 pM with 15% PhiX.
- 🟡 The schematic does not yet draw the cleanups and BluePippin selection as graph separations (deferred).

## Represented branch

The page uses the published N501 Tn5 top oligo, N701 P7 primer and CEP290 guide-323 anchor as concrete examples. Article Table 2 (and Additional file 2 Appendix A) publishes eight Tn5 pooling-index top oligos (OLI6380–OLI6387; Table 2 leaves OLI6387 unnamed, Appendix A names it i5_N508), four locus anchors (OLI6062, OLI6256, OLI6253, OLI6259) and twelve i7 primers (N701–N707, N710–N712, N714, N715). These alternatives change index or target sequence, not the library architecture.
