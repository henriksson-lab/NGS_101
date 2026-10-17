# scTrio-seq — research notes

## Sources read

- 🟢 Hou et al. 2016, DOI `10.1038/cr.2016.23`, PMC article and supplementary methods.
- 🟢 Tang et al. 2009, DOI `10.1038/nmeth.1315`, supplementary oligo table `nmeth.1315_41592_2009_BFnmeth1315_MOESM306_ESM.xls.txt` lines 3–6, for the inherited RNA-amplification chemistry. Declared in `catalogue/extra_sources.tsv` so the fetch is reproducible.
- 🟢 Guo et al. 2013, DOI `10.1101/gr.161679.113`, for the inherited scRRBS branch. Also declared in `catalogue/extra_sources.tsv`.
- 🟢 Tang et al. 2010, DOI `10.1038/nprot.2009.236`, same-author detailed protocol, PMC3847604, for Exonuclease I, TdT, UP2 second-strand synthesis and amine-blocked re-amplification (JATS lines 281–316, 338–436, 487–531). The original 2009 oligo table supplies the exact primer bases.
- 🟢 NEBNext Ultra DNA Library Prep E7370 [official manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7370.pdf), version 9.1 (2023), ligation and PCR sections: adapter and index oligos are supplied in separate NEBNext oligo kits. This later revision cannot establish the exact 2016 kit configuration; the scTrio article does not identify its chosen oligo kit.

## Protocol boundary

🟢 Cytoplasm and nucleus are separated. Cytoplasmic RNA supplies the transcriptome library. Nuclear DNA undergoes MspI reduced-representation digestion, methylated-adapter ligation and bisulfite conversion. The same DNA reads provide both CpG methylation and coarse copy number; there is no third molecular library.

## RNA branch

🟢 The scTrio methods say that the cytoplasmic supernatant is reverse-transcribed with poly(T) primers and amplified as in Tang et al. The inherited supplementary table prints, in spaced triplets, UP1 (`5’-ATA TGG ATC CGG CGC GCC GTC GAC TTT TTT TTT TTT TTT TTT TTT TTT-3’`, i.e. `ATATGGATCCGGCGCGCCGTCGAC` + T24), UP2 (`5’-ATA TCT CGA GGG CGC GCC GGA TCC TTT TTT TTT TTT TTT TTT TTT TTT-3’`, i.e. `ATATCTCGAGGGCGCGCCGGATCC` + T24), and `NH2-UP1`/`NH2-UP2` — the *same* 48 bases carrying a 5′-NH<sub>2</sub>, which is why the model's AUP1/AUP2 are the UP1/UP2 bases rather than distinct oligos. Both poly(T) tails are exactly 24 T (eight printed triplets), counted from the extracted table rather than eyeballed. The page now derives their binding sites from the shared Tang-lineage model.

🟢 scTrio then gel-selects 500–3,000-bp amplified cDNA, shears it to 150–250 bp, and prepares barcoded libraries with NEBNext Ultra DNA Prep for Illumina (E7370). 🔴 The scTrio paper does not identify a matching NEBNext oligo kit or selected index, and the E7370 manual offers multiple separately purchased oligo-kit options. The page therefore does not reuse the obsolete SOLiD P1/P2 endpoint from Tang 2009; it marks the kit-generated Illumina arms as inferred.

## DNA branch

🟢 The nuclear fraction follows scRRBS: MspI digestion, recessed-end fill/dA-tailing, premethylated indexed Illumina adapter ligation, bisulfite conversion, and PCR. 🟡 The shared scRRBS model illustrates a canonical TruSeq adapter; neither the scRRBS source nor the scTrio paper prints the exact ligation-facing oligos or proves that adapter geometry. Its final-library arms and sequencing-primer positions are conditional on that inferred design.

🟢 scTrio Methods ("Single-cell segregation…", PMC JATS line 154): the bottom 3 µl nuclear fraction receives protease and 60 fg unmethylated lambda DNA before scRRBS digestion and bisulfite conversion; the 4 µl supernatant goes into a tube pre-loaded with ERCC spike-in RNA and reverse transcriptase. 🟢 Guo et al. 2013 Methods ("Construction of single-cell RRBS sequencing libraries", JATS line 105): Klenow exo– fill-in/dA-tailing, "Illumina standard premethylated indexed adaptors" with T4 DNA ligase, MethylCode bisulfite conversion with on-column desulfonation and tRNA carrier, then two PCR rounds (PfuTurbo Cx, 25 cycles; Phusion, 22 cycles) and 200–500-bp native PAGE size selection. The graph now shows the protease/lambda release step and the two PCR rounds; the round-1 primer sequences were not in the inspected passage and are inferred inputs.

## Cleanup and read layout

🟢 scTrio Methods ("Sample quality control and library construction", JATS line 158): amplified cDNA is purified with DNA Clean & Concentrator 5 and gel-selected at 500–3,000 bp, which the paper says removes the amplification primers; the graph shows this selection and its discarded gel waste before Covaris S2 shearing. 🟢 All constructed libraries were sequenced 100-bp paired-end on HiSeq 2000 or HiSeq 2500; the model declares Read 1 and Read 2 as 100 cycles. 🔴 Index-read cycle counts are not stated in the inspected passages.

## Bounded pass 2026-10-17 (round 17)

Inspected: Hou et al. 2016 PMC4783472 JATS Methods lines 152–168 and Supplementary Tables S1/S2 captions (lines 525–531); Guo et al. 2013 PMC3847781 JATS Methods lines 103–107; scraper run over the bundle (only hits outside Tang/Guo tables are the four bisulfite-qPCR QC primers, lines 160–166, which do not enter the library). Deferred: Guo 2013 supplementary zip (round-1/round-2 primer identities, QP1/QP2), Tang 2009/2010 re-read, NEBNext E7370 revision contemporary with 2016, and index-read cycles.
