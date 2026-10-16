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
