# scTrio-seq — research notes

## Sources read

- 🟢 Hou et al. 2016, DOI `10.1038/cr.2016.23`, PMC article and supplementary methods.
- 🟢 Tang et al. 2009, DOI `10.1038/nmeth.1315`, including the supplementary oligo table, for the inherited RNA-amplification chemistry.
- 🟢 Guo et al. 2013, DOI `10.1101/gr.161679.113`, for the inherited scRRBS branch.

## Protocol boundary

🟢 Cytoplasm and nucleus are separated. Cytoplasmic RNA supplies the transcriptome library. Nuclear DNA undergoes MspI reduced-representation digestion, methylated-adapter ligation and bisulfite conversion. The same DNA reads provide both CpG methylation and coarse copy number; there is no third molecular library.

## RNA branch

🟢 The scTrio methods say that the cytoplasmic supernatant is reverse-transcribed with poly(T) primers and amplified as in Tang et al. The inherited supplementary table prints UP1 (`ATATGGATCCGGCGCGCCGTCGAC` + T24), UP2 (`ATATCTCGAGGGCGCGCCGGATCC` + T24), and the corresponding 5′-amine-blocked AUP1/AUP2 reagents. The page now derives their binding sites from the shared Tang-lineage model.

🟢 scTrio then gel-selects 500–3,000-bp amplified cDNA, shears it to 150–250 bp, and prepares barcoded libraries with NEBNext Ultra DNA Prep for Illumina (E7370). 🔴 The scTrio paper does not print the proprietary kit adapter or index-primer bases. The page therefore does not reuse the obsolete SOLiD P1/P2 endpoint from Tang 2009; it marks only the later kit-generated Illumina arms as inferred.

## DNA branch

🟢 The nuclear fraction follows scRRBS: MspI digestion, recessed-end fill/dA-tailing, premethylated indexed TruSeq adapter ligation, bisulfite conversion, and PCR. Its adapter structure, final library, and Illumina sequencing-primer sites come from the shared scRRBS model rather than a second hand transcription.
