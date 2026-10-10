# miCLIP for m6A and m6Am

## Sources read

- 🟢 Linder et al., *Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome*, Nature Methods (2015), [doi:10.1038/nmeth.3453](https://doi.org/10.1038/nmeth.3453).
- 🟢 Grozhik et al., *Mapping m6A at individual-nucleotide resolution using crosslinking and immunoprecipitation*, Methods in Molecular Biology (2017), PMCID `PMC5562447`, for the complete experimental path and oligos.

## Molecular path

🟢 Fragmented poly(A) RNA is bound by anti-m6A antibody and irradiated at 254 nm. Protein A/G purification isolates covalent antibody–RNA complexes. A pre-adenylated 3′ RNA adapter is ligated before proteinase K leaves a peptide remnant at the crosslink.

🟢 Reverse transcriptase either reads through the adduct with a characteristic substitution or stops at it. These are alternative molecular outcomes: a mutation marks the crosslinked position in a full cDNA, whereas the end of a truncated cDNA marks an RT arrest.

🟢 First-strand cDNA is circularized, restriction-linearized and PCR-amplified. The sequenced molecule begins `NNN-BBBB-NN`: five random bases surrounding a four-base experimental barcode. The random bases are used together to identify PCR duplicates.

## Evidence boundary

🟡 The page models the inline identifier layout and small-RNA Read 1 landing site. It does not present a reconstructed order-form sequence where the detailed protocol's complete oligo table has not been independently transcribed.
