# miCLIP for m6A and m6Am

## Sources read

- 🟢 Linder et al., *Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome*, Nature Methods (2015), [doi:10.1038/nmeth.3453](https://doi.org/10.1038/nmeth.3453).
- 🟢 Grozhik et al., *Mapping m6A at individual-nucleotide resolution using crosslinking and immunoprecipitation*, Methods in Molecular Biology (2017), [doi:10.1007/978-1-4939-6807-7_5](https://doi.org/10.1007/978-1-4939-6807-7_5), PMCID `PMC5562447`, section 2.3 “Oligonucleotides” read directly from the PMC XML.

## Oligos from the detailed 2017 protocol

🟢 L3 linker (HPLC purified): `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. Cut oligo (PAGE purified): `GTTCAGGATCCACGACGCTCTTC/3ddC/`.

🟢 P5 Solexa PCR primer: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`. P3 Solexa PCR primer: `CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT`.

🟢 Sixteen RT primers share `/5phos/NN[4-base experimental barcode]NNNAGATCGGAAGAGCGTCGTGgatcCTGAACCGC`. The protocol prints barcodes `aacc`, `acaa`, `attg`, `aggt`, `cgcc`, `ccgg`, `ctaa`, `catt`, `gcca`, `gacc`, `ggtt`, `gtgg`, `tccg`, `tgcc`, `tatt`, `ttaa` in RT1–RT16 respectively. RT1 writes its first random position as lowercase `n` (`NnaaccNNN`), while the other rows use uppercase `N`; retain that spelling when transcribing an individual primer. The 5′ phosphate and lowercase `gatc` are source modifications/case notation, not guessed bases.

## Molecular path

🟢 Fragmented poly(A) RNA is bound by anti-m6A antibody and irradiated at 254 nm. Protein A/G purification isolates covalent antibody–RNA complexes. A pre-adenylated 3′ RNA adapter is ligated before proteinase K leaves a peptide remnant at the crosslink.

🟢 Reverse transcriptase either reads through the adduct with a characteristic substitution or stops at it. These are alternative molecular outcomes: a mutation marks the crosslinked position in a full cDNA, whereas the end of a truncated cDNA marks an RT arrest.

🟢 First-strand cDNA is circularized, restriction-linearized and PCR-amplified. The sequenced molecule begins `NNN-BBBB-NN`: five random bases surrounding a four-base experimental barcode. The random bases are used together to identify PCR duplicates.

## Evidence boundary

🟡 The page models the inline identifier layout and small-RNA Read 1 landing site. The complete oligo ordering strings above are verified, but their alignment into the final duplex has not yet been reconciled base by base; the schematic remains a role-level representation.
