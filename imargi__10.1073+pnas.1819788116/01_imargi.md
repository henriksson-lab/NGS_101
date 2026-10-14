# iMARGI — research notes

## Sources read

- 🟢 Yan et al. 2019, DOI `10.1073/pnas.1819788116`, PMC full text and supplement.
- 🟢 Wu et al. 2019, DOI `10.1038/s41596-019-0229-4`, PMC detailed protocol, especially the “iMARGI linkers and library PCR primers” section and Fig. 1.

## Construct path

🟢 iMARGI improves MARGI by performing the linker ligations in situ. Fragmented chromatin RNA is joined to one end of an annealed bivalent linker, and the linker’s DNA end is joined to nearby restriction-cut genomic DNA. Reverse transcription converts the RNA side, and paired-end sequencing reports RNA and DNA partners on opposite sides of the bridge.

## Exact published structural oligos

🟢 The detailed protocol prints both linker strands, the BamHI cut oligo and the PCR pair verbatim:

- Linker top: `5′-/Phos/NNAGATCGGAAGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCC*GA/iBiodT/CT-3′`. `NN` reduces RNA-ligation sequence bias; the source defines `C*` as a phosphorothioated C and `/iBiodT/` as internal biotin-dT.
- Linker bottom: `5′-/Phos/GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT-3′`.
- BamHI cut oligo: `5′-TCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCT/3InvdT/-3′`.
- NEBNext universal PCR: `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`.
- NEBNext index PCR: `5′-CAAGCAGAAGACGGCATACGAGAT[6 bp barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-S-T-3′`, retaining the protocol's printed terminal `-S-T` notation.

🟢 The top strand is 5′ adenylated after synthesis and annealed to the bottom strand. The bottom linker strand itself primes first-strand cDNA synthesis; the prior model's separate RT primer was unsupported and is removed. The cDNA is circularized, annealed to the cut oligo, and digested at the linker BamHI site. The split linker halves then serve as PCR landing sites. The run is paired-end 100 cycles. 🔴 The exact six-base sample-index member used in each library is not specified here; the source prints only a degenerate index slot. 🟡 The compact schematic still leaves the strand-resolved post-BamHI junction as a declared geometry limitation.
