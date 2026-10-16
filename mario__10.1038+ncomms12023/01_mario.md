# MARIO — research notes

Defining source: Nguyen *et al.*, [10.1038/ncomms12023](https://doi.org/10.1038/ncomms12023).

## Evidence used

- 🟢 MARIO maps RNA–RNA contacts mediated by a single protein in unperturbed cells.
- 🟢 UV crosslinking preserves RNPs; a biotinylated RNA linker is ligated so a selected chimera has the form RNA1–linker–RNA2.
- 🟢 Streptavidin purification selects linker-containing chimeras and a random barcode supports PCR-deduplication.
- 🟢 The defining paper's Methods section publishes every custom library oligo, rather than merely naming them: the internally biotinylated RNA linker, phosphorothioate-protected linker-complement oligo, activated RT adapter, three RT-primer variants, Cut_oligo, DP5/DP3 and full Illumina PE1/PE2 PCR primers.
- 🟢 The RNA linker is printed as `rCrUrArG/iBiodT/rArGrCrCrCrArUrGrCrArArUrGrCrGrArGrGrA`; the DNA oligo used to remove terminal linker is `T*C*G*C*ATTGCATGGGCTACTAGCAT`.
- 🟢 The pre-adenylated RT adapter is `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. The three 5′-phosphorylated RT primers contain random bases, the source's four-base sample barcode, a BamHI-bearing body, and the adapter-binding end.
- 🟢 The final amplification uses the fully printed Illumina PE PCR Forward Primer 1.0 and Reverse Primer 2.0, followed by paired-end HiSeq sequencing. Standard read-primer binding can therefore be derived from the completed Illumina shell.

## Protocol boundary

MARIO is protein-contact capture, not a psoralen duplex assay; its shared linker is therefore represented separately from PARIS, SPLASH and LIGR-seq.

The oligo table on the schematic preserves the source's internal modifications and all three experimental RT-primer variants. The diagram uses one generic UMI/sample-index-bearing payload, because those variants change identifiers rather than the reaction topology.

- 🟢 **Published discrepancy:** the original Methods section ‘Reverse transcription’ gives the generic RT-primer form `/5Phos/NNXXXXNNNN…` and says the random part has six bases; Supplementary Fig. 1 calls the whole inline barcode 10 nt. The original ‘Oligonucleotides used’ list prints the ES-1 primer `/5Phos/NNAGGTNNN…`, which has five random bases plus a four-base sample index (nine total), while the ES-2 and ES-indirect primers each have six random bases (ten total). The model copies the explicit ES-1 oligo; whether the ordered ES-1 primer actually had five random bases remains uncertain.

## Additional source conflict and linked-source audit

- 🟢 **Published enzyme conflict:** Methods ‘Proximity ligation’ (JATS line 149) says RNA 3′ phosphate is removed using PNK phosphatase activity, but its mixture names T4 PNK (3′-phosphatase minus). The preceding linker-ligation section (line 143) explicitly distinguishes this variant from phosphatase-active PNK. [NEB M0236 product information](https://www.neb.com/en-us/products/m0236-t4-polynucleotide-kinase-3-phosphatase-minus) says the named variant has no 3′ phosphatase activity. The model retains the intended ligation outcome; the exact dephosphorylation reagent remains unresolved.
- The article links its analysis software at [MARIO tools](https://mariotools.ucsd.edu/). The portal identifies v0.3 as the publication implementation; the linked legacy documentation could not be retrieved during this pass. Source coverage is therefore incomplete.
- The RT-primer/circularization procedure cites Huppertz et al., iCLIP, [10.1016/j.ymeth.2013.10.011](https://doi.org/10.1016/j.ymeth.2013.10.011), and terminal-linker removal cites T7-exonuclease and phosphorothioate-protection originals. Those linked method sources still need explicit independent inspection before accepting any missing-information conclusion.
