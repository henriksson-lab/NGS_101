# SPLASH — research notes

Defining source: Aw *et al.*, [10.1016/j.molcel.2016.04.028](https://doi.org/10.1016/j.molcel.2016.04.028). Detailed protocol: [10.3791/55255](https://doi.org/10.3791/55255).

## Evidence used

- 🟢 SPLASH expands to Sequencing of Psoralen crosslinked, Ligated, and Selected Hybrids.
- 🟢 A biotinylated psoralen crosslinks paired RNA; streptavidin enrichment selects those duplexes before proximity ligation.
- 🟢 The crosslink is reversed after ligation and the joined RNA partners are converted into a sequencing library.
- 🟢 The author protocol prints the activated 3′ adapter as `/5rApp/CTGTAGGCACCATCAAT/3ddC/`.
- 🟢 It prints the 96-nt RT primer as `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC/iSp18/CACTCA/iSp18/TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG`. The 3′ terminal 17 nt are complementary to the adapter.
- 🟢 After RT, the source gel-selects a product 96 bases longer than the ligated RNA, circularizes the single-stranded cDNA, and amplifies with reagents called “Universal PCR Primer” and “Index (X) Primer”. It specifies paired-end 150-base sequencing.
- 🔴 Neither the defining article nor the detailed protocol prints the Universal PCR Primer or Index (X) Primer sequences. Without those bases and the resulting circle-junction geometry, the final amplicon and exact read-primer binding sites cannot be reconstructed source-exactly.

## Protocol boundary

The biotin–streptavidin enrichment is the defining distinction from the 2D-gel selection in PARIS and RNase R enrichment in LIGR-seq.

The two disclosed custom oligos are rendered on the schematic page. The missing PCR-primer pair is kept explicit rather than guessed from partial Illumina motifs in the RT primer.
