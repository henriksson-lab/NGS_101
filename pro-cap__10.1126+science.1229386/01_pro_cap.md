# PRO-cap

**Defining paper:** Kwak et al., *Science* (2013), doi:10.1126/science.1229386. The open author manuscript was read from PMCID PMC3974810.

## Molecular boundary

🟢 Precision run-on permits engaged polymerases to incorporate one or very few biotin-NTPs. Streptavidin selection enriches those nascent molecules.

🟢 PRO-cap modifies the same selected pool to recover capped RNA 5′ ends. Enzymatic cap selection therefore changes the informative end from the PRO-seq active-site 3′ end to the initiating 5′ nucleotide.

## Evidence boundary

🟡 The original paper uses standard Illumina RNA adapters and RT/PCR primers for PRO-cap, without printing their bases in the author-hosted supplement. A later author protocol prints a complete PRO-cap oligo set, but its release is not silently substituted into the defining 2013 construct.

## Oligo source audit

🟢 Checked the original Science supplement hosted by the authors (`594009_926d2d49ebf94caa8e6cf707c2f4846a.pdf.txt`, lines 118–147). It prints modified **reverse PRO-seq** adapters, but explicitly says PRO-cap instead uses standard Illumina 3′ and 5′ RNA adapters, RT primer, and PCR/index primers. Those standards are named rather than sequenced in the 2013 supplement; copying the adjacent reverse PRO-seq oligos into PRO-cap would be wrong.

🟢 The authors' later detailed protocol (`PMC5502525_PMC5502525.xml.txt`, lines 2885–2902; [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5502525/)) specifies a concrete PRO-cap implementation: RA3 `UGGAAUUCUCGGGUGCCAAGG-/Inverted dT/`, RA5 `GUUCAGAGUUCUACAGUCCGACGAUC`, RTP `GCCTTGGCACCCGAGAATTCCA`, RP1 `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`, and an RPI-n index-primer family. These are useful source-verified oligos for that 2016 implementation, but the 2013 paper does not establish that every base and modification was identical. The schematic retains role-labelled ends pending a version-specific reconstruction.
