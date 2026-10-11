# BotSeqS

**Source.** Hoang et al., *PNAS* (2016), [doi:10.1073/pnas.1607794113](https://doi.org/10.1073/pnas.1607794113), PMCID PMC5024639.

- 🟢 Libraries used the Illumina TruSeq DNA PCR-Free kit.
- 🟢 The bottleneck dilution is made after adapter ligation and before PCR, at about one haploid-genome equivalent.
- 🟢 For 500 ng input, the typical post-ligation dilution was 100,000-fold. The paper used 34 ng–1 µg input, 350-bp insert selection and 2 × 100-base HiSeq 2500 sequencing.
- 🟢 Read pairs sharing both genomic coordinates and orientation form a family; P5-derived Read 1 and P7-derived Read 2 distinguish Watson and Crick strands.
- 🟢 Native fragment endpoints are endogenous identifiers; BotSeqS adds no random UMI.
- 🟡 The paper says standard Illumina Y adapters but does not specify a dual-index arrangement. The schematic therefore shows the canonical single-index TruSeq molecule.

Dilution restricts the adapter-ligated originals entering PCR. Deep paired-end sequencing then supplies redundant Watson and Crick families for consensus.
