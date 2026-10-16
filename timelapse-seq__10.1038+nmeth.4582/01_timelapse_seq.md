# TimeLapse-seq

**Defining paper:** Schofield et al., *Nature Methods* (2018), doi:10.1038/nmeth.4582. The open author manuscript and supplementary protocol were read from PMCID PMC5831505.

## Molecular boundary

🟢 Cells incorporate 4-thiouridine (4sU) into newly synthesized RNA during a pulse. Oxidation followed by 2,2,2-trifluoroethylamine treatment changes the pairing behavior of 4sU before reverse transcription.

🟢 Treated 4sU positions are recovered as T-to-C mismatches in ordinary RNA-seq reads. The mutation is the temporal label; it is not a UMI or sample barcode.

## Evidence boundary

🟡 TimeLapse chemistry is compatible with several downstream RNA-seq preparations. The schematic therefore fixes the recoding reaction and readout, but labels kit-dependent library ends by role.

🟢 For the transcriptome-wide branch, the defining Online Methods (`nmeth.4582_PMC5831505.xml.txt`, RNA library preparation) use the SMARTer Stranded Total RNA-Seq Pico Input kit after chemical recoding. They do not establish a separate RNA-adapter ligation before kit entry, so the graph passes recoded RNA directly to the kit-based cDNA step. The targeted RT-PCR branch and its printed primers remain a separate implementation.

🟢 The **TT-TimeLapse-seq K562 enrichment branch** has an additional upstream capture step (defining JATS line 92). It cites [Duffy et al. 2015](https://doi.org/10.1016/j.molcel.2015.07.023): MTS-biotin forms a reducible disulfide with 4sU-containing RNA, followed by streptavidin capture. TimeLapse changes the washes, treats the beads with warm TE, cleans the recovered RNA, then chemically recodes and reduces it before SMARTer kit entry. This is a branch of the K562 experiment, not a universal step in the MEF or ordinary K562 libraries. The current graph depicts their shared recoding and kit boundary; it does not show the branch-specific capture chemistry.

🟢 **Published pulse discrepancy:** The original Supplementary Protocol (`nmeth.4582_41592_2018_BFnmeth4582_MOESM3_ESM.pdf.txt`, line 35, step 2) instructs a K562 pulse at **100 mM final 4sU** for 4 hours. The original Supplementary Note (`nmeth.4582_supp_NIHMS931206-supplement-1.pdf.txt`, lines 29–31) says the K562 four-hour experiment succeeded at **100 µM** and gives a successful range of **10 µM–1 mM**; the defining article's JATS reproduces the latter in its supplementary text. The model uses no numerical pulse concentration. For describing the experiment, this page uses the Supplementary Note's 100 µM observation because it fits the published working range; the intended value of the Protocol's 100 mM entry remains uncertain.

🟢 A [November 2016 original-generation SMARTer Pico Mammalian user manual](https://manuals.plus/m/fa0d970016b83cb070dba8c79d0f92b6efcf3bab7ff2ce4ba5361703a491b5dd.pdf) (cat. 635005–635007, original Takara guide at a third-party mirror) describes N6 random priming, LNA template switching, PCR1 addition of dual indexes and Illumina adapters, ZapR/R-Probe ribosomal cDNA depletion, and PCR2. Its Figure 2 and Sections V.A–E specify a ligation-free kit workflow. Its p. 7 prints the available i5/i7 barcode members; the paper does not identify which members were used. The guide labels the SMART Pico Oligos and TSO Mix sequences proprietary, so their exact bases and strand-resolved junctions remain unresolved. Whether the authors used this exact manual revision is also unverified.

## Oligo source audit

🟢 The defining Online Methods give different acquired read lengths for these implementations: targeted MEF Actb/Gapdh libraries used paired-end 75 bp (JATS line 76), transcriptome-wide MEF libraries paired-end 100 bp (line 88), and K562 libraries paired-end 150 bp (lines 92 and 96). The graph leaves one universal run length unspecified because it represents the recoding chemistry shared by these branches.

🟢 The paper's open author manuscript (`nmeth.4582_PMC5831505.xml.txt`, lines 74–88) distinguishes targeted TimeLapse sequencing from transcriptome-wide TimeLapse-seq. The targeted assay reverse-transcribes Actb/Gapdh or MYC with gene-specific primers, amplifies approximately 150-nt products, and adds TruSeq Index adapters. Transcriptome-wide libraries instead use the SMARTer Stranded Total RNA-Seq Pico Input kit. Neither kit's entire adapter/index set is a fixed TimeLapse chemistry constant.

🟢 Supplementary Table S1b (`nmeth.4582_41592_2018_BFnmeth4582_MOESM4_ESM.xlsx.txt`, lines 12–25) prints twelve targeted primer pairs: four mouse Actb, four mouse Gapdh and four human MYC. Every forward primer carries `CTACACGACGCTCTTCCGATCTNNNNN` before its gene-specific arm; every reverse primer carries `CAGACGTGTGCTCTTCCGATCT` before its arm. For example, the printed mActB1 pair is `CTACACGACGCTCTTCCGATCTNNNNN ACTGAGCTGCGTTTTACACCC` and `CAGACGTGTGCTCTTCCGATCT TCCTGAGTCAAAAGCGCCAAAAC` (spaces separate universal and gene-specific arms). The same table gives PCR amplification primers `CAGACGTGTGCTCTTCCGATC` and `CTACACGACGCTCTTCCGATCT`. This is a **targeted assay branch**; the page does not claim these Actb/Gapdh/MYC arms are universal TimeLapse-seq adapters.
