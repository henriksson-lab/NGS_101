# TimeLapse-seq

**Defining paper:** Schofield et al., *Nature Methods* (2018), doi:10.1038/nmeth.4582. The open article and supplementary protocol were read from PMCID PMC5927882.

## Molecular boundary

🟢 Cells incorporate 4-thiouridine (4sU) into newly synthesized RNA during a pulse. Oxidation followed by 2,2,2-trifluoroethylamine treatment changes the pairing behavior of 4sU before reverse transcription.

🟢 Treated 4sU positions are recovered as T-to-C mismatches in ordinary RNA-seq reads. The mutation is the temporal label; it is not a UMI or sample barcode.

## Evidence boundary

🟡 TimeLapse chemistry is compatible with several downstream RNA-seq preparations. The schematic therefore fixes the recoding reaction and readout, but labels kit-dependent library ends by role.

## Oligo source audit

🟢 The paper's open author manuscript (`nmeth.4582_PMC5831505.xml.txt`, lines 74–88) distinguishes targeted TimeLapse sequencing from transcriptome-wide TimeLapse-seq. The targeted assay reverse-transcribes Actb/Gapdh or MYC with gene-specific primers, amplifies approximately 150-nt products, and adds TruSeq Index adapters. Transcriptome-wide libraries instead use the SMARTer Stranded Total RNA-Seq Pico Input kit. Neither kit's entire adapter/index set is a fixed TimeLapse chemistry constant.

🟢 Supplementary Table S1b (`nmeth.4582_41592_2018_BFnmeth4582_MOESM4_ESM.xlsx.txt`, lines 12–25) prints twelve targeted primer pairs: four mouse Actb, four mouse Gapdh and four human MYC. Every forward primer carries `CTACACGACGCTCTTCCGATCTNNNNN` before its gene-specific arm; every reverse primer carries `CAGACGTGTGCTCTTCCGATCT` before its arm. For example, the printed mActB1 pair is `CTACACGACGCTCTTCCGATCTNNNNN ACTGAGCTGCGTTTTACACCC` and `CAGACGTGTGCTCTTCCGATCT TCCTGAGTCAAAAGCGCCAAAAC` (spaces separate universal and gene-specific arms). The same table gives PCR amplification primers `CAGACGTGTGCTCTTCCGATC` and `CTACACGACGCTCTTCCGATCT`. This is a **targeted assay branch**; the page does not claim these Actb/Gapdh/MYC arms are universal TimeLapse-seq adapters.
