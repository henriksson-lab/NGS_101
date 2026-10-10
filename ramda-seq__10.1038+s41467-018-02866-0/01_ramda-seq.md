# RamDA-seq — source notes

Defining paper: Hayashi et al., *Nature Communications* 2018, [doi:10.1038/s41467-018-02866-0](https://doi.org/10.1038/s41467-018-02866-0).

## Sources read

- 🟢 PMC5809388 full-text XML, Methods “RT-RamDA”, “Not-so-random primers”, “RamDA-seq”, and “Sequencing library information”.
- 🟢 Supplementary Note 1 and Supplementary Figure 1 in the paper's main supplementary PDF.

## Evidence trail

- 🟢 The first-strand pool contains 408 individually synthesized mouse NSR hexamers chosen to exclude exact rRNA matches; oligo(dT)18 is also present.
- 🟢 DNase I nicks cDNA in the RNA/cDNA hybrid. RNase-H-minus PrimeScript extends from nicked 3′ ends and displaces downstream cDNA; T4 gene 32 protein promotes displacement and protects displaced products.
- 🟢 Complementary second-strand NSRs and Klenow Fragment (3′→5′ exo−) make double-stranded cDNA.
- 🟢 Scaled Nextera XT tagmentation and 13 PCR cycles produce the original library; the advanced method uses 14 cycles.
- 🟢 Defining datasets used 76-cycle single-read NextSeq 500 sequencing. The protocol has neither cellular barcodes nor UMIs; the authors explicitly identify strand-displacement amplification as preventing their addition.

The 408+408 individual NSR sequences are supplied in Supplementary Data 3 and are not reproduced as an ordering list here.
