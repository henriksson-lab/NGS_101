# RamDA-seq — source notes

Defining paper: Hayashi et al., *Nature Communications* 2018, [doi:10.1038/s41467-018-02866-0](https://doi.org/10.1038/s41467-018-02866-0).

## Sources read

- 🟢 PMC5809388 full-text XML, Methods “RT-RamDA”, “Not-so-random primers”, “RamDA-seq”, and “Sequencing library information”.
- 🟢 Supplementary Note 1 and Supplementary Figure 1 in the paper's main supplementary PDF.
- 🟢 Supplementary Data 3 (`s41467-018-02866-0_supp_41467_2018_2866_MOESM6_ESM.xlsx`),
  the complete first- and second-strand mouse NSR ordering lists.

## Evidence trail

- 🟢 The first-strand pool contains 408 individually synthesized mouse NSR hexamers chosen to exclude exact rRNA matches; oligo(dT)18 is also present.
- 🟢 A representative matched pair is `1st_MmNSR_1` = `AAAGCA` and
  `2nd_MmNSR_1` = `TGCTTT`; the latter is the reverse complement in the published
  Supplementary Data 3 table. The page shows this pair as an ordering example,
  not as a replacement for the 408-member pools.
- 🟢 DNase I nicks cDNA in the RNA/cDNA hybrid. RNase-H-minus PrimeScript extends from nicked 3′ ends and displaces downstream cDNA; T4 gene 32 protein promotes displacement and protects displaced products.
- 🟢 Complementary second-strand NSRs and Klenow Fragment (3′→5′ exo−) make double-stranded cDNA.
- 🟢 Scaled Nextera XT tagmentation and 13 PCR cycles produce the original library; the advanced method uses 14 cycles.
- 🟢 Defining datasets used 76-cycle single-read NextSeq 500 sequencing. The protocol has neither cellular barcodes nor UMIs; the authors explicitly identify strand-displacement amplification as preventing their addition.

The 408+408 individual NSR sequences are supplied in Supplementary Data 3; the
full pool is not reproduced as an ordering list here.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5809388 JATS180–204, standard/advanced/C1 RamDA preparation and sequencing.

Added initial DNase I genomic-DNA digestion before RT and explicit unresolved fragments. Specified 13-cycle standard versus 14-cycle advanced library PCR; candidate Tn5 products no longer presumed discarded. Confirmed single-read76 acquisition.

Deferred: Original NSR tables, pre-/post-tagmentation cleanup, concurrent nicking/RT representation; C1 variant is separate and its JATS 130°C wording needs PDF verification before interpretation. This pass does not establish a complete source review or publication-wide absence.
