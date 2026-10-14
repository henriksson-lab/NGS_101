# i-BLESS — research notes

Biernacka et al., *Communications Biology* (2018), [doi:10.1038/s42003-018-0165-9](https://doi.org/10.1038/s42003-018-0165-9), PMCID: PMC6208412.

## Sources read

- 🟢 Full text and Methods from the PMC article.
- 🟢 Supplementary Information from the PMC Cloud copy (`42003_2018_165_MOESM1_ESM.pdf`), especially Supplementary Table 4.
- 🟢 Rubicon Genomics, *ThruPLEX DNA-seq Kit Manual*, QAM-108-003 (revised September 2016), Appendix 1B and Fig. 9, [archived vendor PDF](https://www.takara.co.kr/file/manual/pdf/ThruPLEX-DNA-seq-Kit-Manual-QAM-108-003-updated-9_16_16.pdf).
- 🟢 [Takara's later ThruPLEX adapter FAQ](https://help.takarabio.com/36168/kb/article/110112/what-are-the-index-adapter-sequences-after-preparation-using-thruplex-dna-libraries) was checked for final-library context, but its newer kit description does not establish the ordered oligos in the 2018 study.

## What makes this i-BLESS rather than BLESS

1. 🟢 Yeast cells are embedded in agarose beads before spheroplasting, lysis and protein removal. The beads permit reagent diffusion while protecting very small yeast nuclei and long DNA from handling-induced breaks.
2. 🟢 Break ends are blunted in the beads and ligated to the same biotinylated proximal hairpin used by BLESS.
3. 🟢 Encapsulated DNA is first sonicated, recovered from agarose and sonicated again to about 400 bp. Streptavidin capture therefore still separates the pre-labelled native end from a later fragmentation end.
4. 🟢 The captured free end is blunted, phosphorylated and ligated to the BLESS distal hairpin. I-SceI linearization, PCR and XhoI digestion retain the proximal and distal 11-nt tags.
5. 🟢 Final sequencing libraries were prepared with the ThruPLEX DNA-seq Kit and run 2×75 bp on HiSeq 2500 and HiSeq 4000.
6. 🟢 The 2016 kit guide says repaired DNA receives a blunt-ended, 5′-blocked stem-loop adapter; extension and cleavage of that adapter occur during PCR. Its single-index kit has an **8-nt Sanger index**. Appendix 1B gives all 12 index *bases* (for example, Tube 1 `ATCACGTT`), and Fig. 9 places the index at the library's 3′ end.
7. 🔴 The paper does not identify the index tube used for any i-BLESS sample. The historical kit guide does not print the ordered stem-loop adapter or complete amplification primer oligos. The later Takara FAQ describes a newer final-library architecture, but cannot establish those 2018 reagent sequences. The schematic therefore marks the final adapter shell inferred and keeps its index as an 8-nt placeholder, while retaining the exact public i-BLESS insert and tag sequences.

Supplementary Table 4 prints the proximal and distal adapters, their modifications, and the two insert-amplification primers. It explicitly states that the adapters are identical to those of original BLESS.
