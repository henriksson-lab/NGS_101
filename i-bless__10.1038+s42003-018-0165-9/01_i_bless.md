# i-BLESS — research notes

Biernacka et al., *Communications Biology* (2018), [doi:10.1038/s42003-018-0165-9](https://doi.org/10.1038/s42003-018-0165-9), PMCID: PMC6208412.

## Sources read

- 🟢 Full text and Methods from the PMC article.
- 🟢 Supplementary Information from the PMC Cloud copy (`42003_2018_165_MOESM1_ESM.pdf`), especially Supplementary Table 4.

## What makes this i-BLESS rather than BLESS

1. 🟢 Yeast cells are embedded in agarose beads before spheroplasting, lysis and protein removal. The beads permit reagent diffusion while protecting very small yeast nuclei and long DNA from handling-induced breaks.
2. 🟢 Break ends are blunted in the beads and ligated to the same biotinylated proximal hairpin used by BLESS.
3. 🟢 Encapsulated DNA is first sonicated, recovered from agarose and sonicated again to about 400 bp. Streptavidin capture therefore still separates the pre-labelled native end from a later fragmentation end.
4. 🟢 The captured free end is blunted, phosphorylated and ligated to the BLESS distal hairpin. I-SceI linearization, PCR and XhoI digestion retain the proximal and distal 11-nt tags.
5. 🟢 Final sequencing libraries were prepared with the ThruPLEX DNA-seq Kit and run 2×75 bp on HiSeq 2500 and HiSeq 4000.
6. 🔴 ThruPLEX's complete adapter oligos are not printed in the defining paper. The schematic therefore marks the final Illumina-compatible shell inferred while keeping the exact public i-BLESS insert and tag sequences explicit.

Supplementary Table 4 prints the proximal and distal adapters, their modifications, and the two insert-amplification primers. It explicitly states that the adapters are identical to those of original BLESS.
