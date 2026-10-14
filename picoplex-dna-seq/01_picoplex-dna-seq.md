# PicoPLEX DNA-seq

## Defining source

🟢 Takara Bio, *PicoPLEX DNA-seq Kit Instruction Manual*, QAM-106-002. This is the named sequencing-ready implementation of PicoPLEX chemistry; SurePlex is a related commercial displacement-DOP-PCR whole-genome-amplification implementation, not silently merged into this page.

## Molecular workflow

🟢 One to ten cells, or 6–60 pg purified DNA, are lysed in the reaction vessel. Proprietary quasi-random primers bind selected genomic sites and support multiple cycles of linear strand-displacement preamplification.

🟢 Exponential PCR then installs a unique pair of 8-nt i5 and i7 indexes. The supplied plate contains 48 unique dual-index combinations and produces an Illumina-compatible library without a separate library-preparation kit. Figure 4 of manual QAM-106-002 prints the i5/i7 index sequences and well map; the current model leaves the chosen plate member variable.

🔴 The quasi-random primers, their common handles and complete amplification/index oligos are proprietary. Consequently the page can state read roles and index lengths but cannot compute sequencing-primer landing positions from published bases.
