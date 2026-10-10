# DMS-MaPseq — research notes

## Sources read

- 🟢 Zubradt et al. 2017, DOI [10.1038/nmeth.4057](https://doi.org/10.1038/nmeth.4057), PMC5508988, Methods and Supplementary Table 1.

## Protocol boundary

The public page follows the genome-wide branch. The paper also describes targeted
gene-specific RT-PCR, a ten-base-UMI targeted branch, and Nextera XT conversion of
targeted amplicons. Those are not silently merged into this library.

## Molecular order

🟢 Cells are treated with 5% DMS, which modifies exposed A/C bases. Total RNA is
Zn2+-fragmented, size selected, dephosphorylated, depleted of rRNA and ligated at its
3-prime end to preadenylated, 3-prime-blocked linker-2. TGIRT-III reverse transcription
at 57 °C reads through DMS adducts and encodes them as substitutions. RNA is hydrolyzed,
the 5-prime-phosphorylated cDNA is circularized with CircLigase, and 9–13 PCR cycles add
TruSeq-indexed Illumina arms. Custom primer oNTI202 produces a 50-nt single-end read.

## Source-transcribed oligos used by the model

- 🟢 linker-2: `5rApp/CACTCGGGCACCAAGGA/3ddC`
- 🟢 RT handle before the internal spacer: `5'/5phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCG`
- 🟢 custom Read 1 primer oNTI202: `GCAGCGACAGGTTCAGAGTTCTACAGTCCGACGATC`

The full RT oligo also contains an internal spacer and a second priming region. The page
shows the region that fixes circle/PCR orientation and leaves spacer chemistry in these
notes.
