# DMS-MaPseq — research notes

## Sources read

- 🟢 Zubradt et al. 2017, DOI [10.1038/nmeth.4057](https://doi.org/10.1038/nmeth.4057), PMC5508988, Methods and Supplementary Table 1 in the publisher's Supplementary Text and Tables PDF.

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
- 🟢 oCJ200-link2 RT primer: `5'/5phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCG/iSp18/CAAGCAGAAGACGGCATACGAGATTCCTTGGTGCCCGAGTG`
- 🟢 library PCR primer oNTI231: `CAAGCAGAAGACGGCATACGA`
- 🟢 indexing PCR primer: `AATGATACGGCGACCACCGAGATCTACACGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNCGACAGGTTCAGAGTTC`
- 🟢 custom Read 1 primer oNTI202: `CGACAGGTTCAGAGTTCTACAGTCCGACGATC`

The internal `/iSp18/` is a chemical spacer, not eighteen bases. The model keeps it as a
named junction in the ordering record while constructing the PCR template only from
polymerase-readable sequence. The six `N` positions in the indexing primer are the
source-declared sample index, not unknown adapter bases.
