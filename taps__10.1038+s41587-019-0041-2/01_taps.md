# TAPS

## 1. What it is

TAPS is a bisulfite-free methylation-sequencing method that converts 5mC and 5hmC, rather than unmodified cytosine, into bases read as thymine.

## 2. Authoritative source

🟢 Liu et al., *Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution*, Nature Biotechnology (2019), DOI `10.1038/s41587-019-0041-2`.

## 3. Construct-changing operations

🟢 TET oxidation produces 5caC; pyridine-borane reduction produces dihydrouracil; PCR reads that product as thymine while leaving unmodified cytosines unchanged.

## 4. Library and read layout

🟢 Supplementary Methods section 3.3 prints a complete single-index TruSeq model library, including TruSeq Universal Adapter, Index 6 (`GCCAAT` as carried in the library), and the final product sequence. The uracil-containing linker used to make the synthetic model is removed with USER after ligation.
🟢 The supplement's separate **genomic sequencing-library procedure** specifies the ligation pair `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` and `/5Phos/GATCGGAAGAGCACACGTCT` without phosphorothioate bonds, followed by indexed PCR using NEBNext Multiplex Oligos for Illumina (Index Set 1, E7335S). It does not identify which index primer was used for a particular genomic sample. The schematic therefore keeps i7 as a source-limited role, rather than carrying over the synthetic model's fixed Index 6.

🟢 The endpoint contains standard TruSeq Read 1, Index 1 and Read 2 primer sites. No i5 index is present in the printed construct.
