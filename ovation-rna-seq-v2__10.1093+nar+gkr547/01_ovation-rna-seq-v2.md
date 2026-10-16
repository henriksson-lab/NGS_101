# Ovation RNA-Seq V2 / Ribo-SPIA

## Defining sources

🟢 Tariq *et al.*, “Whole-transcriptome RNAseq analysis from minute amount of total RNA”, *Nucleic Acids Research* (2011), [doi:10.1093/nar/gkr547](https://doi.org/10.1093/nar/gkr547), is the representative published sequencing application. The NuGEN Ovation RNA-Seq System V2 user guide is the commercial reaction authority.

## Molecular workflow

🟢 Chimeric DNA/RNA primers initiate first- and second-strand synthesis from total RNA. The second-strand product carries an RNA-containing SPIA priming end.

🟢 RNase H cleaves the RNA portion, a new SPIA primer anneals, and a strand-displacing DNA polymerase releases a linear single-stranded cDNA copy. Repeated cycles give linear Ribo-SPIA amplification rather than PCR.

🟢 The defining study tested sheared and unsheared amplified cDNA. Its separate
Illumina library branch used a New England Biolabs library-preparation kit,
then paired 50-nt sequencing. The paper does not name an Encore Rapid kit.

🔴 The audited article does not print the chimeric primer sequences or complete
adapter/index oligos. The cited NuGEN V2 guide, the 2005 SPIA predecessor and
the exact NEB adapter kit still need inspection before absence can be accepted.
The page marks its conventional final Illumina shell as inferred.
