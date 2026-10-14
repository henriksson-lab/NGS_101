# GenIE — source notes

Defining paper: Cooper et al., *Nucleic Acids Research* (2020), DOI
[10.1093/nar/gkaa960](https://doi.org/10.1093/nar/gkaa960).

🟢 GenIE is the paper's name for arrayed editing of endogenous variants followed by
allele-specific targeted cDNA and genomic-DNA amplicon sequencing. The relative cDNA and
genomic-DNA allele counts measure the effect of each edit on expression or splicing.

## Oligo source audit

🟢 The PMC full text, Supplementary Tables S3–S5 in the author-supplied XLSX, and the
supplementary methods PDF were checked. Table S3 prints the common forward and reverse
sequencing-adapter tails and locus-specific gene-specific RT, amplicon forward and reverse
primers for the assayed loci. Tables S4–S5 give the experimental guide and donor oligos.
These are published as locus-specific sets; there is no single universal GenIE RT or
amplicon primer pair. Table S3 also has a separate cDNA forward primer for the SDF4
splicing assay, which should not be silently substituted for the genomic assay primer.

🟡 The current schematic represents a generic cDNA amplicon and its standard indexed
sequencing shell. The source also compares a matched genomic-DNA branch, so the schematic
does not yet cover the complete assay. The paper points to a previously published MiSeq
barcoding method but does not identify one selected index-primer pair for all loci.
