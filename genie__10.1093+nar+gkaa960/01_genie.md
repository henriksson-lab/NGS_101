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

## Adaptor tails and the two-PCR route (2026-10-17 pass)

🟢 Table S3 prints one tail per orientation, "added to the 5′ end of the forward and
reverse primers respectively": forward `{{= illumina.TRUSEQ_READ1 }}` and reverse
`TCGGCATTCCTGCTGAACCGCTCTTCCGATCT`. The forward tail is exactly the TruSeq Read 1 primer; the reverse tail is the Illumina paired-end /
Genome-Analyzer-era arm used by the cited 96-plex barcoding method (ref 8, Kozarewa &
Turner 2011), **not** the current TruSeq Read 2 arm `{{= illumina.TRUSEQ_READ2 }}`. The
schematic previously drew a generic dual-index TruSeq ligation library, including a dA
junction that an amplicon never has; both are now corrected (a repository modelling
mistake, not a publication error).

🟢 Methods ("Genomic DNA isolation", "RNA isolation and reverse transcription",
"Sequencing") describe exactly two PCRs: PCR1 with the tailed locus primers from gDNA
(MagAttract HMW, 250–500 ng, PowerUp SYBR Green, 95 °C 10 min then 30 × 95/57/62 °C, split
over four tubes to avoid founder bias) or from cDNA (SuperScript IV primed by the
gene-specific RT primer at 0.1 µM, Direct-zol/TRI RNA with in-column plus TURBO DNase),
then PCR2 from 1 µl of PCR1 to add the MiSeq indices. The earlier model's separate "outer"
and "nested" locus PCRs were not in the source and have been removed; the matched
genomic-DNA branch is now drawn.

🟢 The gene-specific RT primer lies 30–50 bp outside the amplicon primers, in the
direction opposite transcription, which the paper credits for nascent-RNA sensitivity.
Amplicons are kept below 295 bp for a 150-base paired-end MiSeq run.

🔴 The barcoding primers of ref 8 are not printed, and neither is which amplicon end
receives P5, so the outer flow-cell/index shell and the drawn orientation are inferred and
index cycles are not shown as acquired reads.
