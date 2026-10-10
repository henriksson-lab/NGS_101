# Bulk UMI 5′-RACE AIRR-seq

## Source read

🟢 Gupta et al., *Bulk Sequencing from mRNA with UMI for Evaluation of B-Cell Isotype and Clonal Evolution: A Method by the AIRR Community*, Methods in Molecular Biology (2022), [doi:10.1007/978-1-0716-2115-8_19](https://doi.org/10.1007/978-1-0716-2115-8_19), including Figs. 2 and 5 and the raw-read processing section.

## Molecular path

🟢 Oligo(dT)-primed SMARTScribe RT copies polyadenylated BCR mRNA. Non-templated bases at the completed cDNA end allow the SMART UMI oligo to template-switch, installing a universal PCR handle and a 12-nt UMI on each parental cDNA.

🟢 PCR1 pairs the universal forward primer with separate IgG, IgM, IGK and IGL constant-region reverse primers. PCR2 is semi-nested and adds Illumina Read 1/Read 2 arms plus i5 and i7 sample indexes. This preserves the complete V(D)J region and enough constant region for isotype calling.

## Read layout and uncertainty

🟢 The protocol states that Read 1 starts in the constant region and Read 2 starts upstream of the V region. Read 2 begins with a 12-nt UMI followed by a seven-base linker/template-switch sequence; the first 19 nt may be removed when UMI analysis is not performed.

🔴 The kit primer and SMART UMI oligo sequences are not public. The model therefore uses named, length-preserving role placeholders and does not present them as orderable oligos or claim exact primer-binding geometry. The 12-nt UMI, seven-base linker, dual indexes and their read-side order are published rather than inferred.
