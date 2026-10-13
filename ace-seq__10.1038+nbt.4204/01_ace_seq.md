# ACE-seq — APOBEC-coupled epigenetic sequencing

## Source read

- 🟢 Schutsky *et al.*, “Nondestructive, base-resolution sequencing of 5-hydroxymethylcytosine using a DNA deaminase,” *Nature Biotechnology* (2018), DOI `10.1038/nbt.4204` (defining paper and open author manuscript, PMCID: PMC6453757).

## Molecular interpretation

- 🟢 DNA is sheared to about 300 bp and mixed with methylated lambda and 5hmC-containing T4-phage controls.
- 🟢 T4 β-glucosyltransferase converts 5hmC to glucosyl-5hmC. After heat denaturation and snap cooling, APOBEC3A deaminates C and 5mC while protected 5hmC resists conversion.
- 🟢 PCR and sequencing therefore read original C and 5mC as T, while original 5hmC remains C. No bisulfite treatment is used.
- 🟢 The defining whole-genome experiment prepared converted DNA with the Swift Accel Methyl-NGS kit and used single-end NextSeq sequencing.
- 🟢 Supplementary Table 1 prints four PCR primers for the methylated-lambda and 5hmC-containing T4 spike-in controls. These are shown on the page, but they assay conversion controls; they do not construct or sequence the library.
- 🔴 The proprietary kit does not expose complete library-adapter or sequencing-primer oligos, so those adapter ends and the sequencing entry are shown by molecular role only.

## Final readout

ACE-seq has the same base-call polarity as TAB-seq but achieves it by selective enzymatic deamination after glucosylation, avoiding destructive bisulfite chemistry.
