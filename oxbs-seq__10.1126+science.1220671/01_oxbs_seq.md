# oxBS-seq — oxidative bisulfite sequencing

## Sources read

- 🟢 Booth *et al.*, “Quantitative sequencing of 5-methylcytosine and 5-hydroxymethylcytosine at single-base resolution,” *Science* (2012), DOI `10.1126/science.1220671` (defining paper).
- 🟢 Booth *et al.*, “Oxidative bisulfite sequencing of 5-methylcytosine and 5-hydroxymethylcytosine,” *Nature Protocols* (2013), DOI `10.1038/nprot.2013.115` (open author manuscript, PMCID: PMC3919000).

## Molecular interpretation

- 🟢 Genomic DNA is sonicated, supplied with conversion controls and prepared with the Illumina TruSeq DNA kit before conversion.
- 🟢 Potassium perruthenate selectively oxidizes 5hmC to 5fC. Bisulfite then converts unmodified C and 5fC-derived material to U, while 5mC remains protected.
- 🟢 Uracil-tolerant PCR copies original C and 5hmC as T; original 5mC reads as C.
- 🟢 A matched ordinary BS-seq aliquot reads 5mC + 5hmC as C. 5hmC abundance is calculated from the BS-minus-oxBS difference; it is not directly distinguished within one oxBS molecule.

## Final readout

The final construct is a single-index TruSeq bisulfite library. Paired or single-end genomic reads carry the converted-base signal; the index read identifies the library.
