# DRIPc-seq — strand-specific sequencing of R-loop RNA

## Sources read

- 🟢 Sanz *et al.*, “Prevalent, dynamic, and conserved R-loop structures associate with specific epigenomic signatures in mammals,” *Molecular Cell* (2016), DOI `10.1016/j.molcel.2016.05.032` (defining paper and open author manuscript, PMCID: PMC4955522).
- 🟢 Sanz and Chédin, *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1` (stepwise protocol, PMCID: PMC6615061).
- 🟢 Its open XML text (`PMC6615061_PMC6615061.xml.txt`) was checked at the reagents list (lines 606–621) and the example adapter-trimming FASTA (lines 2310–2325).

## Molecular interpretation

- 🟢 Gentle extraction, five-enzyme restriction digestion and S9.6 capture are shared with DRIP-seq.
- 🟢 DNase I removes the DNA component after capture. The retained R-loop RNA is reverse-transcribed with a mixture of poly(dT) and random hexamers.
- 🟢 Second-strand synthesis substitutes dUTP for dTTP. After cDNA sonication, end repair, dA tailing and indexed-adapter ligation, UNG removes the uracil-containing strand before PCR.
- 🟢 Because only the first cDNA strand survives amplification, read orientation identifies the RNA strand that participated in the original hybrid.
- 🟢 The detailed protocol specifies Illumina TruSeq Single index adapters and gives library-amplification primers P5 `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′` and P7 `5′-CAAGCAGAAGACGGCATACGAGAT-3′`. Its analysis example gives adapter-trimming motifs `ACACTCTTTCCCTACACGACGCTGTTCCATCT` and `CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT`; these are trimming references, not evidence for the complete indexed ligation oligos or their physical end configuration.
- 🔴 The exact selected TruSeq Single index oligos and the instrument's sequencing primer are not printed in the defining paper or detailed protocol; the page therefore still shows the indexed ligation ends by role.

## Final readout

The library sequences RNA-derived cDNA rather than the captured genomic DNA. That change provides higher-resolution, strand-specific R-loop mapping and makes DRIPc-seq molecularly distinct from DRIP-seq.
