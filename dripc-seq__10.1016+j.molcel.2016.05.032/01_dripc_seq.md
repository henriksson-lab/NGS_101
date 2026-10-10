# DRIPc-seq — strand-specific sequencing of R-loop RNA

## Sources read

- 🟢 Sanz *et al.*, “Prevalent, dynamic, and conserved R-loop structures associate with specific epigenomic signatures in mammals,” *Molecular Cell* (2016), DOI `10.1016/j.molcel.2016.05.032` (defining paper and open author manuscript, PMCID: PMC4955522).
- 🟢 Sanz and Chédin, *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1` (stepwise protocol, PMCID: PMC6615061).

## Molecular interpretation

- 🟢 Gentle extraction, five-enzyme restriction digestion and S9.6 capture are shared with DRIP-seq.
- 🟢 DNase I removes the DNA component after capture. The retained R-loop RNA is reverse-transcribed with a mixture of poly(dT) and random hexamers.
- 🟢 Second-strand synthesis substitutes dUTP for dTTP. After cDNA sonication, end repair, dA tailing and indexed-adapter ligation, UNG removes the uracil-containing strand before PCR.
- 🟢 Because only the first cDNA strand survives amplification, read orientation identifies the RNA strand that participated in the original hybrid.
- 🔴 Complete indexed-adapter oligo sequences are not printed; the page therefore shows those ends by role.

## Final readout

The library sequences RNA-derived cDNA rather than the captured genomic DNA. That change provides higher-resolution, strand-specific R-loop mapping and makes DRIPc-seq molecularly distinct from DRIP-seq.
