# SNARE-seq — research notes

## Sources read

- 🟢 Chen et al. 2019, DOI `10.1038/s41587-019-0290-0`, PMC full text, Supplementary Protocol and Supplementary Table 4.
- 🟢 The Teichlab SNARE-seq construct page, checked as a secondary reconstruction.

## Construct path

🟢 Nuclei are tagmented before droplets. The phosphorylated Nextera-R1-rc-polyA oligo anneals to Nextera-R1-bk and exposes a poly(A) tract, allowing the tagmented chromatin product to use the same poly(dT) bead capture system as RNA. 🟢 A 12-base cell barcode is therefore shared between RNA and accessibility products; the adjacent 8-base UMI belongs to the RNA readout.

🟢 The RNA branch uses template switching and cDNA amplification. The chromatin branch selectively amplifies tagmented genomic products. Supplementary Table 4 is the authority for the two modified splint oligos and TSO shown on the page.
