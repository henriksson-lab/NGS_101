# UMI-4C — research notes

## Citation

Schwartzman O et al. *UMI-4C for quantitative and targeted chromosomal contact profiling.* Nature Methods 13, 685–691 (2016). [doi:10.1038/nmeth.3922](https://doi.org/10.1038/nmeth.3922)

## Sources read

- 🟢 The defining paper's publisher record and abstract.
- 🟢 The authors' deposited experimental description for [GSM2054403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2054403).
- 🟢 The defining paper's Supplementary Table 4 (MOESM330) prints the forked adapter, universal reverse primer and locus-specific example bait primers; its Supplementary Figure 1 gives library-construction context. The source bundle and oligo scraper were checked against the original supplement.

## 1. What it is

🟢 UMI-4C makes quantitative, viewpoint-targeted profiles from a linear 3C library. It uses random sonication endpoints as molecule-specific identifiers: the “UMI” is a physical genomic coordinate, not a synthesized N-base barcode.

## 2. Construct-changing steps

1. 🟢 Crosslink, perform high-concentration DpnII digestion and proximity ligation, then reverse crosslinks.
2. 🟢 Sonicate the 3C DNA to roughly 700–900 bp, end-repair, dA-tail and dephosphorylate 5′ ends.
3. 🟢 Ligate indexed Illumina adapters to one end, heat-denature and select the strand carrying the ligated adapter.
4. 🟢 Amplify with two nested viewpoint PCRs and a universal adapter primer. The second viewpoint primer adds the remaining Illumina end.

## 3. Final library and reads

🟢 Paired-end sequencing associates the viewpoint/reporting contact with its random sonication endpoint. Reads with the same contact and endpoint can be collapsed as PCR duplicates.

🟡 The endpoint is represented with structured molecular-identifier metadata so it is searchable as a UMI while its note preserves that it is coordinate-derived.

🟢 Supplementary Table 4 prints both forked-adapter oligos: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` and `(PO4)-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC(NNNNNN)ATCTCGTATGCCGTCTTCTGCTTG`. The six-base window belongs to the adapter sample index; it is distinct from the genomic sonication-endpoint UMI. The table also prints Universal Primer 2 (`CAAGCAGAAGACGGCATACGA`) and example viewpoint primers, including G1p1 upstream `CATCCCCTCAACTCCACCTTAC` and downstream `CTCAACTCCACCTTACCCACCAT` (with the printed Illumina prefix).

🟡 Bait sequences still depend on the chosen locus. The page shows source-printed G1p1 examples while the generic contact construct keeps the viewpoint as a role. The precise annealed fork and dA-tail ligation junction are not yet verified from a source figure, so that graph geometry remains unresolved.
