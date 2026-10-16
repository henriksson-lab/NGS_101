# Tri-C — research notes

## Citation

Oudelaar AM et al. *Single-allele chromatin interactions identify regulatory hubs in dynamic compartmentalized domains.* Nature Genetics 50, 1744–1751 (2018). [doi:10.1038/s41588-018-0253-2](https://doi.org/10.1038/s41588-018-0253-2)

## Sources read

- 🟢 The defining author manuscript and Methods in [PMC6265079](https://pmc.ncbi.nlm.nih.gov/articles/PMC6265079/).
- 🟢 The later Capture-C protocol's Tri-C module in [PMC7613269](https://pmc.ncbi.nlm.nih.gov/articles/PMC7613269/), used as corroboration.
- 🟢 The defining paper's Supplementary Tables 5 and 7 (publisher supplementary PDF), its PMC supplementary figure/table files, and the earlier bioRxiv version were checked for oligos. The cited [NG Capture-C predecessor](https://doi.org/10.1038/nmeth.3664) was checked for inherited library construction.
- 🟢 The defining article also cites the author's [2018 Tri-C Protocol Exchange method](https://doi.org/10.1038/protex.2018.113); its [2019 revision](https://doi.org/10.21203/rs.2.1650/v2) was checked independently for reagent identity and reaction order.

## 1. What it is

🟢 Tri-C targets multiway rather than only pairwise chromatin contacts. A single physical 3C concatemer must retain the viewpoint and at least two reporter fragments, so its two proximity-ligation junctions represent a concurrent single-allele contact.

## 2. Construct-changing steps

1. 🟢 Fix cells in 2% formaldehyde, digest with NlaIII and proximity-ligate in nuclei.
2. 🟢 NlaIII was selected because it produces small fragments; the defining experiment reports a median mouse-genome fragment size of 131 bp.
3. 🟢 Reverse crosslinks and sonicate to about 450 bp. About half of these fragments can contain multiple ligation junctions.
4. 🟢 Install Illumina TruSeq adapters with NEBNext DNA Library Prep reagents and perform seven cycles of indexing PCR.
5. 🟢 Perform double oligonucleotide capture with a biotinylated probe targeted to the middle of the short viewpoint restriction fragment.

## 3. Final library and reads

🟢 The defining experiments used 250-bp or 150-bp paired-end Illumina sequencing. The molecular unit is viewpoint–reporter 1–reporter 2 (or more), not three independently observed pairwise reads.

🟢 Random sonication ends serve as unique molecular identifiers for PCR duplicate removal.

🟡 Capture-probe bases are locus-specific and therefore remain role tokens on the generic page.

🟢 Published Supplementary Table 5 supplies the genomic coordinates and 120-bp length of four Tri-C capture oligos, but does not print their base strings. Supplementary Table 7 prints two PCR primers for the separate C-Trap Nanopore validation assay; they are not Tri-C library or capture primers. The preprint numbers these as Supplementary Tables 3 and 5.

🔴 The selected 2018 TruSeq/NEBNext adapter product and sample-index pair, and the exact adapter-end ligation geometry, still need a matching historical vendor guide. The linked NG Capture-C method gives the reagent families and end-repair/dA-tailing order, but does not identify a unique Tri-C lot or index pair.

🟢 The 2018 Protocol Exchange method specifies NEBNext E6040S/L, indexed oligo sets E7500S/L and E7335S/L, end repair, dA-tailing, adapter ligation, USER treatment, then Herculase II indexing PCR. These kit identities narrow the missing vendor evidence but do not disclose the selected index oligo in each experimental library or the adapter duplex bases. The later revision repeats these reagent families and introduces Ultra II as a low-input alternative; that alternative is not evidence for the original experiment.

🟡 **Published probe-concentration discrepancy.** The 2018 Nature Genetics Methods (Tri-C experimental procedure) states 2.9 nM for the *total* equimolar pool; its directly cited 2018 Protocol Exchange method (Capture enrichment, “Biotinylated capture oligonucleotides”) specifies a working solution with *each individual oligonucleotide* at 2.9 nM. The earlier bioRxiv v1 Methods instead calls the pool 13 fmol. The latter is an amount, so it cannot be compared with concentration without the volume and number of probes. The 2019 revised method also says 2.9 nM per oligo. The model uses no concentration; the 2018 experiment's exact mixture remains uncertain.

🟡 **Published kit and dilution discrepancies.** The 2018 Protocol Exchange materials list (p. 2; extracted line 50) specifies NEBNext E6040S/L, while its library-preparation heading (p. 5; line 135) calls the procedure E7370. The 2019 revision's materials list (p. 3; line 71) and procedure heading (p. 9; lines 223–227) both specify E6040; it separately offers E7645 Ultra II for limited input (p. 9; lines 230–235). We use the 2019 revision to identify the internally consistent E6040 procedure, without assuming that every 2018 library used one kit lot. Both protocol versions also say that dissolving 10 nmol of probe in 3.46 µl yields 2.9 µM (2018 p. 7; lines 228–233; 2019 p. 14; lines 383–389). That recipe computes to about 2.9 mM, a thousandfold higher concentration. No numeric probe concentration is used in the model pending clarification.
