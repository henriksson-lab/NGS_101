# Named pooled-reporter protocol backlog

Candidates adjacent to TRIP: pooled assays in which sequencing links a reporter, tested
sequence, perturbation or genomic context to expression or molecular abundance. This is a
review list, not a commitment to add every method. It deliberately excludes papers that do
not give their workflow a name, in accordance with `DESIGN.md`.

Catalogue coverage was checked against `catalogue/ours.tsv` on 2026-10-12. None of the
candidates below is currently covered.

## Closest to TRIP and general reporter assays

| Protocol | Defining source | Distinct library-preparation or readout trick | In `ours.tsv`? |
|---|---|---|---|
| MPRA | [Melnikov et al., 2014, *Journal of Visualized Experiments*](https://doi.org/10.3791/51719) | Synthesized candidate regulatory sequences are coupled to transcribed sequence tags; tag counts from reporter cDNA are compared with the input plasmid library. This is the useful canonical barcoded, episomal MPRA page. | No |
| STARR-seq | [Arnold et al., 2013, *Science*](https://doi.org/10.1126/science.1232542) | The candidate enhancer is cloned into the reporter transcript's 3′ UTR and therefore identifies itself in RNA; no separate candidate-to-barcode lookup is needed. | No |
| UMI-STARR-seq | [Neumayr et al., 2019, *Current Protocols in Molecular Biology*](https://doi.org/10.1002/cpmb.105) | After reporter-specific second-strand synthesis, a linear PCR introduces a UMI in place of the i7 index before any exponential amplification; the UMI is read as Index 1. | No |
| lentiMPRA | [Gordon et al., 2020, *Nature Protocols*](https://doi.org/10.1038/s41596-020-0333-5) | Candidate elements and random barcodes are delivered by lentivirus; the 5′-UTR barcode is counted in RNA and matched DNA libraries to normalize integrated reporter abundance. | No |
| LV-MPRA | [Maricque et al., 2017, *Nucleic Acids Research*](https://doi.org/10.1093/nar/gkw942) | A distinct early lentiviral MPRA implementation places reporters into chromosomal context, allowing integrated activity to be contrasted with episomal activity. | No |
| SuRE | [van Arensbergen et al., 2017, *Nature Biotechnology*](https://doi.org/10.1038/nbt.3754) | Random 0.2–2-kb genomic fragments are cloned upstream of a 20-bp reporter barcode; paired-end library decoding associates fragment and barcode before cDNA barcode counting. | No |
| TRACE | [Tchasovnikarova et al., 2022, *Molecular Cell*](https://doi.org/10.1016/j.molcel.2021.11.035) | A barcode in the lentiviral U3 region is copied into both LTRs after integration, leaving one copy close enough to flanking genomic DNA for restriction/ligation-mediated site mapping and another in the reporter transcript for expression measurement. | No |
| scMPRA | [Zhao et al., 2023, *Nature Genetics*](https://doi.org/10.1038/s41588-022-01278-7) | Each reporter carries a candidate barcode and a high-complexity molecule barcode; reporter RNAs acquire a shared 10x cell barcode, linking regulatory activity to cell identity. A U6-driven second copy identifies constructs even when the tested promoter is silent. | No |
| CiBER-seq | [Lobel et al., 2025, *Genome Biology*](https://doi.org/10.1186/s13059-025-03610-w) | CRISPRi guides are linked to matched reporter and normalizer RNA barcodes in a single-copy genomic cassette, so bulk barcode sequencing reports a molecular phenotype rather than guide abundance alone. | No |

## Named STARR/MPRA variants with distinct chemistry

| Protocol | Defining source | Distinct library-preparation or readout trick | In `ours.tsv`? |
|---|---|---|---|
| CapStarr-seq | [Vanhille et al., 2015, *Nature Communications*](https://doi.org/10.1038/ncomms7905) | Hybridization capture first enriches chosen genomic regions; captured fragments become the self-transcribing inserts, avoiding oligo synthesis while restricting library complexity. | No |
| STAP-seq | [Arnold et al., 2017, *Nature Biotechnology*](https://doi.org/10.1038/nbt.3739) | Candidate fragments occupy the core-promoter position rather than the enhancer/3′-UTR position. Their exact RNA 5′ ends measure transcription initiation and enhancer responsiveness. | No |
| ChIP-STARR-seq | [Barakat et al., 2018, *Cell Stem Cell*](https://doi.org/10.1016/j.stem.2018.06.014) | ChIP-enriched DNA is cloned as the self-transcribing candidate pool, directly testing the regulatory activity of fragments selected by transcription-factor binding or histone marks. | No |
| mSTARR-seq | [Lea et al., 2018, *eLife*](https://doi.org/10.7554/eLife.37513) | Genomic fragments are inserted into a CpG-free STARR vector, then the same plasmid library is split into M.SssI-methylated and sham-treated branches before DNA/RNA sequencing. | No |
| ATAC-STARR-seq | [Hansen and Hodges, 2022, *Genome Research*](https://doi.org/10.1101/gr.276766.122) | Tn5-isolated accessible fragments are cloned directly into the STARR reporter. Sequencing the matched plasmid DNA and reporter RNA retains accessibility and footprint information as well as regulatory activity. | No |

## Adjacent named reporter protocols

These are not enhancer assays, but their construct and readout geometry belongs to the same
barcoded-reporter design space and would exercise reusable models beyond TRIP.

| Protocol | Defining source | Distinct library-preparation or readout trick | In `ours.tsv`? |
|---|---|---|---|
| PTRE-seq | [Cottrell et al., 2018, *Nature Communications*](https://doi.org/10.1038/s41467-017-02745-0) | Synthetic 3′-UTR regulatory elements are represented by multiple internal barcodes; barcode libraries are prepared from total RNA and ribosomal/polysomal fractions to separate abundance, initiation and translation effects. | No |
| Vex-seq | [Adamson et al., 2018, *Genome Biology*](https://doi.org/10.1186/s13059-018-1437-x) | Test exons plus flanking intronic sequence are linked to a transcript barcode; a 10-nt UMI is added during RT, and sequencing compares variant representation before and after splicing. | No |

## Suggested order for later work

1. TRACE, because it is the most direct named successor/relative of TRIP and adds lentiviral
   LTR barcode duplication and a different integration-site recovery geometry.
2. MPRA, STARR-seq and lentiMPRA, which establish the three major reusable reporter models:
   separate barcode, self-transcribing insert and integrated barcode.
3. SuRE, UMI-STARR-seq and mSTARR-seq, which add fragment–barcode decoding, pre-amplification
   molecule counting and a methylated/unmethylated branch.
4. scMPRA and CiBER-seq, which connect reporter barcodes to cell identity or perturbation.
5. The input-selection variants (CapStarr-seq, ChIP-STARR-seq and ATAC-STARR-seq), STAP-seq,
   then the post-transcriptional reporters PTRE-seq and Vex-seq.
