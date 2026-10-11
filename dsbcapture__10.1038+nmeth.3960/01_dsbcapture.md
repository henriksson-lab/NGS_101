# DSBCapture — research notes

Lensing et al., *Nature Methods* (2016), [doi:10.1038/nmeth.3960](https://doi.org/10.1038/nmeth.3960), PMCID: PMC5045719.

## Sources read

- 🟢 Full author manuscript and Methods from PMC.
- 🟢 Main Figure 1 and Supplementary Figure 1c from the PMC Cloud article package. These print the modified P5/P7 duplexes and the PCR primers.
- 🟢 Supplementary figures and tables from the same package.

## Molecular interpretation

1. 🟢 DSBs in fixed nuclei are blunt-ended and A-tailed. A biotinylated, phosphorothioate-protected modified P5 duplex with a 3′ T overhang is ligated at the native break.
2. 🟢 Lambda exonuclease removes excess adapter before DNA extraction. The genomic DNA is then sonicated to 200–500 bp and biotin-positive fragments are captured on streptavidin beads.
3. 🟢 The sonication-created distal end is blunt-ended, A-tailed and ligated to the modified P7 duplex. Thus P5 is uniquely attached to the native-break side and P7 to the later random end.
4. 🟢 PCR F and an indexed PCR R complete a single-index Illumina library. Figure 1d states that the sequencing primer anneals to the P5 adapter and the first sequenced base is the break site.
5. 🟢 EcoRV pilot libraries were sequenced 2×50 bp on MiSeq; AID-DIvA and NHEK libraries were sequenced 2×75 bp on NextSeq 500. Only proximal Read 1 was retained for break mapping; paired-end data helped distinguish PCR duplicates and compare with BLESS.

## Sequence orientation

The figures draw duplex bottom strands left-to-right as 3′→5′. The ordering list on the public page reverses those bottom-strand characters into actual 5′→3′ synthesis order. It does not copy the diagram text backwards.

`B` in the source denotes biotin, `P` a terminal phosphate and `*` a phosphorothioate bond. PCR R (AD06) is one of three printed example index primers; AD12 and AD19 differ only in the six-base barcode.
