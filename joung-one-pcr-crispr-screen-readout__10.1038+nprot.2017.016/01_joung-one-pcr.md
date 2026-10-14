# Joung one-PCR pooled CRISPR screen readout

Joung et al. amplify an integrated sgRNA cassette directly into a sequencing-ready single-index Illumina library in one PCR. The forward primer binds hU6; the reverse primer binds the invariant sgRNA scaffold, making the geometry independent of cPPT placement.

## Source read

- 🟢 Joung et al. 2017, *Nature Protocols*, [doi:10.1038/nprot.2017.016](https://doi.org/10.1038/nprot.2017.016), especially the NGS-Lib-Fwd and NGS-Lib-Rev primer sets and the 80-cycle Read 1 plus 8-cycle index run.

🟢 Table 3 prints **ten staggered NGS-Lib-Fwd primers** and says all ten are mixed for
guide-library sequencing. It also prints **eight indexed NGS-Lib-KO-Rev primers** for
the GeCKO branch; the number used depends on sample multiplexing. The molecular page
shows one published forward/reverse member as a representative molecule, not the whole
primer mix. The SAM reverse family in the same table binds a different scaffold and is
outside this GeCKO readout.

## Distinguishing feature

One PCR simultaneously enriches the integrated guide and adds the complete P5/P7 library arms. It is not interchangeable with workflows where PCR1 first enriches the locus and PCR2 later adds indexed adapters.
