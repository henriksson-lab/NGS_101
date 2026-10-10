# M3-seq — research notes

Defining source: McNulty et al., 2023, [doi:10.1038/s41564-023-01462-3](https://doi.org/10.1038/s41564-023-01462-3).

## Sources read

- 🟢 PMC full text and Methods.
- 🟢 Supplementary Figure 1 workflow and supplementary oligo/probe tables.

## Molecular path

🟢 Round-one random RT installs BC1 and UMI in fixed bacteria. Round-two 10x ATAC droplet ligation installs BC2. BC1+BC2 is the cell identity.

🟢 After second-strand synthesis the library is tagmented, PCR-tagged with a T7 promoter, transcribed, hybridized to rRNA probes and RNase-H depleted. A second RT plus PCR creates the sequenced cDNA library.

🟡 Commercial bead bases are not asserted.
