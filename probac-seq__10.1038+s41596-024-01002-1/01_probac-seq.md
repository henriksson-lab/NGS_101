# ProBac-seq — source notes

Defining protocol: Samanta et al., *Nature Protocols* (2024), [doi:10.1038/s41596-024-01002-1](https://doi.org/10.1038/s41596-024-01002-1).

- 🟢 PMC protocol text and oligo tables were read.
- 🟢 Gene-specific probes have a 23-bp PCR/circularization handle, about 50 bp of transcript-complementary sequence and a 17-bp HindIII-containing extender.
- 🟢 Three RCA rounds amplify the pool; HindIII relinearizes it. Isothermal extension adds a random 12-nt probe UMI and poly(A)30.
- 🟢 Fixed bacteria bearing hybridized probes are encapsulated with 10x 3′ v3 beads. Six-cycle in-GEM PCR adds the cell barcode to DNA probes without reverse transcription.
- 🟢 Two later PCRs scale the library and add Illumina adapters. The stated NextSeq 2000 run uses i7 8 cycles and Read 1 119 cycles with 30% PhiX.
