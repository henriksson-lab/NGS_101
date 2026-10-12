# ProBac-seq — source notes

Defining protocol: Samanta et al., *Nature Protocols* (2024), [doi:10.1038/s41596-024-01002-1](https://doi.org/10.1038/s41596-024-01002-1).

## Sources read

- 🟢 The PMC author manuscript, including the complete stepwise protocol and figure legends.
- 🟢 Supplementary Table 1, containing the ordered *B. subtilis* probe pool.
- 🟢 Supplementary Table 2, containing the RCA, blocked-extension and library-PCR primers.

## Probe preparation

- 🟢 Each ordered protoprobe is 5′-`AAGCTTTGCTCGAAGAAACTATG`–target-specific region–`GTATGGACCGGGTCATC`-3′. The pool members do not all have identical target-region lengths; the public drawing uses the protocol's nominal 50 nt.
- 🟢 RCA-1 and RCA-3 use 5′-`TTCGAGCAAAGCTTGATGACCCGGTCC`-3′. RCA-2 uses its reverse complement, 5′-`GGACCGGGTCATCAAGCTTTGCTCGAA`-3′. Each primer bridges the two probe ends for circularization, primes phi29 and creates the HindIII-containing duplex needed to release unit-length products.
- 🟢 Three successive ligation–RCA–HindIII cycles alternate product orientation. The third releases the strand complementary to bacterial mRNA.
- 🟢 The blocked extension template is 5′-T30–V–N11–`AAGCTTGCTCGAAGAAACTATG`-/3Phos/-3′. Its twelve degenerate positions and T30 tract are copied onto the protoprobe as a 12-nt probe UMI followed by poly(A)30.
- 🟢 The two supplementary tables disagree by one T at this arm: every ordered probe begins `AAGCTTTGCT…`, whereas the extension template ends `AAGCTTGCT…`. The public model retains both transcriptions and does not draw them as an exact fully paired duplex.

## Cellular barcoding and final library

- 🟢 Fixed bacteria bearing hybridized probes are encapsulated with 10x 3′ v3 beads. The released gel-bead oligo supplies partial Read 1, a 16-nt cell barcode, a 12-nt GEM UMI and poly(dT)30. Six-cycle in-GEM PCR copies bound DNA probes without reverse transcription.
- 🟢 PCR-1 and PCR-2 use the 18-nt InDrop reverse primer 5′-`GATGACCCGGTCCATACA`-3′; PCR-2 uses the source-listed 10x forward primer 5′-`CTACACGACGCTCTTCCGATCT`-3′.
- 🟢 PCR-3 adds the full P5/Read 1 arm and a P7/index/Read 2 arm. Supplementary Table 2 provides sixteen distinct eight-base i7 primers.
- 🟢 The stated NextSeq 2000 run acquires Index 1 for 8 cycles and Read 1 for 119 cycles with 30% PhiX. Read 1 contains, in order, the 10x cell barcode, GEM UMI, poly(dT), probe UMI, extension handle and the beginning of the probe identity. The analysis splits this long read and substitutes the probe UMI for the GEM UMI when counting molecules.
- 🟡 The schematic uses the standard 10x v3 gel-bead architecture documented by the same kit named in the protocol; the commercial bead oligo itself is not reprinted in Supplementary Table 2.
