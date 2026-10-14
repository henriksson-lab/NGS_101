# Pseudo-seq

**Defining source:** Carlile *et al.*, *Nature* (2014), doi:10.1038/nature13802.
The PMC full text and Methods were read.

## Evidence used

- 🟢 CMC derivatization followed by alkaline reversal leaves Ψ–CMC and a parallel mock
  treatment supplies the background control.
- 🟢 Methods publish CIP/PNK treatment, 80–140-nt size fractions, the modified 3′ adapter,
  RT primer, CircLigase step, both PCR primers and HiSeq 2000 endpoint.
- 🟢 Selected cDNAs terminate at the RT block one nucleotide 3′ of pseudouridine.

## Published library oligos

The open author manuscript's “Pseudo-seq Library Preparation” Methods
(PMCID PMC4224642) prints these oligos. Modifications and the unknown six-base
sample index are retained as written.

| Role | 5′→3′ ordered oligo | Evidence |
|---|---|---|
| 3′ RNA adapter | `/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/` | 🟢 Methods; described there as pre-adenylated, although the printed IDT string says `/5Phos/` |
| RT primer | `/5Phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT/iSp18/CACTCA/iSp18/GCCTTGGCACCCGAGAATTCCA` | 🟢 Methods |
| PCR forward | `AATGATACGGCGACCACCGA` | 🟢 Methods |
| PCR reverse | `CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` | 🟢 Methods; `X` is a variable sample index |

🟡 The paper reports an Illumina HiSeq 2000 run and no separate custom
sequencing-primer oligo. Its PCR primers and adapter establish the library
chemistry, but the current role-level final construct does not yet compute an
exact sequencing-primer binding coordinate.

## Construct logic

Pseudo-seq reads a reaction endpoint rather than the modified base itself. The informative
coordinate is the end of a truncated cDNA, and enrichment over the −CMC control distinguishes
the chemical signal from ordinary RT stops.
