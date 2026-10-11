# Pseudo-seq

**Defining source:** Carlile *et al.*, *Nature* (2014), doi:10.1038/nature13802.
The PMC full text and Methods were read.

## Evidence used

- 🟢 CMC derivatization followed by alkaline reversal leaves Ψ–CMC and a parallel mock
  treatment supplies the background control.
- 🟢 Methods publish CIP/PNK treatment, 80–140-nt size fractions, the modified 3′ adapter,
  RT primer, CircLigase step, both PCR primers and HiSeq 2000 endpoint.
- 🟢 Selected cDNAs terminate at the RT block one nucleotide 3′ of pseudouridine.

## Construct logic

Pseudo-seq reads a reaction endpoint rather than the modified base itself. The informative
coordinate is the end of a truncated cDNA, and enrichment over the −CMC control distinguishes
the chemical signal from ordinary RT stops.
