# DMS-seq

**Defining source:** Rouskin *et al.*, *Nature* (2014), doi:10.1038/nature12894.
The PMC full text and Methods were read.

## Evidence used

- 🟢 DMS modifies solvent-accessible A and C and blocks reverse transcription.
- 🟢 Methods specify Zn²⁺ fragmentation, 60–80-nt RNA selection, PNK end repair,
  modified linker-1, 25–50-nt stopped-cDNA selection, CircLigase, 8–10 PCR cycles and HiSeq 2000.
- 🟢 In-vivo, refolded in-vitro and denatured controls distinguish cellular structure and bias.
- 🟡 PCR and sequencing-primer sequences remain role-level because the paper does not print them.

## Construct logic

The library preserves a termination coordinate rather than a full RNA fragment. A read end
immediately before A or C reports a DMS adduct; relative counts measure accessibility.
