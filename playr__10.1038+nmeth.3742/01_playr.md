# PLAYR — source notes

**Defining source:** Frei et al., *Nature Methods* (2016),
[doi:10.1038/nmeth.3742](https://doi.org/10.1038/nmeth.3742).

**Detailed procedure:** Frei et al., *Nature Protocols* (2019),
[doi:10.1038/s41596-018-0120-8](https://doi.org/10.1038/s41596-018-0120-8).

## Molecular architecture

- 🟢 Two gene-specific probes bind adjacent RNA regions. Each carries a poly(A)
  spacer plus tails that recruit one target-specific insert and one universal backbone.
- 🟢 The insert and backbone are 5′ phosphorylated. T4 DNA ligase seals them into a
  reporter circle only when the two RNA probes provide the proximity scaffold.
- 🟢 Probe 2 contributes backbone-binding then insert-binding sequence; Probe 1
  contributes insert-binding then backbone-binding sequence. This reversed order lets
  the insert bridge the inside while the universal backbone closes around the outside.
- 🟢 Phi29 amplification is initiated by the 3′-OH of one member of the probe pair;
  there is no separately added RCA primer.
- 🟡 The schematic uses Probe 1 as the sequence-resolved representative. Its terminal
  reporter tail has one exact site spanning the circular insert/backbone boundary, which
  is drawn before the resulting single-stranded concatemer.
- 🟢 A fluorescent or metal-tagged detection oligo binds the repeated insert sequence,
  making different insert identities separate RNA channels.
- 🟢 Supplementary Table 3 publishes 27 reporter systems. The page shows system 1 as
  a concrete sequence example and leaves the gene-specific RNA-binding arms symbolic.

PLAYR is read by flow cytometry or mass cytometry together with antibody measurements.
It does not create an NGS library.
