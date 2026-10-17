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
- 🟢 Nature Methods Supplementary Table 1 (Springer MOESM647, sheet "PLAYR sequences")
  and Nature Protocols Table 1 publish 27 reporter systems. The page shows system 1 as
  a concrete sequence example and leaves the gene-specific RNA-binding arms symbolic.
  (Model correction 2026-10-17: an earlier note called this "Supplementary Table 3".)
- 🟢 System 1 as written in both tables: Probe 1 tail `AAAAAAAAAACTCAGTCGTGACACTCTT`,
  Probe 2 tail `AAAAAAAAAAGACGCTAATATCGTGACC` (target-specific sequence is added 5′ of
  the poly(A) spacer), insert `P-ACGACTGAGTTTGGTCACGAT`, universal backbone
  `P-ATTAGCGTCCAGTGAATGCGAGTCCGTCTAGGAGAGTAGTACAGCAGCCGTCAAGAGTGTC`, detection
  `Z-ACGACTGAGTTTGGTCACGAT`. P = 5′-phosphate; Z = 5′-C6-thiol (metal X8 polymer
  conjugate, mass cytometry) or fluorophore (flow cytometry) per the Nature Protocols
  Table 1 footnote. MOESM647 also lists a universal `Backbone Detection`
  `Z-CAGTGAATGCGAGTCCGTCT` (a backbone substring, so it binds the RCA product); it is
  not drawn.

## Reaction order and conditions

🟢 Frei 2016 Online Methods "PLAYR protocol"; Nature Protocols Reagent setup and
Procedure "PLAYR and intracellular antibody staining" (Steps 13–32):

1. PFA fixation (1.6%, 10 min) and ice-cold methanol permeabilization; optional
   antibody staining with BS3 crosslinking before methanol.
2. Probes pre-heated 90 °C 5 min, chilled, hybridized at 100 nM each in 1× SSC/PVSA/RVC/
   salmon-sperm-DNA buffer, 1 h 40 °C; washes; 4× SSC post-hybridization wash 20 min 40 °C
   (optional Pd barcoding for mass cytometry around this step).
3. Insert/backbone hybridization, 100 nM each in PBS(T) 1× SSC, 37 °C 30 min; two washes.
4. T4 DNA ligase (Thermo), 30 min room temperature (Nature Protocols: 21 °C, 5 U/mL);
   two washes.
5. Phi29 DNA polymerase (Thermo), 30 °C, 2 h (flow) or 4 h (mass), up to 16 h/overnight.
6. Detection oligo 30 min 37 °C. Nature Methods: 5 nM (flow) or 10 nM metal-conjugated
   (mass); Nature Protocols reagent setup: 1 nM of each metal-conjugated probe. This is a
   later protocol revision, not modelled; the page draws no concentrations.

The page draws insert/backbone hybridization and ligation as separate states
(model correction 2026-10-17; they were previously merged into one step).

PLAYR is read by flow cytometry or mass cytometry together with antibody measurements.
It does not create an NGS library.
