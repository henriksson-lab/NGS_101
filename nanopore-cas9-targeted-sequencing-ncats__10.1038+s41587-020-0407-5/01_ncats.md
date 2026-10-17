# nCATS — research notes

## Defining source

Gilpatrick T, Lee I, Graham JE, et al. *Targeted nanopore sequencing with Cas9-guided adapter ligation.* Nature Biotechnology 2020. [doi:10.1038/s41587-020-0407-5](https://doi.org/10.1038/s41587-020-0407-5), PMCID PMC7145730.

Sources read (bounded re-read 2026-10-17: JATS Results line 58, Online Methods Guide RNA design to Sequencing lines 102–118; protocols.io V.2 steps 7–21): primary-paper Results, Figure 1 legend and Online Methods (“Cas9 enrichment and library preparation”); PMC7145730 Supplementary Table S1 (`NIHMS1548342-supplement-Sup_Tab1.xlsx`) for the locus-specific guide RNA names, 20-base guide sequences and target positions; the linked author protocols.io Cas9-enrichment protocol. Table S1 is the authority for selecting the relevant guide pair and is separate from the undisclosed ONT adapter.

## Evidence record

- 🟢 About 3 µg HMW genomic DNA is dephosphorylated with Quick CIP for 10 min at 37 °C, then heated for 2 min at 80 °C. Existing DNA ends therefore cannot efficiently accept sequencing adapters.
- 🟢 Preassembled Cas9/gRNA RNP makes fresh targeted cuts. These newly produced ends retain the 5′ phosphate needed for ligation.
- 🟢 Supplementary Table S1 prints each guide RNA's target-specific 20-base sequence and genomic target site. These guide oligos determine where the native DNA is cut; they are not incorporated into the sequenced library.
- 🟢 The same tube receives dATP and Taq polymerase (NEB M0267): 20 min at 37 °C for cleavage, then 5 min at 72 °C for A-tailing (article Online Methods).
- 🟢 **Published discrepancy:** the cited protocols.io V.2 (10 Sep 2019, step 9) gives 37 °C for **15 min** for Cas9 cleavage, versus 20 min in the published Online Methods. Both are the authors' documents and neither says which time each experiment used; the page shows both and the graph fixes only the reaction order.
- 🟢 RNPs: an equimolar crRNA pool (IDT) is duplexed with tracrRNA (IDT 1072532) at 95 °C 5 min, then 10 pmol duplex is combined with 10 pmol HiFi Cas9 Nuclease V3 (IDT 1081060) in CutSmart for 20 min at room temperature. Guide pairs cut complementary strands flanking the region.
- 🟢 Results: after cleavage Cas9 remains bound on the 5′ side of the gRNA, so adapters ligate preferentially onto the 3′ side of the cut; flanking guides therefore give coverage of both strands. The graph's discard reason for outward cut ends now states this instead of attributing their exclusion only to old dephosphorylated ends.
- 🟢 LSK109 sequencing adapters and Quick Ligase are added for 10 min at room temperature. A 0.3× AMPure XP cleanup with two Long Fragment Buffer washes and elution in 15 µL ONT Elution Buffer follows (protocols.io steps 13–20 first dilute the ligation with an equal volume of TE). The graph now draws this cleanup; sequencing buffer, loading beads and tether are then added, and samples run on R9.4.1 MinION or Flongle flow cells.
- 🟢 The native, PCR-free DNA supports simultaneous long-range sequence, structural-variant and CpG-methylation observation.
- 🟡 The page draws a target bracketed by two guide sites because that is the common whole-region design shown in Figure 1. The chemistry also supports one-sided reads and multiplexed guide sets.
- 🔴 LSK109 adapter bases and motor attachment are proprietary; roles are shown without invented sequence.

## Selectivity principle

The guide does not capture DNA by hybridization. Selectivity is encoded chemically: old breaks are 5′-OH after CIP, whereas Cas9 creates new 5′-phosphorylated, subsequently dA-tailed ends. Adapter ligation is therefore concentrated at chosen cut sites.
