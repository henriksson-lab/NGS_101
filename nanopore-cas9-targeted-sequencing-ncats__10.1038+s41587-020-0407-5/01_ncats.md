# nCATS — research notes

## Defining source

Gilpatrick T, Lee I, Graham JE, et al. *Targeted nanopore sequencing with Cas9-guided adapter ligation.* Nature Biotechnology 2020. [doi:10.1038/s41587-020-0407-5](https://doi.org/10.1038/s41587-020-0407-5), PMCID PMC7145730.

Sources read: primary-paper Results, Figure 1 legend and Online Methods (“Cas9 enrichment and library preparation”).

## Evidence record

- 🟢 About 3 µg HMW genomic DNA is dephosphorylated with Quick CIP for 10 min at 37 °C, then heated for 2 min at 80 °C. Existing DNA ends therefore cannot efficiently accept sequencing adapters.
- 🟢 Preassembled Cas9/gRNA RNP makes fresh targeted cuts. These newly produced ends retain the 5′ phosphate needed for ligation.
- 🟢 The same tube receives dATP and Taq polymerase: 20 min at 37 °C for cleavage, then 5 min at 72 °C for A-tailing.
- 🟢 LSK109 sequencing adapters and Quick Ligase are added for 10 min at room temperature. A 0.3× AMPure cleanup with long-fragment buffer follows.
- 🟢 The native, PCR-free DNA supports simultaneous long-range sequence, structural-variant and CpG-methylation observation.
- 🟡 The page draws a target bracketed by two guide sites because that is the common whole-region design shown in Figure 1. The chemistry also supports one-sided reads and multiplexed guide sets.
- 🔴 LSK109 adapter bases and motor attachment are proprietary; roles are shown without invented sequence.

## Selectivity principle

The guide does not capture DNA by hybridization. Selectivity is encoded chemically: old breaks are 5′-OH after CIP, whereas Cas9 creates new 5′-phosphorylated, subsequently dA-tailed ends. Adapter ligation is therefore concentrated at chosen cut sites.
