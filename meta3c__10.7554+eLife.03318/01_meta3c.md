# meta3C — research notes

## Sources read

- 🟢 Marbouty et al. 2014, DOI `10.7554/eLife.03318`, eLife PDF/XML and Materials and methods.

## Construct path

🟢 Mixed microbial cells are crosslinked before lysis. The bacterial demonstration uses HpaII (`C^CGG`); the yeast mixture uses DpnII (`G^ATC`). Crosslinked restriction fragments are ligated under dilute conditions, de-crosslinked and sheared. 🟢 Unlike Hi-C, meta3C deliberately omits junction enrichment, leaving enough ordinary paired-end shotgun sequence for assembly while contact chimeras provide species grouping and scaffolding information.

## Sequencing-library oligos

🟢 Table 1 of the defining paper prints both custom Illumina PE adapter oligos for each organism mixture. The yeast and river-sediment pair is MM70 `GTANNNNNNAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG` and MM71 `ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNTACT`. The bacterial pair is MM76 `CAGNNNNNNAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG` and MM77 `ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNCTGT`. The six random bases on each oligo are used to collapse PCR duplicates; the fixed terminal tags separate the libraries after sequencing. These custom adapter bases were omitted from the former role-only schematic.

🔴 The paper says that the library was processed using Illumina Paired-End DNA Sample Prep Kit PE-930-1001 and amplified by test and large-scale PCR, but it does not print the kit PCR primer oligos. The exact ligation-end geometry of the custom adapter solution is also not specified by Table 1, so the schematic retains that limitation.
