# RIC-seq — research notes

Defining source: Cai *et al.*, [10.1038/s41586-020-2249-1](https://doi.org/10.1038/s41586-020-2249-1).

## Evidence used

- 🟢 RIC-seq globally profiles RNA–RNA spatial interactions in situ.
- 🟢 Crosslinking preserves RNA–protein assemblies, fragmented RNA ends are labeled and spatially proximal RNAs are ligated.
- 🟢 Crosslinks are reversed and chimeric molecules are enriched and converted to a paired-end sequencing library.

## Protocol boundary

RIC-seq captures spatial proximity including protein-mediated contacts; it is not represented as a psoralen-selected base-paired duplex assay.

## Oligo audit

The [defining Nature article](https://doi.org/10.1038/s41586-020-2249-1) is the authority for this catalogue entry. Its fetched **Supplementary Table 9**, `MOESM12_ESM.xlsx.txt`, lines 248–262, directly lists the library PCR primers: 🟢 PE1.0 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` and twelve P7 index primers. Each P7 primer has `CAAGCAGAAGACGGCATACGAGAT` followed by a source-listed sample index and `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`. The twelve printed index windows are `AACGTGAT`, `ACAGGAGC`, `ACCGAGCT`, `ACCTTGAA`, `ACACATCG`, `ACGCCTAA`, `ACTGGTCA`, `CGAAAC`, `GATCTGG`, `TTGACTG`, `GGAACTG` and `TGACAT`; these have unequal lengths as printed, so they must not be flattened to a six-base placeholder without selecting a concrete primer.

The original authors' later [RIC-seq protocol for infected cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC11245969/) attributes the following adapter and PCR oligos to Cao *et al.*'s [2021 RIC-seq protocol](https://doi.org/10.1038/s41596-021-00524-2). These are 🟢 for the later protocol; the 2020 Supplementary Table 9 confirms the PCR arms but does not itself identify the two adapter oligos.

| Reagent in the later protocol | Sequence as printed, 5′→3′ |
|---|---|
| PEI adapter oligo A | `/5Phos/GATCGGAAGAGCACACGTCT` |
| PEI adapter oligo B | `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` |
| P5 primer | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` |
| P7 index primer | `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` |

Here `/5Phos/` is a 5′ phosphate and the six `N` bases are the index window. The later protocol also specifies an N6 random primer for first-strand synthesis, pCp-biotin for RNA-end labeling, dUTP in second-strand synthesis followed by USER digestion, and an A-tailed DNA intermediate before PEI adapter ligation. Its adapter and primer table is in the *Oligonucleotides* section; molecular steps are in *Adapter preparation* and *Strand-specific library construction*.

The current generic RIC-seq schematic omits these source-reported steps and exact PCR arms. Its shared factory needs a separate model revision after the original adapter's applicability is resolved.
