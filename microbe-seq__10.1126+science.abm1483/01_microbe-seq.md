# Microbe-seq — research notes

Defining source: Zheng et al., 2022, [doi:10.1126/science.abm1483](https://doi.org/10.1126/science.abm1483).

## Sources read

- 🟢 Linked bioRxiv preprint 10.1101/2020.12.14.422699, supplementary Methods (bead synthesis, barcoding and sequencing) and Supplementary Table S6 oligo workbook. The published *Science* supplementary material was not available through the source fetcher; preprint oligos are attributed to that version.

## Molecular path

🟢 One cell is lysed per droplet. Its genome is MDA-amplified before a tagmentation-reagent droplet is merged. A later barcode-bead droplet supplies a microbe-specific barcode and PCR primer.

🟢 The preprint's Table S6 actually prints **all 96 barcode-1 and 384 barcode-2 synthesis oligos**. The first barcode window has 24 members each of lengths 8, 9, 10 and 11 nt; barcode 2 is 8 nt. Thus the combined bead barcode is 16–19 nt, and the current 16-nt schematic token represents the shortest member only. The table's first examples are:

| Role | 5′→3′ oligo, spaces separate constant arms and barcode |
|---|---|
| Barcode-1 oligo-1 | `AAGGCGTCACAAGCAATCACTC AAACAAAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAG` |
| Barcode-2 oligo-1 | `TCTCTTATACACATCTGACGCTGCCGACGA AAACAAAC AAGGCGTCACAAGCAATCACTC` |

🟢 The supplementary Methods describe an acrydite-modified hydrogel oligo bearing a photocleavable linker, two split-and-pool extension rounds, the `P5PE1` and Nextera i7 library PCR primer roles, and a required custom Read-1 primer. The barcoding PCR takes place in droplets after MDA and Nextera tagmentation.

🔴 Table S6 does not print the initiating hydrogel oligo, `P5PE1`, or custom Read-1 primer sequence. Their exact bases and the linker chemistry beyond its named role remain open. The published *Science* supplement should be checked before the preprint oligo list is treated as final-version authoritative.
