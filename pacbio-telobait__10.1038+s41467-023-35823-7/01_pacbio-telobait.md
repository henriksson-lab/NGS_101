# PacBio telobait HiFi sequencing

## 1. What it is

PacBio telobait sequencing ligates phased, barcoded and biotinylated baits to native human telomere overhangs, affinity-purifies intact chromosome ends, then makes SMRTbell HiFi libraries.

## 2. Sources read

- 🟢 Tham et al. (2023), [doi:10.1038/s41467-023-35823-7](https://doi.org/10.1038/s41467-023-35823-7), full text and Supplementary Data 1–2.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), comparison of targeted long-read methods.

## 3. Telobait oligos

🟢 One published forward bait is `/5Phos/ATCGACGGTTCAAGATGCCAGATGCACGGAGCAGAATTCTTTTAC/3Bio/`. Its six reverse partners terminate in `CCCTAA`, `ACCCTA`, `AACCCT`, `TAACCC`, `CTAACC`, or `CCTAAC`, covering all phases of the human repeat. The bait carries a barcode, an EcoRI release site, and 3′ biotin.

## 4. Reaction order

1. 🟢 RsaI and HinfI trim the genome without cutting canonical TTAGGG repeats.
2. 🟢 A six-phase telobait pool anneals and is ligated at native telomere overhangs; T4 DNA polymerase fills the gap.
3. 🟢 Streptavidin capture and washing select tagged chromosome ends.
4. 🟢 EcoRI releases the selected fragments.
5. 🟢 PacBio SMRTbell preparation converts them to closed HiFi templates.

## 5. Final library and sequencing

A PacBio primer/polymerase complex enters at the hairpin primer site and repeatedly traverses the selected telomere fragment and its reverse complement.
