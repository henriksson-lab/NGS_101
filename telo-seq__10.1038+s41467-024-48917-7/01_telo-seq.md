# Telo-seq

## 1. What it is

Telo-seq selectively places an Oxford Nanopore motor adapter at the native human telomere end, so reads begin on the telomeric C-rich strand and proceed toward the subtelomere.

## 2. Sources read

- 🟢 Schmidt et al. (2024), [doi:10.1038/s41467-024-48917-7](https://doi.org/10.1038/s41467-024-48917-7), full text and Supplementary Table 7.
- 🟢 Oxford Nanopore, [TELO-seq requirements](https://nanoporetech.com/document/requirements/TELO-seq), current commercial implementation (additional, not defining).
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), comparative review.

## 3. Oligos and reaction order

🟢 Six 5′-phosphorylated telorettes share `AGCAATACGTAACTGAACGAAGT` and carry the six phases of an 18-nt telomere-complement arm. S1 is `ACTTCGTTCAGTTACGTATTGCTAGCAAT`.

1. 🟢 Telorettes anneal to the G-rich overhang and ligate to the C-rich chromosome strand.
2. 🟢 EcoRV cuts genomic DNA at blunt sites.
3. 🟢 Klenow exo− dA-tails those new internal blunt ends, reducing concatemer and wrong-end sequencing-adapter ligation.
4. 🟢 S1 anneals to the telorette and Oxford Nanopore AMII adapter; ligation completes the sequencing end.

## 4. Final library and sequencing

There is no sequencing primer. The telomere-side motor controls entry of the C-rich strand, producing a telomere-to-subtelomere read.
