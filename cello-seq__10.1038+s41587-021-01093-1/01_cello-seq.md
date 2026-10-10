# CELLO-seq — research notes

## Sources read

- 🟢 Berrens et al. 2022, DOI `10.1038/s41587-021-01093-1`, PMC article, Methods and Supplementary Table 5.
- 🟢 Berrens et al. 2025, DOI `10.1038/s41596-025-01203-2`, PMC detailed protocol and supplementary materials.

## Construct path

🟢 Poly(A)-primed, high-temperature RT uses a 3′-amino-blocked TSO. After Exonuclease I cleanup and denaturation, a pre-annealed splint is ligated to first-strand cDNA at 55 °C. It adds the cellular barcode and a 22-base UMI designed as repeated RYN classes to avoid long homopolymers. 🟢 High PCR duplicate counts are deliberate: noisy ONT reads sharing one long UMI are grouped to form a consensus suitable for locus-specific transposable-element assignment.
