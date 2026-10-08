# TeloPCR-seq

## 1. What it is

TeloPCR-seq attaches an anchor to native *S. pombe* chromosome termini, amplifies each end with a subtelomeric primer, and reads the resulting telomere amplicons as PacBio circular libraries.

## 2. Sources read

- 🟢 Bennett et al. (2016), [doi:10.1002/1873-3468.12444](https://doi.org/10.1002/1873-3468.12444), abstract and methods metadata.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), method overview.

## 3. Construct and reaction order

1. 🟢 CircLigase attaches a single-stranded anchor at the native telomere terminus.
2. 🟢 A chromosome-end-specific subtelomeric primer and the anchor primer amplify the complete telomere tract.
3. 🟢 The amplicon is converted into a circular PacBio template.

The exact anchor and primer order-sheet bases are not in the locally accessible article text. They remain `N`/named regions in the model rather than inferred sequence.

## 4. Final library and sequencing

The selected telomere amplicon lies between two hairpin adapters. A PacBio sequencing primer binds a hairpin primer site and the polymerase makes repeated passes around the closed template.
