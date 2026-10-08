# NEBNext Low-bias Small RNA

## 1. What it is

A gel-free small-RNA library protocol that uses randomized splint ligation to reduce sequence-dependent adaptor ligation bias.

## 2. Sources and evidence

- 🟢 NEB E3420 manual v2.0, September 2025, and [product page](https://www.neb.com/en-us/products/e3420-nebnext-low-bias-small-rna-library-prep-kit).
- 🟢 Captured RNAs are under 120 nt and require a 5′ monophosphate and 3′ hydroxyl.
- 🟢 The order is 3′-adaptor ligation, excess-adaptor removal, simultaneous 5′ ligation plus modification of the 3′ adaptor into an RT substrate, RT, bead selection and UDI PCR.
- 🔴 NEB does not print the adaptor or randomized-splint sequences; these remain visibly unknown.

## 3. Molecular path

The protocol directly ligates RNA ends. It therefore does not inherit the dA-tailed DNA junction used by Ultra II DNA workflows. A typical 56-base Read 1 is sufficient for the small-RNA insert.
