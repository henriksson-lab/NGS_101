# NEBNext Ultra II FS DNA PCR-free

## 1. What it is

A PCR-free Illumina DNA workflow combining enzymatic fragmentation/end preparation with full-length, dual-index UMI adaptors.

## 2. Sources and evidence

- 🟢 NEB E7430/E7435 manual v3.0, September 2025, and [product page](https://www.neb.com/en-us/products/e7435-nebnext-ultra-ii-fs-dna-pcr-free-library-prep-with-sample-purification-beads).
- 🟢 For high input, 10 minutes at 37°C targets 350-bp inserts and 8 minutes targets 450-bp inserts; both receive 30 minutes at 65°C.
- 🟢 The current UDI-UMI adaptor contains 8-base i5 and i7 indices plus a 12-base UMI; NEB specifies a 20-cycle Index 1 read comprising i7 then UMI.
- 🔴 NEB withdrew the adaptor sequence as inaccurate and now treats the correct sequence as proprietary. No historical sequence is used here.

## 3. Molecular path

The full-length adaptor already supplies every sequence needed for cluster generation and sequencing. Unlike the USER-hairpin Ultra II page, there is no loop opening and no PCR completion step.
