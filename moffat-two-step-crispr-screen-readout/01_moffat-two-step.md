# Moffat two-step pooled CRISPR screen readout

The Moffat Lab TKOv3 workflow deliberately separates locus enrichment from indexed-library construction. PCR1 enriches the integrated guide region; PCR2 uses nested primers to add complete dual-index TruSeq arms.

## Source read

- 🟢 Moffat Lab, *Pooled CRISPR Screen Protocol v2*, revision 2019-07-23, section 3.5.3 and Tables 1–2, hosted by Addgene. The page follows the LCV2::TKOv3 branch.
- 🟢 The protocol specifies an approximately 600-bp PCR1 product, approximately 200-bp PCR2 product, and a run of 21 dark cycles + 26 Read 1 cycles + 8 cycles each for i7 and i5.

## Boundary

The same document mentions a possible one-PCR alternative. That alternative is not drawn here: the page is specifically the two-step workflow, while the Joung protocol has its own one-PCR page.
