# Moffat two-step pooled CRISPR screen readout

The Moffat Lab TKOv3 workflow deliberately separates locus enrichment from indexed-library construction. PCR1 enriches the integrated guide region; PCR2 uses nested primers to add complete dual-index TruSeq arms.

## Source read

- 🟢 Moffat Lab, *Pooled CRISPR Screen Protocol v2*, revision 2019-07-23, section 3.5.3 and Tables 1–2, hosted by Addgene. The page follows the LCV2::TKOv3 branch.
- 🟢 The protocol specifies an approximately 600-bp PCR1 product, approximately 200-bp PCR2 product, and a run of 21 dark cycles + 26 Read 1 cycles + 8 cycles each for i7 and i5.

🟢 [Table 1 of the Addgene-hosted protocol](https://media.addgene.org/cms/filer_public/44/9f/449f22c6-f818-4af6-8af4-2bab1d20f00a/moffat_lab_crispr_screen_protocol_v2.pdf)
prints the LCV2 PCR1 pair and the full i5/i7 PCR2 primer families. The schematic uses
one D501/D701 index combination as a representative library. Table 2 prints a distinct
pLCKO2 branch, whose different PCR1 reverse primer and i5 family must not be read as
interchangeable with this page's LCV2 construct.

## Boundary

The same document mentions a possible one-PCR alternative. That alternative is not drawn here: the page is specifically the two-step workflow, while the Joung protocol has its own one-PCR page.
