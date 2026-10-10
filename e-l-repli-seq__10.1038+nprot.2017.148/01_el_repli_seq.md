# E/L Repli-seq — research notes

## Citation

Marchal C et al. *Genome-wide analysis of replication timing by next-generation sequencing with E/L Repli-seq.* Nature Protocols 13, 819–839 (2018). [doi:10.1038/nprot.2017.148](https://doi.org/10.1038/nprot.2017.148)

## Sources read

- 🟢 The complete defining protocol in [PMC6044726](https://pmc.ncbi.nlm.nih.gov/articles/PMC6044726/).

## 1. What it is

🟢 E/L Repli-seq compares newly synthesized DNA from early- and late-S populations. It uses a BrdU pulse plus sorting and immunoprecipitation; BrdU is a selection mark, not a base conversion read directly by sequencing.

## 2. Construct-changing steps

1. 🟢 Pulse actively cycling cells for 2 h with 100 μM BrdU.
2. 🟢 Ethanol-fix and sort early- and late-S fractions by DNA content.
3. 🟢 Purify genomic DNA and sonicate it. The protocol processes the two fractions separately.
4. 🟢 End-repair the DNA, dA-tail it and ligate universal Illumina adapters before BrdU immunoprecipitation.
5. 🟢 Denature the adapter-ligated DNA and immunoprecipitate BrdU-containing nascent strands.
6. 🟢 Add indices by PCR, pool and sequence.

## 3. Final library and reads

🟢 These are conventional indexed Illumina genomic libraries. Normalized coverage is calculated separately for early and late fractions and reported as an early/late ratio.

🟡 The page uses canonical compatible Illumina geometry. Kit-specific universal-adapter oligos are not retyped from a commercial reagent.

## 4. Boundary

🟢 This page does not merge high-resolution multi-fraction Repli-seq or copy-number-only Sort-seq into E/L Repli-seq.
