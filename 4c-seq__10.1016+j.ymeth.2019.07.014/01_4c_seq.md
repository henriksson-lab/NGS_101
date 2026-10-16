# 4C-seq — research notes

## Citation

Krijger PHL, Geeven G, Bianchi V, Hilvering CRE, de Laat W. *4C-seq from beginning to end: A detailed protocol for sample preparation and data analysis.* Methods 170, 17–32 (2020). [doi:10.1016/j.ymeth.2019.07.014](https://doi.org/10.1016/j.ymeth.2019.07.014)

## Sources read

- 🟢 The defining open-access protocol, including its abstract, Fig. 1 workflow and library-generation sections.
- 🟢 The authors' institutional full-text mirror linked in `ref/MANIFEST.md`, especially §3.4 and Table 1 for the published first-round tails and universal second-round primers.

## 1. What it is

🟢 4C-seq asks which genomic fragments contact one chosen restriction fragment (the viewpoint). Formaldehyde fixation, primary restriction digestion and proximity ligation first make a conventional 3C library.

🟢 The purified 3C DNA is cut with a different, secondary restriction enzyme. A second dilute ligation circularizes these shorter products. Outward-facing primers in the known viewpoint then inverse-PCR the unknown ligated reporter.

## 2. Construct-changing steps

1. 🟢 Crosslink cells, digest chromatin with a primary cutter and proximity-ligate spatially adjacent restriction fragments.
2. 🟢 Reverse crosslinks and purify the 3C library.
3. 🟢 Digest again with a secondary cutter; this reduces fragment size and defines the second boundary read from each reporter.
4. 🟢 Ligate under dilute conditions to make circular 4C templates.
5. 🟢 Perform viewpoint-specific inverse PCR. The protocol uses two PCR stages; the second completes Illumina-compatible adapters and sample multiplexing.

## 3. Final library and reads

🟢 The genomic read begins from a viewpoint-specific reading primer and proceeds across the primary restriction-ligation junction into the captured reporter. The opposite end is constrained by the secondary cutter.

🟢 The **reading** inverse-PCR primer has the published 5′ tail `TACACGACGCTCTTCCGATCT`, then a 0–5-nt optional spacer and the viewpoint-specific sequence. The **non-reading** primer has the published 5′ tail `ACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, then its viewpoint-specific sequence. Table 1 gives the second-round universal forward primer, 5′-`AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`-3′, and an indexed reverse-primer panel.

🟢 Table 1 contains 24 named reverse-primer rows. The later rows, numbered 13–27, include two additional bases before the six indexing bases; the caption explicitly says these two bases are not used for indexing. The diagram represents the simpler rows 1–12 family without selecting an actual sample index.

🟡 The diagram uses the published DpnII/Csp6I combination as an illustrative branch: Table 3 uses it for the Sox2 and Mccc1 examples. Enzyme choice and genomic primer arms are experiment-specific. The optional spacer, locus sequence and sample index remain symbolic; their unknown lengths are not replaced by arbitrary biological sequence lengths.

🟢 Section 5.2 uses two PCR stages: 16 cycles, 0.8× AMPure cleanup, then 20 cycles and QIAquick purification. Section 6 describes single-end 50-base HiSeq or 75-base MiniSeq/NextSeq reads with a six-base index. The model does not assign physically loaded sequencing primers from PCR-arm identity alone.

🟢 The paper describes amplified self-ligation and undigested-fragment products (note 12, extracted lines 1226–1234). These are drawn as separate possible amplicon classes. Their exclusion from analytical metrics does not demonstrate physical removal from the library.

## 4. Open questions

🟡 A concrete experiment must supply its enzyme pair, viewpoint-specific primer portions, optional spacer and selected index. These are unassigned model choices, not confirmed omissions from the publication.

The source audit remains in progress. Appendix A identifies publisher-hosted supplementary data (the article also mentions supplementary Fig. 1); its actual attachment has not yet been inspected. The linked pipe4C and peakC repositories also remain pending original-source scope checks. The retrieved Illumina pooling chapter is the current v17 document, a later comparator; the paper does not identify a historical revision to assign. No missing-information claim is accepted from this incomplete linked-source inventory.
