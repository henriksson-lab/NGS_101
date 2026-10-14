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

🟡 Only the genomic portions of the inverse-PCR primers vary with the chosen viewpoint; their two tails and the universal second-round forward sequence are fixed in this protocol. The particular second-round reverse primer depends on the selected six-base index from Table 1.

## 4. Open questions

🔴 None affecting the generic topology. A concrete experiment must supply its primary/secondary enzyme pair, viewpoint-specific primer portions and selected index from Table 1.
