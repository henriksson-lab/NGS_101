# GoT — research notes

Defining source: Nam *et al.*, *Nature* (2019), [10.1038/s41586-019-1367-0](https://doi.org/10.1038/s41586-019-1367-0) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6782071/)).

## 1. What it is

Genotyping of Transcriptomes (GoT) adds targeted genotyping to a 10x single-cell transcriptome by making an amplicon that reaches from a chosen transcribed variant to the already attached 10x cell barcode and UMI.

## Sources read

- 🟢 The defining paper's “Genotyping of Transcriptomes (GoT)” and “Circularization GoT” Methods, run configurations and supplementary-primer descriptions were read.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1038/s41586-019-1367-0`; the open manuscript is the readable primary source.

## Protocol variants kept explicit

- 🟢 Standard 10x Chromium 3′ v2, 3′ v3 or 5′ chemistry runs through emulsion break and first-strand-cDNA recovery.
- 🟢 For weakly detected 3′ targets, a gene-specific primer may be spiked into cDNA amplification. For 5′ libraries, a gene-specific RT primer and later additive primer improve recovery.
- 🟢 A small aliquot of amplified cDNA enters GoT; the remainder continues to the standard expression library.

## Oligo inventory

The defining article's [Supplementary Table 5](https://www.nature.com/articles/s41586-019-1367-0) is explicitly titled “List of primers used in GoT and single cell colony assay.” The fetched author-manuscript copy (`PMC6782071`, `sup_tab_5.xlsx.txt`, lines 1–74) supplies the library oligos below. `U` and phosphorothioate `*` are kept exactly as printed. The main text additionally gives the following **circularization GoT** primer architecture:

| Primer | Printed sequence, 5′→3′ |
|---|---|
| Forward #1 | `AGGUCAGTCU-[50bp-upstream-locus-specific]` |
| Reverse #1 | `AGACUGACCUCTACACGACGCTCTTCCGATCT` |
| Forward #2 | `AGGUCAGTCU[3’end-locus-specific]` |
| Reverse #2 | `AGACUGACCU[10bp-downstream-locus-specific]` |

🟢 These are the author manuscript's “Genotyping of Transcriptomes (GoT)” Methods and Extended Data Figure 1b–c, not sequences inferred from reads. The leading uracil-bearing tails form complementary cohesive ends after USER II digestion. They document the named circularization variant. The locus-specific portions vary with the interrogated mutation and must remain target-specific rather than a fictitious universal sequence.

### Linear GoT primers in Supplementary Table 5

| Oligo or role | Printed sequence, 5′→3′ |
|---|---|
| P5generic | `AATGATACGGCGACCACCGAGATCTACAC` |
| RPI-X indexing primer | `CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` |
| Staggered target-primer handles | `CACCCGAGAATTCCA`, `CACCCGAGAATTCCAA`, `CACCCGAGAATTCCATT`, `CACCCGAGAATTCCACAT`, each followed by `[Gene Specific Primer]` |
| SF3B1 targeted cDNA preamplification | `GGGCATTCCTTCTTTATTGC` |
| CALR L367fs*46 target portion | `GCAGCAGAGAAACAAATGAAGGAC` |
| CALR K385fs*47 target portion | `GAAGAAGACAAGAAACGCAAAGAGG` |
| SF3B1 H662D target portion | `TGTGCAAAAGCAAGAAGTCCT` |
| NFE2 E261fs*3 target portion | `TCAACTTGCCGGTAGATGAC` |
| XBP1 splice-site target portion | `GGAATGAAGTGAGGCCAGTG` |

🟢 Each named linear target primer is a staggered handle plus the printed target portion. `XXXXXXXX` is the variable sample index, not eight known bases. The table also prints the 5′-kit SpikeIn Guide RT handle `AGCAAGTGAGAAGCATCGTGTC-[Gene Specific Primer]` and additive primer `AGCAAGTGAGAAGCATCGTG*T*C`; the two terminal `*` marks are phosphorothioate linkages. It gives locus-specific circularization primers and alternative Gibson handles as separate rows. Thus the GoT-specific oligos are largely recovered, while the upstream 10x kit oligos and a concrete choice of index/target remain implementation-dependent.

## Construct-changing steps

1. 🟢 Target PCR starts about 50 bases upstream of the variant and amplifies toward the 10x end, thereby preserving cell barcode and UMI.
2. 🟢 Further amplification/indexing makes the genotype library. It may be sequenced alone or spiked into the expression library.
3. 🟢 “Circularization GoT” is an explicitly named alternative amplification strategy in the same paper; it is noted but not silently substituted in the main linear schematic.

## Final library and read layout

- 🟢 For 3′ v3: Read 1 is 28 cycles, Read 2 is 98 or 130 cycles, and the sample-index read is 8 cycles. The exact barcode/UMI split follows the chosen 10x chemistry.
- 🟢 The cell barcode and UMI are the molecular join key to the parent expression profile.
- 🟡 Exact 10x primer bases are not retyped from a secondary source; the final page records the barcode, UMI, target and read direction by role.
