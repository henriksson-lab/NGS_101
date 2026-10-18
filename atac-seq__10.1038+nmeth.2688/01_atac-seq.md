# ATAC-seq — research notes

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. The model is `tools/atac_seq.py`.

## 1. What it is

ATAC-seq profiles accessible chromatin by allowing adapter-loaded Tn5 to cleave and tag
native chromatin in one reaction.

## 2. Sources read

Buenrostro JD, Giresi PG, Zaba LC, Chang HY, Greenleaf WJ. "Transposition of native
chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding
proteins and nucleosome position." *Nature Methods* 10:1213–1218 (2013).
DOI [10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688).

| Source | Version | Sections read |
|---|---|---|
| `nmeth.2688_PMC3959825.xml.txt` | PMC author manuscript (NIHMS554473) of the published article | Online Methods lines 71–85 (nuclei, transposition, PCR, low-cell-number protocol, library QC); Primary data processing line 99 |
| `nmeth.2688_41592_2013_BFnmeth2688_MOESM200_ESM.pdf.txt` | publisher Supplementary Data and Methods (Springer) | Supplementary Table 1, lines 1–103 |

The publisher's typeset article is paywalled (https://www.nature.com/articles/nmeth.2688)
and was not read; the accepted author manuscript is used for the Methods. 🟡

## 3. Oligos

🟢 Supplementary Table 1 ("A list of ATAC-seq oligos used for PCR"), as printed, no
modifications given:

```
Ad1_noMX        AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG
Ad2.1_TAAGGCGA  CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT
...             Ad2.2–Ad2.24: same frame, 8-nt index window varies
```

🟡 (computed against `lib/`) Ad1_noMX = P5 + s5 + the first 7 nt of the mosaic end, with
no index; Ad2.N = P7 + 8-nt index (printed as the reverse complement of the name, e.g.
`TCGCCTTA` in Ad2.1_TAAGGCGA) + s7 + the first 6 nt of the mosaic end. The libraries are
therefore single-indexed (i7 only). The Methods call these "custom Nextera PCR primers 1
and 2". The transposase is the Illumina Nextera enzyme in 2× TD buffer; its loaded
adapter sequences are not printed in the paper 🔴 and are modelled from the canonical
Nextera s5/s7 + mosaic-end components. 🟡

## 4. Step by step

1. 🟢 **Nuclei**: 50,000 cells spun 500 × g 5 min, washed in 50 µL cold PBS, lysed in cold
   10 mM Tris-Cl pH 7.4, 10 mM NaCl, 3 mM MgCl₂, 0.1% IGEPAL CA-630, spun 500 × g 10 min.
2. 🟢 **Transposition**: 25 µL 2× TD buffer, 2.5 µL Transposase (Illumina), 22.5 µL water;
   37 °C 30 min. Tn5 fragments accessible DNA and transfers its adapters.
3. 🟢 **MinElute purification** (standard protocol). The 500- and 5,000-cell variant uses a
   5 µL transposition and takes it directly into the 50 µL PCR without MinElute.
4. 🟢 **Extension and PCR**: 1× NEBNext PCR master mix, 1.25 µM each primer; 72 °C 5 min,
   98 °C 30 s, then cycles of 98 °C 10 s, 63 °C 30 s, 72 °C 1 min. 🟡 The 72 °C step fills
   the transposition gaps so both tagged ends become primer-binding templates (the source
   gives the step, not the mechanism).
5. 🟢 **qPCR-limited cycling**: 5 cycles, then an aliquot with SYBR Green (0.6×) is run 20
   cycles to choose the additional cycles for the remaining 45 µL; 10–12 cycles in total.
6. 🟢 **Cleanup**: Qiagen PCR cleanup, ~30 nM in 20 µL. No size selection; quantified by
   qPCR (insert distribution 40 bp–1 kb, mean ~120 bp; fragments >2 kb seen).

## 5. Final library and reads

🟡 The productive two-ended molecule is the heterologous s5/s7 Nextera product:
`P5 · s5 · ME · insert · ME' · s7' · i7' · P7'`. Products with the same Tn5 end on both
sides lack one distinct PCR entry and are drawn as discarded off-products 🟡.

🟢 Reads: "34 × 8 × 34 reads from a MiSeq or 50 × 8 × 50 reads on a HiSeq" (two run
configurations, not a discrepancy). The page draws the HiSeq layout. Standard Nextera
Read 1 / Index 1 / Read 2 primers are implied by the adapter design, not named 🟡.

## 6. Open questions

- 🔴 Loaded Tn5 adapter sequences and the historical Nextera kit guide (not fetched).
- 🟡 Whether every Ad2.N index matches the Illumina Nextera i7 set — not checked.
- The publisher typeset article (paywalled) was not compared with the author manuscript.

## Source-owned validation migration, 2026-10-18

Original2013 declaration owns its exact Ad1_noMX/Ad2.1 example and single-eight-base-index Nextera endpoint. Removed batch_ngs and family primer/endpoint dependencies, including unused chromatin recipe import. Staggered pre-extension material remains distinct from the complete gap-filled template; PCR primer footprints are now computed on the actual extended precursor rather than final flow-cell product. End-class candidates lacking intended S5/S7 geometry remain unresolved instead of a physical discard. The50×8×50 HiSeq run is preserved; the34×8×34 MiSeq run remains another source variant.

Source scope: Buenrostro2013 published PMC3959825 Methods PCR line77 and Primary data processing line99; published MOESM200 SupplementaryTable1 beginning lines1–12 (Ad1_noMX and Ad2.1).

Deferred: Exact loaded transposome oligo chemistry and9nt stagger/gap-fill mechanism remain molecular-coverage work; Explicit qPCR aliquot branch and source34×8×34 MiSeq run variant remain separate modelling work. This bounded check does not establish complete molecular validation or an accepted source-absence claim.
