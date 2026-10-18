# GoT — research notes

## Active native migration, 2026-10-18

Selected linear GoT from native 10x3′v2 recovered first-strand cDNA. Native version-specific kit amplification is bounded; its source extra cycle is explicit. SPRIselect recovery precedes physical allocation of3µl for genotyping and the distinct expression remainder. The genotyping aliquot receives source3–4 extra cycles with KAPA/10x primer mix and another cleanup, then SI-PCR plus a locus-specific reverse primer, SPRI, and P5generic/RPI-X PCR. The selected zero-stagger CALR reverse is assembled from the source handle row and gene-specific row; its two components are source-tested. The source P5generic is not placed on a pre-indexPCR product as if SI-PCR were a disclosed P5 primer. Cell/UMI identities retain lineage with native widths unresolved; the RPI sample index is source-positive8bases. Acquired linear3′v2 reads26/index8/98or130 are separate from physical native geometry.

Other 3′v3/5′ inputs and optional spike-ins are distinct experimental variants. Circular GoT remains outside the selected native graph, with its extensive source-known operations and discrepancy summarized publicly: optional SF3B1 pre-enrichment; PCR1;1.3XSPRI; USERII; overnight T4 circularization; lambda exonuclease;1.3XSPRI; PCR2 followed by the analogous conversion/closure/digestion/recovery; PCR3 linearization; SPRI; P5/RPI indexPCR; PE150MiSeq. The Gibson variant uses distinct complementary handles/KAPA and50°C1h assembly. Methods PCR2Fw/Rv U handles are AGGUCAGTCU/AGACUGACCU while Table5 assigns the reverse names; Methods PCR1Rv omits the Table5 terminal phosphorothioate “…CCGATC*T”. Neither ambiguity is silently repaired. Native circular topology/primer assignment requires a follow-up graph; the linear migration is not complete circular-variant coverage.

Digestion, nonamplification and computational filtering never imply physical discard. Actual first/native/targetPCR templates remain separate from their daughter endpoint, and unknown kit stages do not manufacture symmetric primer species.

This bounded pass read 2 originals and reran the mandatory bundle primer scraper. Original/source identity was checked in the external manifest. Historical review status and the two-reader missing-information gate are unchanged.

## Historical observations

The earlier observations below retain their original review scope; current source-owned declarations above replace earlier canonical endpoint and recipe defaults.


Defining source: Nam *et al.*, *Nature* (2019), [10.1038/s41586-019-1367-0](https://doi.org/10.1038/s41586-019-1367-0) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6782071/)).

## 1. What it is

Genotyping of Transcriptomes (GoT) adds targeted genotyping to a 10x single-cell transcriptome by making an amplicon that reaches from a chosen transcribed variant to the already attached 10x cell barcode and UMI.

## Sources read

- 🟢 The defining paper's “Genotyping of Transcriptomes (GoT)” and “Circularization GoT” Methods, run configurations and supplementary-primer descriptions were read.
- 🟢 2026-10-17 bounded pass: PMC6782071 JATS `xml.txt` line 196 (GoT Methods), lines 200–202 (Circularization GoT), line 216 (read screening); `sup_tab_5.xlsx.txt` lines 1–74. Supplementary Information Guide confirms Table 5 is the primer list.
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
| Reverse #1 | `AGACUGACCUCTACACGACGCTCTTCCGATCT` (Table 5 row 21 prints `…CCGATC*T`) |
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

### Published handle discrepancy (circularization GoT)

🟢 Methods (xml.txt line 200) give Fw-primer#2 `AGGUCAGTCU[3’end-locus-specific]` and Rv-primer#2 `AGACUGACCU[10bp-downstream-locus-specific]`; Supplementary Table 5 rows 24–25 list Handle PCR2_fw – Uhang `AGACUGACCU` and Handle PCR2_ rv – Uhang `AGGUCAGTCU`, the opposite assignment. Methods also print Rv-primer#1 without the 3′ phosphorothioate that Table 5 row 21 shows. The circularization branch is not drawn and no choice is modelled; the page caveat records both values.

## Construct-changing steps

0. 🟢 After SPRIselect cleanup of amplified cDNA, 3 µL (3′ v2) or 10 µL (3′ v3, 5′) is set aside for GoT. For 3′ v2 only, cDNA amplification gets one extra cycle and the GoT aliquot is re-amplified 3–4 cycles with KAPA HiFi HotStart ReadyMix and 10x primer mix (not drawn).
1. 🟢 Target PCR (linear GoT): the generic forward SI-PCR oligo (10x Genomics; sequence not printed) anneals partial Read 1 and preserves cell barcode and UMI; a locus-specific reverse primer carries a partial Read 2 handle, a stagger and the gene-specific region. Cycle number is set experimentally (10 for CALR). The “about 50 bases upstream” wording belongs to circularization PCR#1, not the linear primers.
2. 🟢 SPRI purification removes unincorporated primers; a second PCR with P5_generic (= `illumina.P5`, {{= len(illumina.P5) }} nt) and RPI-x completes P7 and adds the sample index. The library is spiked into the 10x library (HiSeq 2500) or sequenced separately (MiSeq).
3. 🟢 “Circularization GoT” is an explicitly named alternative amplification strategy in the same paper; it is noted but not silently substituted in the main linear schematic.

## Final library and read layout

- 🟢 Read 1 is 26 cycles (3′ v2 and 5′) or 28 cycles (3′ v3); Read 2 is 98 or 130 cycles; the sample-index read is 8 cycles. Circularization GoT was sequenced PE150 on MiSeq. The exact barcode/UMI split follows the chosen 10x chemistry.
- 🟢 The cell barcode and UMI are the molecular join key to the parent expression profile.
- 🟡 Exact 10x primer bases are not retyped from a secondary source; the final page records the barcode, UMI, target and read direction by role.
