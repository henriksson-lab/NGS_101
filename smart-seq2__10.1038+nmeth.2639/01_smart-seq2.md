# SMART-seq2 — oligo source review

Picelli S, Björklund ÅK, Faridani OR, Sagasser S, Winberg G, Sandberg R. "Smart-seq2 for
sensitive full-length transcriptome profiling in single cells." *Nat Methods*
2013;**10**(11):1096–1098. doi:[10.1038/nmeth.2639](https://doi.org/10.1038/nmeth.2639)

> **Evidence marking.** 🟢 verbatim from a source read · 🟡 derived or inferred ·
> 🔴 not published in any source read here.

This note exists to settle one question: where does every base on the SMART-seq2 schematic
come from? The short answer is that all of them are transcriptions of primary documents,
with one two-base exception that the page now draws as inferred.

## Sources read

- 🟢 **Picelli et al. 2013 Supplementary Table 2** (Springer
  `41592_2013_BFnmeth2639_MOESM301_ESM.xlsx`), "List of all template switching
  oligonucleotides tested" — thirteen TSOs with their 5' and 3' modifications.
- 🟢 **Picelli et al. 2013 Supplementary Table 4** (`…MOESM303…`), "Detailing variants of
  the Smart-seq2 protocol" — the row named `Smart-seq2` uses TSO **`rGrG+G`**, 1 µl of
  10 µM, no bead purification before pre-amplification, KAPA HiFi HotStart ReadyMix.
- 🟢 **Picelli et al. 2013 Supplementary Tables 1 and 3** (`…MOESM300…`, `…MOESM302…`) —
  per-library reaction tables. They identify the oligo-dT only by its commercial name,
  *SMARTer dT30VN*, and the PCR primer only as the kit's *ISPCR primers*.
- 🟢 **Hagemann-Jensen et al. 2020 Supplementary Table 1**, condition 1 ("Smartseq2 —
  Smart-seq2 as of Nature Protocol") — prints that version's oligo-dT
  `AAGCAGTGGTATCAACGCAGAGTAC(T)30VN`, TSO `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G` and
  forward PCR oligo `AAGCAGTGGTATCAACGCAGAGT` in full.
- 🟢 **Hahaut et al. 2022 Supplementary Note 2** — its *Smart-seq2* block prints
  `Smart dT30VN` = `/5Biosg/AAGCAGTGGTATCAACGCAGAGTAC(T)30VN` and `ISPCR` =
  `/5Biosg/AAGCAGTGGTATCAACGCAGAGT`, both credited to Picelli et al. 2013.
- 🟢 **Ramsköld et al. 2012** (PMC3467340) — the identical SMART-seq CDS primer
  `AAGCAGTGGTATCAACGCAGAGTACT(30)VN`, which is where this oligo-dT comes from.
- 🟡 **Illumina adapter-sequence documentation**, via `lib/nextera.py` and
  `lib/illumina.py`, for ME, s5, s7, P5, P7 and the Nextera (XT) N5xx/N7xx index primers.

🔴 **Still out of reach.** The defining paper has no open full text: `get_sources.py`
records it as `(manual) nmeth.2639 full text … no open full text found`. Its *Nature
Protocols* companion, doi:10.1038/nprot.2014.006, is the primary ordered-oligo table for
the version everyone runs and is now requested in `catalogue/extra_sources.tsv`. The
Clontech/Takara SMARTer kit documentation, which would be the primary source for the
*SMARTer dT30VN* and *ISPCR* reagents under those names, is also absent.

This review supersedes an earlier provenance note that credited part of the exact oligo
bases to a secondary drawing of the SMART-seq family. That is not where they come from,
and that drawing is not evidence for anything on these pages.

## Oligos

```
oligo-dTVN   5'- AAGCAGTGGTATCAACGCAGAGT AC TTTTT...TTTTT VN -3'   (30 T)   🟢
TSO          5'- AAGCAGTGGTATCAACGCAGAGT ACAT rGrG+G -3'                   🟢 / 🟡 (see below)
ISPCR        5'- AAGCAGTGGTATCAACGCAGAGT -3'                               🟢
```

`AAGCAGTGGTATCAACGCAGAGT` is the 23-nt **SMART / ISPCR handle**, shared by the oligo-dT,
the TSO and the PCR primer — which is the whole architectural point of the method. `rG`
is a ribonucleotide and `+G` the locked nucleic acid that is SMART-seq2's headline
improvement over SMART-seq.

### The one disagreement between sources: AC or ACAT?

🟡 Two published forms of the TSO linker exist, and they differ by two bases:

| Form | Printed in | TSO |
|---|---|---|
| `AC` | Picelli et al. 2013 Supplementary Table 2, row `rGrG+G` — the row Supplementary Table 4 names as *Smart-seq2*'s TSO | `AAGCAGTGGTATCAACGCAGAGT` **`AC`** `rGrG+G` |
| `ACAT` | Hagemann-Jensen et al. 2020 Supplementary Table 1, condition 1, labelled "Smart-seq2 as of Nature Protocol" | `AAGCAGTGGTATCAACGCAGAGT` **`ACAT`** `rGrG+G` |

Supplementary Table 2 is internally consistent about this: three of its thirteen oligos
(`C6 Amino`, `rGrGrG`, `rG3p`) carry `ACAT`, and the other ten, `rGrG+G` among them, carry
`AC`. So the defining paper's own chosen oligo has the shorter linker, and the longer one
belongs to the *Nature Protocols* version of the protocol published a year later.

The model keeps the longer, `ACAT` form, because that is what the protocol as practised
uses and what `lib/smartseq2.TSO_LINKER` supplies to the SMART-seq2 derivatives
(scM&T-seq, scMT-seq, scNMT-seq, scCAT-seq). The page therefore draws **`AC` as
established and the trailing `AT` as inferred**, with all three required signals: a
preamble caveat naming both forms, an `INFERRED — …` caption on the template-switch step,
and `<inf>` around the two bases in the oligo listing and the step panels. No base is
presented as published by the defining paper when it is not.

🟡 This two-base question has no consequence for the final library: the handle region the
TSO contributes is internal to the amplified cDNA and is tagmented away, so it does not
appear in the sequenced construct at all. It matters only for ordering the oligo.

## What the page asserts, by role

| Role | Status | Where it comes from |
|---|---|---|
| oligo-dT (dT30VN with SMART handle) | 🟢 | Ramsköld 2012 CDS primer; Hagemann-Jensen 2020 Suppl. Table 1; Hahaut 2022 Suppl. Note 2 |
| TSO handle + `AC` + `rGrG+G` | 🟢 | Picelli 2013 Suppl. Tables 2 and 4 |
| TSO linker's trailing `AT` | 🟡 | Hagemann-Jensen 2020 Suppl. Table 1 (Nature Protocols form); drawn inferred |
| ISPCR primer | 🟢 | Hagemann-Jensen 2020 Suppl. Table 1; Hahaut 2022 Suppl. Note 2 |
| Nextera ME / s5 / s7, read and index primer entry | 🟢 | Illumina adapter sequences, via `lib/nextera.py` |
| P5 / P7 flowcell ends | 🟢 | Illumina adapter sequences, via `lib/illumina.py` |
| i5 / i7 index *members* | 🔴 | never asserted: drawn as 8-base index windows, as Nextera XT specifies |

## Open roles

- 🔴 The chosen i5/i7 index pair per well is not established and is not drawn.
- 🟡 The TSO linker's trailing `AT` rests on the *Nature Protocols* companion rather than
  the defining paper; fetching doi:10.1038/nprot.2014.006 would turn it 🟢 for that
  version of the protocol, but not for the 2013 paper.
- 🔴 The commercial identity of *SMARTer dT30VN* and of the kit *ISPCR primers* under
  those names is not confirmable from any source read here — only that oligos with these
  sequences are what the later primary sources print for this method.

## Chemistry and steps

The chemistry, the step order and the final library are the same as SMART-seq's, with the
LNA TSO, betaine and higher MgCl₂ as the optimisations. They are written up once, in
`01_smart-seq.md` in the `smart-seq__10.1038+nbt.2282` directory, which this note does not
repeat.
