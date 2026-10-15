# SMART-seq3xpress — oligo source review

Hagemann-Jensen M, Ziegenhain C, Sandberg R. "Scalable single-cell RNA sequencing from
full transcripts with Smart-seq3xpress." *Nat Biotechnol* 2022;**40**(10):1452–1457.
doi:[10.1038/s41587-022-01311-4](https://doi.org/10.1038/s41587-022-01311-4)

> **Evidence marking.** 🟢 verbatim from a source read · 🟡 derived or inferred ·
> 🔴 not published in any source read here.

## Sources read

- 🟢 **The defining paper's own Methods**, open access under CC BY as PMC9546772. It
  prints the two TSOs side by side: "Original Smartseq3 TSO
  (5′-Biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′; IDT). Improved TSO
  (5′-Biotin-AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG-3′; IDT)", and both PCR primers,
  "Smart-seq3 forward primer (5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG-3′)" and
  "Smart-seq3 reverse primer (5′-ACGAGCATCAGCAGCATACGA-3′)". It also gives the UMI-read
  pattern, `ATTGCGCAATG` with up to two mismatches.
- 🟢 **The authors' protocols.io protocol** `smart-seq3xpress-bwh4pb8w` — the ordered oligo
  table, with supplier, purification and concentration:
  `/5BiosG/ACGAGCATCAGCAGCATACGA(T)30VN`, `/5BiosG/AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG`,
  `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G` and `ACGAGCATCAGCAGCATAC*G*A`, with a
  footnote that `*` marks phosphorothioate bonds.
- 🟢 **Springer supplementary files** `41587_2022_1311_MOESM1/2_ESM.pdf` — supplementary
  figures; no additional oligo table.
- 🟢 **Hagemann-Jensen et al. 2020 Supplementary Table 1** for the shared SMART-seq3 core,
  recorded in that protocol's note.
- 🟡 **Illumina adapter-sequence documentation**, via `lib/nextera.py` and
  `lib/illumina.py`, for ME, s5, s7, P5, P7 and the index-primer architecture.

Every modelled base is a transcription of one of the above. No secondary redrawing of the
SMART-seq family is used as evidence here, and none is needed: the paper and the author
protocol each print the complete oligo set independently.

## Oligos

```
Smartseq3_OligodT30VN  5'- /5BiosG/ ACGAGCATCAGCAGCATACGA TTTTT...TTTTT VN -3'   (30 T)  🟢
Smartseq3xpress_TSO    5'- /5BiosG/ AGAGACAGATTGCGCAATG [8-bp UMI] WW rGrGrG -3'         🟢
Fwd_PCR_primer         5'- TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G -3'            🟢
Rev_PCR_primer         5'- ACGAGCATCAGCAGCATAC*G*A -3'                                   🟢
```

🟢 The only oligo difference from SMART-seq3 is the **`WW` spacer** between the UMI and the
G tail — the "Improved TSO" above. `AGAGACAG` is the 3'-terminal 8 nt of the Nextera
mosaic end, so the stock Read 1 primer (`s5 + ME`) ends exactly where the 11-bp tag
`ATTGCGCAATG` begins; the forward PCR primer rebuilds the whole s5 arm as it amplifies.
`W` is A or T, so the spacer is a two-base degenerate window, not a fixed sequence.

🟡 Everything else that differs from SMART-seq3 is reaction engineering rather than
chemistry: sub-microlitre volumes under an overlay, SeqAmp or NEBNext polymerase options,
and post-pre-amplification clean-up variants. None of it changes a modelled base.

## What the page asserts, by role

| Role | Status | Where it comes from |
|---|---|---|
| biotinylated oligo-dT with its own handle | 🟢 | protocols.io ordered table; shared with SMART-seq3 |
| UMI-bearing TSO with `WW` spacer | 🟢 | defining paper's Methods, "Improved TSO"; protocols.io |
| 8-base UMI | 🟢 | `NNNNNNNN` in both sources |
| 11-bp 5' tag `ATTGCGCAATG` | 🟢 | both sources, and the read-parsing pattern |
| forward / reverse PCR primers, incl. terminal phosphorothioates | 🟢 | protocols.io ordered table |
| Nextera ME / s5 / s7, read and index primer entry | 🟢 | Illumina adapter sequences, via `lib/nextera.py` |
| P5 / P7 flowcell ends | 🟢 | Illumina adapter sequences, via `lib/illumina.py` |
| i5 / i7 index *members* | 🔴 | never asserted; see below |

## Open roles

- 🔴 **Index members and length.** The author protocol orders *custom* Nextera index
  primers and recommends a 10-bp barcode designed with the `DNABarcodes` R package
  (minimum Levenshtein distance 3, no homopolymers ≥ 3, no i5 starting `AC` on NovaSeq
  v1), noting that 8 bp — Illumina's length — works as well. The page draws an 8-base
  index window and names no member; the 8-vs-10 choice is a protocol option, not a fact
  about this chemistry.
- 🟡 The drawn tagmentation and gap-fill geometry is the shared Nextera back end from
  `lib/nextera.py`, not an xpress-specific transcription.

## Everything else

The library architecture, the two fragment classes (5' fragments carrying tag + UMI versus
internal fragments) and the read layout are identical to SMART-seq3's and are written up
once, in `01_smart-seq3.md` in the `smart-seq3__10.1038+s41587-020-0497-0` directory.
