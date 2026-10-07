# SMART-seq and SMART-seq2 — quick reference

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not published.

**SMART-seq** — Ramsköld D, Luo S, Wang Y-C, *et al.* "Full-length mRNA-Seq from single-cell
levels of RNA and individual circulating tumor cells." *Nat Biotechnol* 2012;**30**(8):777–782.
doi:[10.1038/nbt.2282](https://doi.org/10.1038/nbt.2282)

**SMART-seq2** — Picelli S, Björklund ÅK, Faridani OR, Sagasser S, Winberg G, Sandberg R.
"Smart-seq2 for sensitive full-length transcriptome profiling in single cells."
*Nat Methods* 2013;**10**(11):1096–1098. doi:[10.1038/nmeth.2639](https://doi.org/10.1038/nmeth.2639)

Structure transcribed from
<https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html>.

---

## 1. What it is

Full-length single-cell mRNA sequencing. One cell per well, no barcoding, no UMI — plate-based
and low-throughput, but it reads the *whole* transcript, which is what makes it the method of
choice for isoform and allele-level work where 3'-counting methods cannot help.

🟢 SMART-seq2 is an optimisation of SMART-seq, not a redesign: the authors ran **457
optimisation experiments**. Two changes matter:

1. **An LNA at the TSO's 3' G.** Exchanging the last guanylate for a locked nucleic acid
   raises the stability of the three-base-pair TSO contact, making template switching
   markedly more efficient.
2. **Betaine plus higher MgCl₂.** Betaine is a methyl-group donor and destabilises secondary
   structure; combined with more Mg²⁺ it improves yield on structured templates.

The oligo designs are otherwise the same, and **the final libraries are identical**.

## 2. Oligos

```
oligo-dTVN   5'- AAGCAGTGGTATCAACGCAGAGT AC TTTTT...TTTTT VN -3'   (30 T)
TSO          5'- AAGCAGTGGTATCAACGCAGAGT ACAT rGrG+G -3'
ISPCR        5'- AAGCAGTGGTATCAACGCAGAGT -3'
```

`AAGCAGTGGTATCAACGCAGAGT` is the **SMART / ISPCR handle** (23 nt) — shared with Drop-seq,
SPLiT-seq and the Takara SMARTer kits. `+G` is the LNA; `rG` are ribonucleotides.

Back end is stock Nextera XT — see [`../ref/concepts/tn5-tagmentation.md`](../ref/concepts/tn5-tagmentation.md).

## 3. The key architectural fact

**The oligo-dT primer, the TSO and the PCR primer all carry the same handle.** After template
switching, both ends of the cDNA are identical, so a **single primer** (ISPCR) amplifies it.

This is elegant and it has one direct consequence: amplification is **semi-suppressive**.
A molecule with the same handle at both ends forms a pan-handle hairpin that competes with
primer annealing, which suppresses short products and keeps amplification length-biased
toward full-length cDNA. The same property later makes two of the three tagmentation products
non-amplifiable.

It also means **nothing distinguishes a 5' fragment from an internal one** — the symmetry that
makes amplification simple is exactly what SMART-seq3 has to break to get UMIs.

## 4. Steps

| # | Step | Note |
|---|---|---|
| 1 | Anneal oligo-dTVN to the poly(A) tail; reverse transcribe with **MMLV** | the `VN` anchor sets the start at the poly(A) junction |
| 2 | MMLV's terminal transferase adds untemplated **CCC** at the 5' end | see [`template-switching.md`](../ref/concepts/template-switching.md) |
| 3 | TSO anneals to the CCC and the polymerase switches template | the LNA G is the SMART-seq2 improvement |
| 4 | **ISPCR single-primer amplification** | semi-suppressive, as above |
| 5 | **Nextera tagmentation** of the amplified cDNA | 9-bp gap; three product classes, only s5/s7 amplifies |
| 6 | **72 °C gap fill-in** — the first cycle of the Nextera PCR | not an extension step; skipping it loses the library |
| 7 | PCR with N/S5xx and N7xx index primers | adds P5, i5, i7, P7 |

## 5. Final library

```
5'- P5 [8-bp i5] s5 ME --- cDNA --- revcomp(ME) revcomp(s7) [8-bp i7] revcomp(P7) -3'
```

🟢 Verbatim, and reproduced character-for-character by `tools/smartseq.py`
(asserted in `tools/selftest.py`):

```
AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
XXXXXXXX...XXXXXXXX
CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
```

**No UMI. No cell barcode** — the cell's identity is the well it came from, recorded by the
i5/i7 index pair.

## 6. Sequencing

Standard Nextera. Read 1 primes on `s5 + ME`; Read 2 on `s7 + ME`; the index reads use
`revcomp(ME) + revcomp(s7)` and `revcomp(ME) + revcomp(s5)`. Reads land anywhere in the
transcript, in either orientation — there is no positional information in the construct.
