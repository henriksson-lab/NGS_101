# SMART-seq — quick reference

> SMART-seq2's own oligo source review lives in `01_smart-seq2.md`, in the
> `smart-seq2__10.1038+nmeth.2639` directory. This note covers SMART-seq and the
> SMART-seq2 comparison it needs.

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not published.

**SMART-seq** — Ramsköld D, Luo S, Wang Y-C, *et al.* "Full-length mRNA-Seq from single-cell
levels of RNA and individual circulating tumor cells." *Nat Biotechnol* 2012;**30**(8):777–782.
doi:[10.1038/nbt.2282](https://doi.org/10.1038/nbt.2282)

**SMART-seq2** — Picelli S, Björklund ÅK, Faridani OR, Sagasser S, Winberg G, Sandberg R.
"Smart-seq2 for sensitive full-length transcriptome profiling in single cells."
*Nat Methods* 2013;**10**(11):1096–1098. doi:[10.1038/nmeth.2639](https://doi.org/10.1038/nmeth.2639)

Sources read (all primary; no secondary redrawing of this family is used as evidence):

- 🟢 Ramsköld et al. 2012, PMC3467340 full text and the NIHMS379463 supplement — prints the
  CDS oligo-dT primer verbatim, and *names* the SMARTer II A oligo and the kit's IS PCR
  primer without printing either.
- 🟢 Picelli et al. 2013 Springer supplementary files — Supplementary Table 2 ("List of all
  template switching oligonucleotides tested"), Supplementary Table 4 (which TSO
  *Smart-seq2* uses), and Supplementary Tables 1 and 3 (the oligo-dT and ISPCR primer by
  commercial name only).
- 🟢 Hagemann-Jensen et al. 2020 Supplementary Table 1, condition 1 ("Smartseq2 /
  Smart-seq2 as of Nature Protocol") — the Nature Protocols Smart-seq2 oligo-dT, TSO and
  forward PCR primer, printed in full.
- 🟢 Hahaut et al. 2022 Supplementary Note 2 — a second verbatim print of the same handle
  oligo-dT and of the ISPCR primer, attributed to Picelli et al. 2013.
- 🟡 Illumina's adapter-sequence documentation, via `lib/nextera.py` and `lib/illumina.py`,
  for ME, s5, s7, P5, P7 and the N5xx/N7xx index primers.

🔴 Not available: the Clontech/Takara SMARTer kit documentation. The SMARTer II A TSO and
the kit's IS PCR primer are commercial oligos that no source read here prints.

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

The indexed Nextera endpoint drawn here is the same architecture, but the original
SMART-seq paper does not disclose the full SMARTer II A template-switch oligo sequence.
Its exact equivalence to the later SMART-seq2 oligo therefore remains unverified. 🔴

## 2. Oligos

🟢 Ramsköld et al. print the original CDS primer verbatim as
`AAGCAGTGGTATCAACGCAGAGTACT(30)VN` — this is the SMART-seq oligo-dT, and the same
sequence reappears, printed in full, in Hagemann-Jensen et al. 2020 Supplementary
Table 1 and in Hahaut et al. 2022 Supplementary Note 2.

🟡 The TSO and ISPCR lines below are the later SMART-seq2 / Nature Protocols design, as
printed in Hagemann-Jensen et al. 2020 Supplementary Table 1 (TSO) and Hahaut et al. 2022
Supplementary Note 2 (ISPCR). Their equivalence to the commercial SMARTer II A reaction
that SMART-seq actually ran is **not** established by any source; the SMART-seq page
draws all three regions dotted-inferred, with the preamble caveat and an
`INFERRED — …` template-switch caption.

```
oligo-dTVN   5'- AAGCAGTGGTATCAACGCAGAGT AC TTTTT...TTTTT VN -3'   (30 T)   🟢
TSO          5'- AAGCAGTGGTATCAACGCAGAGT ACAT rGrG+G -3'                   🟡
ISPCR        5'- AAGCAGTGGTATCAACGCAGAGT -3'                               🟡
```

🟡 The closest thing to the kit TSO in any source read here is Picelli et al. 2013
Supplementary Table 2's all-DNA `dGCGGG` oligo,
`AAGCAGTGGTATCAACGCAGAGTACGCGGG`, which that table lists *alongside* "SMARTer Oligo IIA"
as a separate reagent — so it is a designed mimic, not a transcription of the kit oligo.

`AAGCAGTGGTATCAACGCAGAGT` is the **SMART / ISPCR handle** (23 nt) — shared with Drop-seq,
SPLiT-seq and the Takara SMARTer kits. `+G` is the LNA; `rG` are ribonucleotides.

Back end is stock Nextera XT — see [`../ref/concepts/tn5-tagmentation.md`](../ref/concepts/tn5-tagmentation.md).

🔴 The original paper names a SMARTer II A oligo and a commercial cDNA amplification
kit without printing those complete ordered oligos. It also reports both Nextera and
Covaris/adapter-ligation library preparation; this page follows the Nextera branch.

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

🟡 Derived, not quoted: `tools/smartseq.py` assembles this from the canonical P5/P7
(`lib/illumina.py`) and ME/s5/s7 (`lib/nextera.py`) constants, and `tools/selftest.py`
asserts the assembly rather than restating a published string. The handle region
contributed by the TSO does not survive tagmentation, so no inferred base reaches the
final library:

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
