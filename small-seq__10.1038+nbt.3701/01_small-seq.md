# Small-seq — single-cell small-RNA sequencing

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published · ✅ computed and asserted in `tools/selftest.py`.

**Small-seq** — Hagemann-Jensen M, Abdullayev I, Sandberg R, Faridani OR. "Small-seq for
single-cell small-RNA sequencing." *Nat Protoc* 2018;**13**(10):2407–2424.
doi:[10.1038/s41596-018-0049-y](https://doi.org/10.1038/s41596-018-0049-y)

Original method: Faridani OR *et al.* "Single-cell sequencing of the small-RNA
transcriptome." *Nat Biotechnol* 2016;**34**(12):1264–1266.
Pipeline: <https://github.com/eyay/smallseq>.

Source text used here: `pdf/hagemann.txt` (the PDF lives in `pdf/`, which is gitignored).
**Supplementary Table 1 is not in that text**, which is why one oligo below is 🔴.

---

## 1. Why this protocol is in the repo

Every other protocol here is a TruSeq-DNA or Nextera library with a different front end.
Small-seq is not. It brings four pieces of vocabulary the repo did not have:

| | New thing | Where else it turns up |
|---|---|---|
| 1 | **TruSeq *Small RNA* adapters** (RA5 / RA3) instead of the TruSeq DNA adapter | every small-RNA kit; the read-1 landing site is a different sequence entirely |
| 2 | **Sequential ligation** — 3' adapter, then destroy the leftovers, then 5' adapter | NEBNext Small RNA, QIAseq miRNA, CLIP/iCLIP variants |
| 3 | **A UMI in the 5' adapter**, ligated on, and the *first* thing sequenced | contrast SMART-seq3 (UMI in a TSO) and CRISPR-MIP (UMI on a probe) |
| 4 | **A masking oligonucleotide** — rRNA depletion by occupying an end, not by pulldown | any ligation protocol with one dominant unwanted species |

No polymerase is involved in attaching either adapter. There is no tagmentation, no
template switching, no oligo-dT, and no rGrGrG tail — all asserted ✅.

## 2. Oligos

🟢 All five verbatim from *Reagent setup* (p. 2414), which states: *"All oligonucleotides
are listed in the 5′ to 3′ direction."* Line-wrap spaces in the extracted text are closed
up; nothing else is changed.

```
RA5   5'- NH2-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC rHrHrHrHrHrHrHrH rCrA -3'
RA3   5'- rApp TGGAATTCTCGGGTGCCAAGG ddC -3'
RTP   5'- biotin CCTTGGCACCCGAGAATTCC rA -3'
RP1   5'- AATGATACGGCGACCACCGAGATCTACAC GTTCAGAGTTCTACAGTCCGA -3'
SRX   🔴 not in the main text — "see Supplementary Table 1"

mask  5'- ATCGGCAAGCGACGCTCAGACAGGCGTAGCCCCGGGAGGAACCCGGGGCCGCAAGTGCGTTCGAAGTGTCGATGAT biotin -3'
```

Lengths ✅: RA5 26 nt of adapter + 8 UMI + 2 linker = **36 nt**; RA3 **21 nt**;
RTP **21 nt**; RP1 **50 nt**; mask **76 nt**.

### The four of them interlock exactly

✅ All four relationships are computed, not quoted:

- **`RTP == revcomp(RA3)`**, base for base. That is not a coincidence, it is the
  mechanism of step 3: RTP anneals to leftover RA3 to make a blunt 21-bp duplex, which is
  what lambda exonuclease wants.
- **`RP1 == P5 + RA5[:21]`**, where `P5` is the canonical flow-cell sequence from
  `lib/illumina.py` (a vendor document, not this paper). So RP1 is the small-RNA P5 arm:
  29 nt of P5 welded to the first 21 nt of the 5' adapter.
- **RP1 stops 5 nt short of RA5's 3' end.** The missing five are `CGATC`. This matters —
  see [`02_sequencing_primer.md`](02_sequencing_primer.md).
- **RA3 begins `TGG`**, which is why the pipeline's precursor filter is called
  `remove_reads_with_genomic_TGG.py` 🟡 — reads from longer RNAs whose *genomic* sequence
  happens to continue `TGG…` get falsely trimmed into "small RNAs".

### Modifications, and what each is for

| Oligo | Modification | Purpose |
|---|---|---|
| RA5 | 5'-NH₂ (aminolink C6) | 🟢 blocks the 5' end, so RA5 cannot be a ligation *acceptor* — no 5'-adapter concatemers |
| RA5 | 3'-OH, all-RNA | 🟢 T4 RNA ligase 1 substrate |
| RA3 | 5'-rApp (**chemically** pre-adenylated) | 🟢 carries its own activation, so no ATP is needed; 🟢 "Enzymatic adenylation of the oligo after synthesis may reduce the efficiency" |
| RA3 | 3'-ddC | 🟢 blocked, cannot extend or ligate onward |
| RTP | 5'-biotin | 🟢 "protects it against the exonuclease" |
| RTP | 3'-terminal rA | 🟢 stated; no reason given 🔴 |
| mask | 3'-biotin | 🟢 "to block its end from undergoing extension" |

### Temperatures that fall out of the sequences

✅ Nearest-neighbour Tm (SantaLucia 1998, `lib/chemdraw.tm`, 0.5 µM / 50 mM Na⁺):

| Duplex | Tm | Published temperature it explains |
|---|---|---|
| RTP · RA3 (21 bp) | **59.8 °C** | 🟢 PCR 1 anneals at **60 °C** — within 0.2 °C |
| RP1 full length (50 nt) | 71.5 °C | 🟢 PCR 2 anneals at **67 °C**, below it |
| RP1's 21-nt annealing foot | 54.0 °C | the first cycle of PCR 1 is the hard one; after that RP1 binds full-length |
| masking oligo (76 nt) | 82.9 °C | 🟢 lysis is **20 min at 72 °C** — the mask stays annealed throughout |

🟡 This is the cleanest internal consistency in the paper: in PCR 1 the reverse primer is
the *leftover RT primer*, a bare 21-mer, and the annealing temperature is its Tm. In PCR 2
both primers carry full-length tails, and the annealing temperature goes up 7 °C.

## 3. The masking oligonucleotide

🟢 *"interference from the highly abundant 5.8S ribosomal RNA (rRNA) is avoided by using a
masking oligonucleotide that targets the 3′ end of the rRNA molecule"*, and 🟢 it *"binds
to the 3′ end of the highly abundant 5.8S rRNA, preventing ligation of the 3′ adapter to
these (interfering) rRNA molecules."*

The mechanism is subtractive by *omission*: an rRNA whose 3' end is inside a 76-bp duplex
cannot accept RA3, and without RA3 there is no RT primer site, so it never enters the
library. Nothing is pulled down and nothing is destroyed. 🟢 The authors attribute the
method's sensitivity to exactly this: *"probably achieved by the use of oligonucleotide
masking of rRNA instead of commonly used rRNA pulldown removal strategies."*

✅ What the mask must be annealing to — its reverse complement, 76 nt:

```
ATCATCGACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGCTTGCCGAT
```

🔴 **The paper never prints the rRNA.** Checking this base-by-base needs a 5.8S reference
sequence, which is *data* and so does not live in this repo (README, *Data policy*). It is
recorded above so the next person can BLAST it rather than retype it. ✅ The mask covers
about half of a 156-nt 5.8S, and ✅ it is not an Illumina sequence in disguise.

🟢 Two scope limits, both stated: *"The current rRNA-masking oligonucleotide is designed
for human and mouse cells"*; and *"no other rRNA species appeared in substantial amounts in
the final library, so we did not design any additional masking oligos for them."*

🟢 Omitting the mask is also the positive control for a cell having been deposited: a peak
appears at **270–290 bp** in cell-containing wells. See §6 — that number is reproduced by
the segment arithmetic.

## 4. Steps

| # | Paper steps | What happens | Note |
|---|---|---|---|
| 1 | 1–7 | FACS one cell into 3 µl lysis buffer: 0.13% Triton X-100, RNase inhibitor, **2 µM masking oligo** | 🟢 plate storable at −80 °C for 6 months |
| 2 | 8–9 | **72 °C, 20 min** | 🟢 lyses, unfolds RNA, anneals the mask; ruptures the nucleus too (snoRNAs are recovered) |
| 3 | 10–12 | **3' ligation**: RA3, **T4 RNA ligase 2 truncated KQ**, 8% PEG8000, 30 °C 6 h → 4 °C ≥10 h | 🟢 no ATP, so RNA–RNA ligation is impossible; 🟢 works through 2′-O-methyl 3' ends (piRNAs, plant miRNAs) |
| 4 | 13–15 | **Digestion of unligated RA3**: 5′-deadenylase then **lambda exonuclease**, with RTP added | 🟢 deadenylase converts 5'-rApp → 5'-P; RTP anneals; exonuclease eats the phosphorylated strand; RTP's biotin protects it |
| 5 | 16–18 | **5' ligation**: RA5, **T4 RNA ligase 1 + ATP**, 37 °C 1 h | 🟢 needs a 5'-phosphate, which capped mRNA lacks — so mRNA is excluded by chemistry |
| 6 | 19–21 | **RT**: SuperScript II, 42 °C 1 h, in **Taq buffer** | 🟢 Taq buffer because carried-over MgCl₂ is already in excess; RTP was added back at step 4 |
| 7 | 22–24 | **PCR 1**: RP1 only + leftover RTP as reverse primer, 13 cycles, anneal **60 °C** | ✅ = Tm(RTP) |
| 8 | 25–31 | **PCR 2**: RP1 + one **SRX** index primer per well, 13 cycles, anneal **67 °C** | 🟢 1 µl of PCR 1 as template |
| 9 | 32–35 | Pool, purify, optional size selection (Pippin **130–160 bp**, or 10% TBE gel cut **120–200 bp**) | 🟢 target **145–155 bp** |
| 10 | 36–38 | Sequence, **single read ≥51 bp** | see [`02_sequencing_primer.md`](02_sequencing_primer.md) |

The ordering of §4 steps 3–5 *is* the protocol. 🟢 *"If this procedure is not implemented,
the adapter dimer fraction will dominate the sequence reads."* Shared write-up:
[`../ref/concepts/small-rna-ligation.md`](../ref/concepts/small-rna-ligation.md).

## 5. The UMI

🟢 *"we introduced a custom-designed unique molecule identifier (UMI) into the 5′ adapter"*.

- 🟢 **8 positions, each `rH` = rA, rU or rC.** Not `N`: *"The H mix was chosen over the
  'N' mix … because it results in less unwanted primer annealing."* ✅ 3⁸ = **6,561**
  UMIs, the paper's own number — against 4⁸ = 65,536 for an 8-mer N.
- 🟢 Followed by **`rCrA`**, which *"serves as a linker between the small RNA insert and
  the UMI."*
- 🟢 **"Sequencing starts from the UMI"**, so it is read before any biology: ✅ offset
  **55** in the final library, immediately after P5 + RA5.
- 🟢 Trim 8 + 2 = **10 nt** before mapping; ✅ a 51-bp read leaves 41, and the paper's own
  worked example of a 76-bp read leaves 66. Both reproduced.
- 🟢 Deduplication uses the **adjacency method** (UMI-tools), requiring Hamming ≥ 2.

🟡 **Where the UMI sits is the whole design.** In SMART-seq3 the UMI rides on a TSO and
only labels molecules that template-switched. Here it is ligated onto *every* molecule that
got a 5' adapter, and is read first — so every read carries a count. The cost is that a
UMI placed before the insert cannot distinguish molecules by position: counting depends on
🟢 *"reads that have exactly the same sequence and length"* plus the UMI.

🟡 6,561 is small. The paper calls it *"in principle, high enough to count all the
molecules transcribed from a single locus"*; for an abundant miRNA in a single cell that is
an assumption, not a result, and saturation would bias counts downward — exactly the
failure mode addressed by lineage-UMI collision models.

## 6. Final library, and where the sizes come from

```
5'- P5 · RA5 · [8-nt UMI] · CA · small RNA · RA3 · ?? · [8-bp i7] · revcomp(P7) -3'
```

✅ Published left arm, base for base, 86 bp:

```
AATGATACGGCGACCACCGAGATCTACAC GTTCAGAGTTCTACAGTCCGACGATC HHHHHHHH CA … TGGAATTCTCGGGTGCCAAGG
```

🔴 **The right arm is not published.** SRX's sequence is in Supplementary Table 1. What
*is* published: 🟢 8-bp barcodes, 🟢 192 of them, 🟢 *"modified from standard Illumina
TruSeq small RNA index primers"*. Its 3' end must anneal where RTP anneals — over RA3 —
so what it *adds* beyond RA3 is P7 + index + some linker.

🟡 That added length can be bracketed from two independent published sizes, which is done
in `tools/smallseq.py` as a flippable constant (`SRX_LINKER_NT`):

| Observation | Published | Implied right arm |
|---|---|---|
| miRNA cut window, insert ≈ 22 nt | 🟢 145–155 bp | 37–47 nt |
| 5.8S peak with the mask omitted, insert ≈ 156 nt | 🟢 270–290 bp | 28–48 nt |

✅ The windows overlap (they did not have to): **37–47 nt**. A TruSeq-small-RNA-shaped arm
— P7 (24) + index (8) + a 12-nt link — is **44 nt**, inside both. Taking that:

| Insert | Library | Published check |
|---|---|---|
| 0 (adapter dimer) | **130 bp** | — |
| 18 nt (shortest analysed) | 148 bp | — |
| 22 nt (typical miRNA) | **152 bp** | 🟢 cut target 145–155 bp ✅ |
| 40 nt (longest "small RNA") | 170 bp | 🟢 profile 100–300 bp ✅ |
| 156 nt (5.8S, mask omitted) | **286 bp** | 🟢 peak 270–290 bp ✅ |

🟡 **And this explains the paper's own complaint.** It lists as a limitation 🟢 *"the
abundance of adapter dimers in the sequencing reads, despite the effects of the digestion
step"*. The dimer is **130 bp** and the Pippin window 🟢 *starts* at 130 bp; the gel cut
🟢 120–200 bp contains it outright. The shortest real library is only **18 bp** away. On a
3% agarose cassette that separation is not achievable, so size selection can *reduce*
dimers but never remove them — which is exactly what the paper reports.

## 7. What this chemistry cannot do

- 🟢 **Ligation bias at both ends.** *"Small-seq is, however, biased toward certain sequence
  compositions at the ends of small RNAs"* — intrinsic to ligation; the random UMI
  *"may reduce the bias at one end."*
- 🟢 **It captures degradation products**, *"specifically those from tumors and primary
  cells"*, which map and annotate like genuine small RNAs.
- 🟢 **It captures contaminants from the enzymes**, so a no-cell control is mandatory and
  its Bioanalyzer trace *"looked almost identical"* to a single cell's.
- 🟢 **Plate scale only**: 192 samples per index set, one cell per well. The authors name
  the fix themselves — *"adapt dual-sample indexing and nano-well plates"*.
- 🟢 **Sensitivity**: *"40% of the miRNA genes detected in 1 μg of RNA can be detected in
  single cells"*, at 🟢 1–2 M reads per cell.

## 8. Unresolved

1. 🔴 **SRX.** Sequence in Supplementary Table 1, absent from the extracted text. The
   8-bp-index, 192-primer, "modified from TruSeq small RNA" description is all there is.
2. 🔴 **The sequencing primer is never named.** See
   [`02_sequencing_primer.md`](02_sequencing_primer.md) — this is the one that matters.
3. 🔴 **The mask's rRNA target is never printed.** Its reverse complement is recorded in
   §3; verifying it needs a reference sequence, which is data.
4. 🔴 **Why RTP's 3'-terminal base is a ribonucleotide (rA)** is stated but not explained.
   🟡 A 3'-rA is a plausible RNase-H2 / cleavable handle, but the paper does not say so and
   nothing downstream uses it.
5. 🟡 **5.8S is taken as 156 nt** to do the size arithmetic. The paper does not give a
   length. 150–160 all land inside the published 270–290 bp window, so nothing in §6
   depends on the exact figure.
