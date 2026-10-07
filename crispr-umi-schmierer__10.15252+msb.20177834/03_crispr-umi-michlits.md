# CRISPR-UMI (Michlits) — a 10-nt barcode cloned in before the guide

> 🟢 verbatim from the paper · 🟡 derived · ✅ computed from the deposited vectors · 🔴 not published.

**Michlits G, Hubmann M, Wu S-H, Vainorius G, Budusan E, Zhuk S, Burkard TR, Novatchkova M,
Aichinger M, Lu Y, Reece-Hoyes J, Nitsch R, Schramek D, Hoepfner D, Elling U.**
"CRISPR-UMI: single-cell lineage tracing of pooled CRISPR–Cas9 screens."
*Nat Methods* 2017;**14**(12):1191–1197. doi:[10.1038/nmeth.4466](https://doi.org/10.1038/nmeth.4466)
PMID 29039415. Paywalled; PDF in `../pdf/`. Step-by-step protocol at Protocol Exchange (ref 51).

---

## 1. ⚠ Three methods, same name

| | **Michlits (this page)** | **Schmierer** ([`01_`](01_crispr-umi-schmierer.md)) | **CRISPR-MIP** ([`../crispr-mip__10.1101+2024.03.28.587082/`](../crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md)) |
|---|---|---|---|
| label | **10-nt barcode**, cloned in a *separate step before* the guide | **6-bp RSL**, in the same oligo as the guide | **13-nt UMI**, on a capture probe |
| what counts as the UMI | 🟢 the **sgRNA–barcode pair** | the RSL | the probe's UMI |
| labels | a transduced **cell lineage** | a transduced **cell lineage** | a captured **molecule** |
| changes the vector | **yes** | **yes** | no |
| bias strategy | **enrich 10³–10⁴× before PCR**, so few cycles are needed | conventional PCR | replace the PCR entirely |

🟡 Michlits and Schmierer are the same *idea* — a lineage label in the library — reached
independently and published two months apart in 2017. They differ in label length, cloning
order, and above all in how they fight PCR bias.

## 2. The conceptual point 🟢

> "the depletion limit in CRISPR screens applies to analyses of populations… Evaluation of
> independently derived single-cell clones as biological replicates in a pooled CRISPR
> screen overcomes the depletion limit by discerning clones with homozygous LOF alleles that
> are depleted from clones with other (non-LOF) alleles that survive."

The argument is about **editing-outcome heterogeneity**, not just counting noise: one sgRNA
produces LOF alleles in some cells and in-frame or neutral alleles in others, and a
population read-count averages over both. Scoring per clone makes the readout binary and
"**a true reflection of biological effect, independent of editing efficiency**".

🟢 Two requirements are stated explicitly, and the second is easy to miss:

1. each infected cell must be individually trackable by a UMI;
2. **the population must be carried through a strong bottleneck**, so that at most one cell —
   and therefore one editing outcome — per UMI remains.

Without (2) a UMI marks a mixture of editing outcomes and the binary logic collapses.

## 3. Library construction — barcode first, then guide 🟢

**Step 1 — the barcode.** A gBlock carrying the Illumina i7 primer binding site, a **10-bp
random stretch** and the P7 adaptor is cloned into the retroviral backbone (gBlock cut
EcoRI + MfeI; vector cut XbaI + MfeI, a 1.5-kb stuffer removed). Ten parallel ligations,
electroporation complexity ≈ **1 million**.

```
acgatgagcagagccagaaccagaaggaacttgactctaga
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC NNNNNNNNNN gtcctcatctgagagctactcatcaacggt
ATCTCGTATGCCGTCTTaTGCTTG TTAATTAA GAATTC ctggacga
\__ Illumina i7 primer site __/ \_ 10-nt _/                \__ P7 __/ \_PacI_/ \EcoRI/
```

🟢 Note the lower-case `a` in the P7 adaptor: *"we changed C to A in the P7 adaptor sequence
to eliminate a BbsI restriction site in the adaptor for library cloning, but reintroduced
the C during PCR"* — the readout primer carries the canonical P7, so the final library is
correct. ✅ The deposited vectors do carry the modified version.

**Step 2 — the guides.** sgRNA subpools (CustomArray, amplified 10 cycles) are cloned into
the barcoded vector by **BbsI Golden Gate** (20 cycles of 37 °C 5 min / 16 °C 5 min).

> 🟢 **This order is what makes the method work.** Each ligation event joins a guide to a
> *different* barcode, so "each sgRNA–barcode pair represented a UMI". Coverage:
> **954–8,776 independent cloning events per sgRNA**; overall complexity **83.5 million**,
> which the authors note exceeds the number of clones any screen will assay.

**Guides are excluded** 🟢 if they contain `GAAGAC` (BbsI), `GTCTCC`, `CTCGAG` (XhoI),
`CGTCTC` (BsmBI) or `GAGACG`; start with `AAGAC`; or end with `CTCGA`. 🟡 The BbsI and BsmBI
exclusions are forced by the two-step Golden Gate cloning — a guide containing either site
would be cut during assembly.

**Library**: 6,560 mouse genes, 4 sgRNAs/gene (5 for druggable), 112 non-targeting controls.
Library skew 4-fold between the 10th and 90th percentiles.

## 4. The readout — enrich first, then a short PCR 🟢

This is the distinctive part, and the reason to read this paper even if you use another method.

> "in a genomic DNA prep, the sgRNA cassette makes up only about **0.1 p.p.m.** of the total
> DNA. We improved this ratio by **three to four orders of magnitude** by flanking our sgRNA
> cassette with **PacI sites** and carrying out a size-selective precipitation with digested
> genomic DNA."

| | |
|---|---|
| 1 | lyse 170 million cells per condition; phenol-extract gDNA |
| 2 | digest with **PacI** |
| 3 | **size-selective precipitation** (SPRI beads) — the cassette lands on a **589-bp PacI fragment**, the control locus on a 7.7-kb one |
| 4 | PCR, **200 × 50 µl reactions per sample** |

🟢 *"This enrichment minimised the number of PCR cycles required for amplification, thereby
reducing PCR amplification biases."* ✅ Verified: the deposited pLenti-UMI has **exactly two
PacI sites**, giving a 596-bp fragment — the paper's 589 bp for its own construct.

**Readout primers** 🟢

```
fwd  AATGATACGGCGACCACCGAGATCTACAC NNNNNN CGAGGGCCTATTTCCCATGATTCCTTC
     \________ Illumina P5 ______/ \ 6-bp / \______ anneals in U6 ______/
                                  experimental index
rev  CAAGCAGAAGACGGCATACGAGAT ACCGTTGATGAGTAG
     \_____ Illumina P7 _____/ \__ anneal __/
```

🟢 *"Our vector design enabled **direct sequencing of the PCR product**"* — the PCR output is
already a complete P5→P7 library; there is no separate library prep. 🟡 Note this is **not**
PCR-free; it is a single, short PCR on a template enriched 10³–10⁴-fold beforehand.

qPCR primers (579-bp amplicon): fwd `AATGATACGGCGACCACCGAGATCTACACGAGTGGCGAGGGCCTATTTCCCATGATTCCTTC`,
rev as above.

## 5. Sequencing 🟢

```
A5 — [6-bp index] — U6 — sgRNA — i7 site — [10-nt barcode] — A7
```

HiSeq 2500, **single-read 50, dual-indexing**. The sgRNA is read with a **10×-concentrated
custom read primer**:

```
CGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG
```

🟢 *"it is necessary to obtain at least **10 bp for index 1 (barcode)** and **6 bp for
index 2 (experimental index)**."*

> ⚠ **The paper contradicts itself here.** Figure 2c labels them the other way round —
> "Index 2: barcode, Index 1: experimental index". 🟡 The Methods text is the one to follow:
> in the construct the 10-nt barcode sits next to **A7**, and the index read adjacent to P7
> is Index 1 by definition, so barcode = Index 1. The figure legend's labels are swapped.

## 6. The deposited vectors ✅

| | Addgene | bp | guide cloning | scaffold |
|---|---|---|---|---|
| **pLenti-UMI** | [#222694](https://www.addgene.org/222694/) | 10,694 | **BsmBI**, 2 sites, `CACC`/`GTTT` | **F+E** |
| **pRetro-UMI** | [#222686](https://www.addgene.org/222686/) | 6,991 | **BbsI**, 2 sites, `CACC`/`GTTT` | **F+E** |

Deposited by Ulrich Elling, cited to this paper. Both carry the built-in i7 adapter, the
BbsI-modified P7 and the flanking PacI sites. ✅ The CRISPR-MIP extension arm (in U6) is
present in both; its ligation arm is not.

> ⚠ **These are 2024 deposits and differ in detail from the 2017 construct.** ✅ The gBlock's
> `gtcctcatctgagag` linker is absent, and the slot between the i7 adapter and P7 holds a
> 6-nt placeholder rather than the 10-nt barcode. Architecture matches; exact sequence does
> not. **Do not compute the paper's numbers from these maps** — get the construct you
> actually have.

🟡 Both use the **F+E scaffold**, not the AU-flip of the Schmierer vector and not the
original — a third scaffold variant across the three UMI-screening methods.

## 7. What to take from it

- 🟢 **The PacI enrichment is the transferable idea.** Flanking the cassette with a rare
  cutter and size-selecting before PCR is cheap, needs no new chemistry, and attacks the same
  bias CRISPR-MIP attacks with padlocks. It is worth considering for any screen readout.
- 🟢 **The bottleneck requirement is a real experimental constraint**, not a detail: without
  it, one UMI tags a mixture of editing outcomes.
- 🟡 **The UMI is the sgRNA–barcode pair**, not the barcode alone — 10 nt gives 10⁶
  combinations, far short of 83.5 million on its own.
