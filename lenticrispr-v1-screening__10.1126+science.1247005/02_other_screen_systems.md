# Other pooled-screening systems — and why their primers don't interchange

> 🟢 verbatim · 🟡 derived · ✅ computed from a real map · 🔴 not published.
> Companion to [`01_lenticrispr_gecko.md`](01_lenticrispr_gecko.md).

---

## 1. The organising fact: four labs, four answers to "which end carries the index"

```
GPP        P5(U6, UNINDEXED)      + P7(cPPT/EFS, 8-nt i7)
Weissman   P5(mU6, 6-nt index)    + P7(scaffold, UNINDEXED)
Bassik     P7(mU6, UNINDEXED)     + P5(scaffold, 6-nt index)
Sabatini   P5(scaffold)           + P7(U6, 5-6 nt barcode)
Chen 2015  P5(U6, 8-nt barcode)   + P7(cPPT, 8-nt barcode)      <- DUAL, up to 144 samples
```

🟡 Nothing about this is arbitrary: each lab's choice follows from its vector's scaffold and
promoter, so a readout is effectively welded to its library.

## 2. Scaffolds are the hard incompatibility

| Scaffold | Where | Diagnostic head |
|---|---|---|
| **Original (Cong/Ran)**, 76 nt | lentiCRISPR v1/v2, lentiGuide-Puro ✅, Sabatini libraries | `GTTTTAGAGCTAGAAATAGCAAG` |
| **F+E** (A-U flip + extended stem) | Broad GPP vectors, pXPR_011 ✅ | `GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG` |
| **CR1 / cr3** (Chen 2013 optimised) | Weissman CRISPRi/a, Bassik pMCB320 | `GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAG` |
| **SAM**, 136 nt with 2× MS2 | lentiSAMv2, lenti sgRNA(MS2) | MS2 stem-loops inserted into the scaffold |

🟢 The two custom Read-1 primers in common use are **mutually exclusive by scaffold**: the
Sabatini primer's reverse complement is an exact prefix of the *original* scaffold and is
absent from CR1; the Bassik primer's is an exact prefix of *CR1* and absent from the
original. So two libraries both loosely called "Sabatini-lineage" cannot share a readout.

🟢 The Weissman vectors (#60955, #66217) match **none** of the usual probes — mouse U6 plus
the CR1 scaffold ✅ — and clone with **BstXI + BlpI**, not BsmBI. Their two `CGTCTC` hits sit
inside PuroR and TagBFP and are useless for cloning. ✅

## 3. Single-cell screens: two ways to read a guide out of a cell

**CROP-seq** puts the guide cassette in the **3' LTR**, which is duplicated during reverse
transcription, so the guide ends up inside a Pol II mRNA and is captured by poly(A) like any
transcript. No modification to the bead chemistry; the cost is that the guide must be
enriched by PCR from the cDNA. ✅ CROPseq-Guide-Puro (#86708) carries the original scaffold,
annotated on the **minus strand** of the deposited map — a motif search that only checks the
plus strand will wrongly report it absent.

**Direct capture Perturb-seq** instead modifies the **scaffold** to carry a capture sequence
(CS1/CS2) complementary to a sequence on the 10x gel bead, so the guide is captured directly
rather than via poly(A). Different trade: needs the right bead chemistry, but no enrichment PCR.

## 4. Systems worth knowing

- **Broad GPP libraries** — Brunello (KO), Brie (mouse), Calabrese (CRISPRa), Dolcetto
  (CRISPRi), Humagne/Inzolia (Cas12a). All read out with the ARGON/KERMIT/BEAKER scheme in
  `01_*.md`. 🟡 Humagne uses 23 nt spacers, Inzolia 20 nt.
- **SAM** (Konermann 2015) — 136 nt scaffold with two MS2 loops; its reverse readout primer
  spans the junction into the second MS2 insertion, so it is **absent from all GeCKO
  vectors** and the two reverse sets cannot be swapped. 🟢
- **Sabatini/Wang** — the only protocol here needing **two** custom sequencing primers
  (Read 1 *and* index). 🟢
- **CHyMErA** — Cas9 + Cas12a hybrid guide on one transcript; solves the diversity problem
  with **dark cycles** (29 then 20) rather than a stagger. 🟢
- **Prime- and base-editing sensor screens** — a synthetic target locus travels in the same
  amplicon as the pegRNA, so the edit and the guide that caused it are read together. Gould
  & Sánchez-Rivera use a **custom four-read Illumina recipe**, one primer per read. 🟢
- **In vivo (Chen 2015)** — dual 8 bp barcodes on both P5 and P7 for up to 144 samples, and
  a nested first PCR because tissue gDNA is dirtier. ~200 µg gDNA and ≥67 PCR wells *per
  mouse per organ* 🟢 — well count scales with animals × tissues, not library size.

## 5. Practical rule

Before adopting any published primer pair, **compute the product against your actual
vector**. `lib/plasmid.py` does this:

```python
from plasmid import read_genbank, amplify
from crispr import clone_guide
v = clone_guide(read_genbank("ref/plasmids/addgene-52961_lentiCRISPRv2.gb"), "GTCGCTGAGTACTTCGAAAT")
amplify(v, fwd, rev)[0].length
```

✅ That is how we found that BEAKER's site is present in lentiGuide-Puro but gives a 554 bp
product, and that KERMIT on lentiCRISPRv2 "works" only in the sense of producing a 12.8 kb
wrap-around. Site presence is the wrong test; product length is the right one.

## 6. Not published 🔴

Cuella-Martin's 6 bp UDI sequences; Sánchez-Rivera's 93 nt improved scaffold (described,
never printed); the exact tracrRNA nucleotides of pRDA_256/pRDA_078; Adamson 2016 CR2/CR3;
10x's kit oligos; Chen 2015's PCR#2 cycle number; pMCB320 and pRDA_550 full sequences.
