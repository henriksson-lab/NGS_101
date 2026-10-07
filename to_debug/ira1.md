# ira1 — debugging the lentiGuide/CRISPR-MIP readout

> **Status: likely root cause found.** The core facility said *"It's the Lenti-Puro
> plasmid, it's not yet in the preprint"*. That is almost certainly **not** plain
> lentiGuide-Puro but **pLenti-Puro-AU-flip-3xBsmBI** from the Schmierer CRISPR-UMI
> paper — see §7, which changes the answers in §3 and §4. Sections 1–6 stand as the
> analysis of the maps we were actually given.

# What the two maps and the primers actually do

All numbers below are computed from the two `.dna` files in this folder (they are GenBank
under another extension) using `lib/plasmid.py`. Reproduce with
`python3 to_debug/check.py`.

---

## 1. Both maps you sent are EMPTY vectors

Each still has **two BsmBI sites 1885 bp apart** — the filler is in. So any size read off
these maps is the *uncloned* size, which is the root of the confusion below.

| | cloned (guide in) | empty (filler in) |
|---|---|---|
| **lentiCRISPRv2** PCR1 | **288 bp** | 2148 bp |
| **lentiGuide-Puro** PCR1 | **557 bp** | 2417 bp |
| **either**, PCR2 | **351 bp** | 2211 bp |

## 2. Why you can't see the lentiGuide product in the sequence

**The lentiGuide-Puro PCR1 amplicon crosses the origin of the map.** `CRISPR_PCR1-F` binds
at 7959 and `CRISPR_PCR1-R` at 108 on a 10,183 bp circle, so the product runs off the end of
the file and back to the start. Looking at the map linearly in ApE, it looks impossible.
On lentiCRISPRv2 it does not wrap (2814 → 4877), which is why that one reads normally.

## 3. The primers do work on lentiGuide-Puro — but PCR1 is a different size

`CRISPR_PCR1-R` anneals in the **EF-1α / EFS leader**, and that element sits at a different
distance from U6 in the two vectors:

```
cloned lentiCRISPRv2    U6 --- guide --- scaffold --- EFS      PCR1 = 288 bp
cloned lentiGuide-Puro  U6 --- guide --- scaffold --- EF-1a    PCR1 = 557 bp
```

**PCR2 is 351 bp on both**, which matches your gel. That is not luck: `CRISPR_PCR1-R`
carries the cPPT sequence as a **non-templated 5' tail**, so PCR1 *grafts* the PCR2 reverse
primer site onto the product. PCR2 then runs U6 → grafted tail, and both of those are
vector-independent — the stretch that differs between the vectors lies outside the PCR2
amplicon. So `~350 bp PCR2` is correct and expected on either backbone.

> **The single most useful measurement you can make:** size PCR1 precisely.
> **557 bp = correctly cloned lentiGuide library. ~2.4 kb = filler still present**, i.e.
> uncut or re-ligated empty backbone in the prep. A band ">1 kb" points at the second.
> On the attached gel, lanes 1–4 show a doublet and lanes 5–6 only the upper band — worth
> cutting both bands and sequencing, since the model predicts 351 bp and 2211 bp are the
> two products PCR2 can give.

## 4. Why the MIP probes failed on plasmid

Three reasons, and **none of them is "plasmids need special probes"**. The first two are
computed from your own maps.

### (a) The ligation arm does not exist in lentiGuide-Puro

The probe's ligation arm is `AGCTAGGTCTTGAAAGGAGTGGG`. Immediately after the scaffold
terminator:

```
lentiCRISPRv2     TTTTTTGAATTCGC AGCTAGGTCTTGAAAGGAGTGGG AA     <- arm present
lentiGuide-Puro   TTTTTTAAGCTTGGCGTAACTAGATCTTGAGACAAATGGC      <- arm absent
```

The probe was designed against lentiCRISPRv2 and **cannot close on lentiGuide-Puro at all**.
The extension arm is fine in both (it sits in U6); it is the ligation arm that is missing.

### (b) On an empty vector the gap is 1972 nt, not 112

The arms sit either side of the cloning site, so the filler is *inside* the captured region.
Uncloned, the polymerase would have to fill **1972 nt** in a 1 h extension. It will not, so
no circle forms and the exonuclease destroys everything.

### (c) A supercoiled plasmid barely denatures *(inference, not computed)*

A covalently closed circular plasmid has topologically interlinked strands: heat it and it
reanneals essentially instantly on cooling, so the arms never get a window in which to bind.
Genomic DNA is already fragmented, which is why the same probe works there. If you want to
run MIPs on plasmid, **linearise first** with a single cutter outside the captured region
(or nick it) before the hybridisation step.

## 5. If you want a lentiGuide-compatible probe

Keep the extension arm, replace the ligation arm. Candidates scanned from the real map,
23 nt, unique, and with a T<sub>m</sub> above the extension arm's 54.1 °C as the design
requires:

| candidate ligation arm | Tm | gap once cloned |
|---|---|---|
| `CCGAGTCGGTGCTTTTTTAAGCT` | 58.8 | 81 nt |
| `GTCGGTGCTTTTTTAAGCTTGGC` | 58.6 | 85 nt |
| `GGTGCTTTTTTAAGCTTGGCGTA` | 57.5 | 88 nt |

Verified end to end on the cloned vector: one capture, gap 89 nt, 215 nt circle, spacer
present. ⚠ These overlap the scaffold's 3' end and the `AAGCTT` (HindIII) just past it, so
they are *partly* lentiGuide-specific.

**A better option if you want one probe for every guide vector:** put the ligation arm
entirely **inside the scaffold**, which is identical in lentiCRISPR v1, lentiCRISPRv2 and
lentiGuide-Puro. That is exactly the trick that makes the Joung readout primer
vector-independent while KERMIT and BEAKER are not. The captured region shrinks to roughly
the +1 G plus the spacer plus part of the scaffold, which is all the readout needs.

## 6. Note on the tooling

Finding this exposed a bug in `lib/plasmid.py`: an amplicon ending on the **last base of a
linear template** had its end coordinate taken modulo the length, sending it to 0, and was
silently dropped. That is why PCR2 first appeared to give no product. Fixed, with
regression checks (`lib/checks.py`, "PCR prediction"). Circular templates were unaffected,
so no earlier result on this project changed.


---

## 7. The core facility's plasmid is probably not the one we modelled

> *"It's the Lenti-Puro plasmid, it's not yet in the preprint"* — core facility

Taken together with the CRISPR-UMI work, this points at
**pLenti-Puro-AU-flip-3xBsmBI**, the vector in

**Schmierer B, Botla SK, Zhang J, Turunen M, Kivioja T, Taipale J.** "CRISPR/Cas9 screening
using unique molecular identifiers." *Mol Syst Biol* 2017;**13**(10):945.
doi:[10.15252/msb.20177834](https://doi.org/10.15252/msb.20177834) · PMC5658704 · CC BY

🟢 Methods, verbatim: it "was created by modifying lentiGuide-Puro (a gift from Feng Zhang,
Addgene #52963) by replacing the sequence `gttttagagctagaaatagcaagttaaaa……TTTTTT` with
`gtttAagagctagaaatagcaagttTaaa……TTTTTTcgtctct` to create an AU-flip (Chen *et al*, 2013)
and an additional BsmBI site downstream of the tracrRNA."

**So it is lentiGuide-Puro with three changes**, and each one breaks something different:

| change | consequence |
|---|---|
| **AU-flip scaffold** — 2 substitutions ✅ (positions 5 `T`→`A` and 26 `A`→`T`) | 🔴 **anything annealing to the scaffold's 5' half no longer matches.** ✅ `GTTTTAGAGCTAGAAATAGCAAG` is absent from this vector |
| **a third BsmBI site** after the terminator | cloning behaves differently; a 2-site digest model is wrong |
| **the library insert carries the Illumina Read-2 adapter and a 6-bp RSL** | 🔴 **the readout PCR is a different design entirely** — see below |

### The insert, as published 🟢

```
ggctttatatatcttgtggaaaggacgaaacaccg [20-nt spacer] gtttaagagctagaaatagcaagtttaaa…tttttt
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC nnnnnn aagcttggcgtaactagatcttgagacaaa
\___ Illumina Read-2 / Index-1 adapter ___/ \RSL/ \____ lentiGuide, past HindIII ____/
```

**The Illumina adapter is built into the construct**, and the 6-bp RSL sits exactly where an
i7 index would — so it is read by the Index 1 read. You do not add the Read-2 site by PCR on
this vector; it is already there. Any readout designed to graft one on (which is exactly what
`CRISPR_PCR1-R` does) is solving a problem this vector does not have, and will give a
different product.

### What this predicts for our observations ✅

Rebuilding the vector from the real lentiGuide-Puro map and applying the published edits:

| | plain lentiGuide-Puro | pLenti-Puro-AU-flip-3xBsmBI |
|---|---|---|
| PCR1 with `CRISPR_PCR1-F/R`, cloned | 557 bp | **596 bp** (+39) |
| Schmierer's own readout (PCR1→2→3) | 249 bp | **288 bp** (+39) |
| original scaffold `GTTTTAGAGCTAGAAATAGCAAG` | present | **absent** |
| MIP ligation arm | absent | **absent** |
| Illumina Read-2 adapter in the vector | no | **yes** |

The +39 is the same 39 bp in both rows — the 33-nt Illumina adapter plus the 6-nt RSL that
Schmierer inserted. ✅ 288 bp is the figure the paper itself quotes, reproduced from the
rebuilt map, so this row is a validated prediction rather than a model artefact.

### ⚠ The primers in `Primers combined from Martin.xlsx` ARE Schmierer's PCR1 pair

This is the strongest evidence yet that the plasmid is Schmierer-lineage. ✅ Comparing the
workbook primers with the Schmierer Methods:

```
CRISPR_PCR1-F  AATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCG
Schmierer         GGACTATCATATGCTTACCGTAACTTGAAAGTATTTCG   PCR1_FW
                  ^^^ identical over all 38 nt; ours just has 3 extra 5' bases

CRISPR_PCR1-R  TCTACTATTCTTTCC CCTGCACTGT TGTGGGCGATGTGCGCTCTG
               |<-- 5' tail -->|          |<-- 20 nt, anneals -->|
Schmierer      CTTTAGTTTGTATGTCTGTTGCTATTATGTCTACTATTCTTTCC   PCR1_REV
                                            ^^^^^^^^^^^^^^^ the same 15 nt
```

So `CRISPR_PCR1-F` **is** Schmierer's `PCR1_FW`, and `CRISPR_PCR1-R`'s non-templated 5' tail
**begins with the 15 nt that form Schmierer's `PCR1_REV` 3' end**. 🟡 The natural reading is
that this primer set was adapted from the Schmierer protocol — whoever built it kept the
forward primer, and reused part of the reverse primer as a grafted tail instead of as an
annealing sequence.

That makes it very likely the plasmid was supplied *together with* this primer set, i.e. both
came from the Schmierer protocol, and the mismatch is that §1–6 modelled them against
lentiCRISPRv2/lentiGuide instead.

### ✅ Correction: the published readout is NOT "one PCR off the built-in adapter"

An earlier version of this note said that. It is wrong on both counts. The published readout
is **three nested PCRs** (14 / 19 / 14 cycles, 200 µg gDNA across 40 parallel PCR1 reactions),
and ✅ **none of its primers anneal in the built-in Illumina adapter** — checked against all
six. `PCR2_FW` anchors in U6 and `PCR2_REV` downstream of the RSL, so the adapter sits
*inside* the amplicon rather than serving as a primer site.

The adapter's job is on the flowcell, not in the PCR: ✅ its reverse complement is the
standard TruSeq Read-2 / i7 index-read primer sequence (less its terminal T), and it sits
immediately 5' of the RSL. That is what lets the 6-bp RSL be read as an **index read** rather
than needing a read of its own. 🟢 The paper: "`20 + 6 + 6` cycles … where **i7 reads the RSL
and i5 the illumina sample index**" — note that assignment, it is the reverse of the usual
habit and an easy way to lose the lineage labels at demultiplexing.

### What to do next

1. **Sequence the plasmid**, or at minimum Sanger across the scaffold and the region just
   past the terminator. Two diagnostics settle it immediately:
   `GTTTTAGAGCTAGAAATAGCAAG` (original) vs `GTTTAAGAGCTAGAAATAGCAAG` (AU-flip), and whether
   `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC` is present.
2. **Ask the facility for the actual map.** We have been reasoning against lentiGuide-Puro,
   and §1–6 show that even *that* map was the empty vector.
3. **Just size PCR3.** Running Schmierer's own three-step readout gives **288 bp** on his
   vector and **249 bp** on plain lentiGuide-Puro. Both templates amplify — which is why the
   mix-up survives a gel — but the 39 bp difference is unambiguous. On lentiCRISPRv2 the
   readout instead **fails at PCR2** (`PCR2_REV`'s site is absent), so a missing PCR2 product
   is itself informative. All three primer pairs are in `../crispr-umi-schmierer__10.15252+msb.20177834/tools/crisprumi.py`.
4. A MIP probe for it **has been designed** — ligation arm
   `AAGCTTGGCGTAACTAGATCTTGAGA`, gap 142 nt, 271 nt circle, capturing the guide *and* the
   RSL. See `../crispr-umi-schmierer__10.15252+msb.20177834/01_crispr-umi-schmierer.md` §6. Note the probe backbone must drop
   its own Read-2 site, because this vector already contains that sequence.

See [`../crispr-umi-schmierer__10.15252+msb.20177834/01_crispr-umi-schmierer.md`](../crispr-umi-schmierer__10.15252+msb.20177834/01_crispr-umi-schmierer.md)
for the full chemistry.

> ⚠ **Name collision worth being explicit about.** "CRISPR-UMI" means two different things.
> **Schmierer** clones the UMI *into the library*, so it labels a transduced cell lineage and
> counting RSLs counts clones. **CRISPR-MIP** carries the UMI *on a capture probe* applied to
> gDNA afterwards, so it labels a captured molecule. Both are "CRISPR + UMI"; only one of
> them changes the plasmid. If the facility supplied a Schmierer-lineage vector while we were
> designing MIP probes against lentiCRISPRv2, that alone explains a lot.

## 8. ✅ "Addgene Brunello plasmid provided by the client" — which pool?

> *"Addgene Brunello plasmid provided by the client"* — client email

This resolves the vector question, because **Brunello is distributed in two different
backbones**, and the readout behaves completely differently in each.

🟢 Addgene, "Human CRISPR Knockout Pooled Library (Brunello)", Doench *et al* 2016
*Nat Biotechnol* (Broad GPP, Root/Doench). Library: **19,114 genes, 76,441 sgRNAs,
1,000 non-targeting controls**, 3rd-generation lentiviral.

| pool | backbone | Cas9? |
|---|---|---|
| **#73179** | **lentiCRISPR v2** (#52961) | one-vector — **expresses Cas9** |
| **#73178** | **lentiGuide-Puro** (#52963) | two-vector — **does NOT**; needs lentiCas9-Blast (#52962) separately |

> ⚠ **The numbering is a trap.** 73178 is the *lentiGuide-Puro* pool and 73179 is the
> *lentiCRISPRv2* pool — the lower number is the one *without* Cas9. "Brunello" alone does
> not identify the plasmid; always quote the pool number.

### ✅ What each pool predicts, computed from the real maps in this folder

| | **#73179** (lentiCRISPRv2) | **#73178** (lentiGuide-Puro) |
|---|---|---|
| PCR1 with `CRISPR_PCR1-F/R` | **288 bp** | **557 bp** |
| PCR2 | 351 bp | 351 bp |
| MIP **extension** arm (U6) | present | present |
| MIP **ligation** arm | present | 🔴 **ABSENT** |
| MIP capture | gap 112 nt ✅ | 🔴 **impossible** |

### 🔴 If it is #73178, the MIP probe cannot work — at all

The core facility's *"It's the Lenti-Puro plasmid"* points at **#73178**. In that backbone the
probe's **ligation arm does not exist anywhere in the vector**, so the padlock has nothing to
ligate to and no circle can form. This is not an efficiency problem to be tuned — it is a
hard incompatibility, and it explains a complete absence of product rather than a weak one.
See §4(a), which found the same thing before the pool was identified.

⚠ **And PCR2 is 351 bp on both**, which is why the PCR2 gel looked normal throughout and gave
no warning. Only **PCR1 discriminates**: 288 bp vs 557 bp. That remains the single most
useful measurement (§2).

> 📄 **2026-10-04 attempt:** the padlock redesign, the Broad protocol and the proposed PCR
> primers are assessed in
> [`ira_attempt_20261004/assessment.md`](ira_attempt_20261004/assessment.md). Headline: the
> new **MS3** padlock for lentiGuide-Puro checks out end to end (272 bp library, as designed),
> but the proposed second PCR builds a **TruSeq Small RNA** library (⚠ primer compatibility
> unverified), and the workbook primers named "Brunello" are the **lentiCRISPRv2** ones.

### What to do

1. **Ask the client for the Addgene pool number**, not the library name. One number settles it.
2. **Size PCR1.** 288 bp → #73179, probe should work. 557 bp → #73178, probe cannot.
3. If it is **#73178**, either
   - re-design the padlock against lentiGuide-Puro — the guide-only option in
     [`../gcbias/datasets/crispr-umi-schmierer.md`](../gcbias/datasets/crispr-umi-schmierer.md)
     uses ligation arm `TCAACTTGAAAAAGTGGCACCGA`, which sits inside the scaffold and therefore
     works on v1, v2, lentiGuide-Puro *and* the Schmierer vector; or
   - switch to the conventional PCR readout, which works on either backbone.
4. 🟡 Note the two-vector pool also requires a separate Cas9 source. If the client supplied
   #73178 without `lentiCas9-Blast`, there may be a second, independent problem — no editing
   at all — which would look like "the screen did nothing" rather than "the readout failed".

🟡 Note also that this is a *different* lineage from the Schmierer vector discussed in §7. Both
are lentiGuide-Puro derivatives, which is why "Lenti-Puro" was ambiguous; Brunello #73178 is
unmodified lentiGuide-Puro, whereas the Schmierer vector adds the AU-flip, a third BsmBI site
and the RSL/Illumina cassette. ✅ They are distinguishable by PCR1 size: 557 bp (#73178) vs
596 bp (Schmierer), and by whether `GTTTTAGAGCTAGAAATAGCAAG` is present (yes / no).

---

## 9. ✅ MIP worked on gDNA but not plasmid — why, and the fix

New information from Ira: cells carrying the Brunello library **and** the plasmid it was
made from. MIP gave a product on genomic DNA (band "a bit weak"), nothing on plasmid.
`Brunello-P5-xx`/`Brunello-P7-xxx` PCR gave a product on **both**.

### The probe is correctly targeted — this settles §8

The PCR result on plasmid resolves the pool ambiguity of §8. Computed from the real maps:

| pool | vector | Brunello-P7 product | MIP lig arm | capture gap |
|---|---|---|---|---|
| **#73179** | lentiCRISPRv2 | **287 bp** | **present** | **112 nt** ✅ |
| #73178 | lentiGuide-Puro | 556 bp | absent | none |

Their stock is **#73179 (lentiCRISPRv2)**. Both arms are present and the gap is the
validated 112 nt, so the MIP failure is **not** the §4(a)/§4(b) sequence problem.
It is physical — and it is the one already predicted in **§4(c)**.

### Why supercoiling kills it

A plasmid straight from a prep is **covalently closed circular**: both strands are
unbroken, so the linking number is topologically fixed and the strands cannot be
separated from each other, only locally unwound. Consequences for our protocol:

1. 94 °C transiently melts base pairs, but nothing becomes a free single strand — the
   two strands stay interlocked.
2. The −0.1 °C/s ramp to 60 °C is ~6 min. Re-pairing of two strands held in the same
   molecule is **intramolecular**: first-order, concentration-independent, and far
   faster than a probe at nM finding its site.
3. By the time the overnight 60 °C hybridisation starts, the target is double-stranded
   again. No stable probe-target duplex → no gap-fill → no circle → exonuclease
   destroys everything → no product.

Genomic DNA is linear and sheared, so its strands *do* separate, and re-annealing is
**bimolecular** — second-order and slow at genomic complexity. The probe wins. That is
the entire difference between the two substrates.

PCR works on both because it never needs a *persistently* single-stranded target: 98 °C
for 10 s each cycle is enough for a primer to extend, and the product is linear from
cycle 1 onward. **So the PCR-vs-MIP comparison on plasmid is not evidence about the
probe — it is evidence about topology.**

### Fix: linearise the plasmid first

An enzyme must cut the plasmid but **not** inside the probe footprint — ext arm (19 nt)
+ captured region (112 nt) + lig arm (23 nt) = **154 bp**. In the cloned #73179 map
(13,013 bp) that footprint spans the GenBank origin, at 12,973–154.

Screened a 29-enzyme panel against the vector, and separately against all **77,441**
real Brunello guides (a site landing in the 20-nt spacer destroys the footprint for
*that* library member):

| enzyme | site | vector cuts | guides hit | % of library |
|---|---|---|---|---|
| **NotI** | GCGGCCGC | **1** | **10** | **0.013%** |
| PacI | TTAATTAA | 1 | 3 | 0.004% |
| SpeI | ACTAGT | 1 | 141 | 0.18% |
| XbaI | TCTAGA | 1 | 216 | 0.28% |
| BamHI | GGATCC | 1 | 389 | 0.50% |
| HindIII | AAGCTT | 4 | 237 | 0.31% |

**Use NotI** (or PacI): single cut, 8-bp site, ≥99.98% of the library untouched.
Unusable — they cut inside the footprint in the vector itself: **EcoRI** (GAATTC),
**NheI** (GCTAGC).

Worth noting both ways: NotI's site is GC-rich so the 10 guides it drops are GC-rich
ones — the exact population our GC-bias work is about — and PacI's AT-rich site drops
AT-rich ones. At 10 and 3 guides out of 77,441 neither matters, but pick NotI for the
library and remember it is not a GC-neutral loss in principle.

**Also:** `HindIII`/`BamHI` — the enzymes the preprint already uses to digest gDNA —
both leave the footprint intact here, so they would work if that is what is on the
bench. The "a bit weak" gDNA band is itself consistent with incomplete accessibility;
digesting the gDNA too should improve it.

### Should it work?

Yes, after linearisation. Nothing about #73179 is incompatible with the probe. Predict:
NotI-linearised plasmid gives a MIP product comparable to gDNA, and because a plasmid
prep is pure target with no 3 Gb of competing sequence, it should be **stronger** than
the gDNA band, not weaker. If linearised plasmid *still* fails, the problem is upstream
of hybridisation (probe stock, ligase, exonuclease) and the plasmid/gDNA contrast is a
red herring.

### Caveat on the numbers above

Footprint coordinates were computed with a placeholder 20-nt spacer. The vector-side
conclusions (which enzymes cut where, the 154 bp span, the 112 nt gap) are independent
of the spacer; the per-guide column is computed from the real Brunello sequences.

---

## 10. ✅ Step-1 PCR size: expected vs. "~1 kb" observed

Observation from Ira: **first PCR of two on the Brunello library gives ~1 kb.**

First, which pair. The workbook's `CRISPR_PCR1-F/R` is already in the repo as
`V2ADAPTOR_F`/`V2ADAPTOR_R` (`lenticrispr-gecko-screen__10.1126+science.1247005/tools/crisprscreen.py:55`) — the **Zhang "v2
adaptor"** step 1. It exists because lentiCRISPRv2 *does* contain the cPPT site Shalem's
original reverse anneals to, but **upstream** of U6, so Shalem's pair points away from the
guide and only yields a 12.7 kb product the long way round the circle. The v2 adaptor
fixes this by grafting the cPPT site on as a non-templated 5' tail instead.

> **Schematics:** the two-step system is drawn base-by-base in
> [`../lenticrispr-gecko-screen__10.1126+science.1247005/crisprscreen.html`](../lenticrispr-gecko-screen__10.1126+science.1247005/crisprscreen.html) &mdash;
> **Part 3, panels (7)&ndash;(10)**: why a second PCR exists, the reverse primer's grafted
> tail, the round-1 and round-2 constructs for both vectors, and the size table. The
> one-step systems are panels (5)&ndash;(6) of the same page.

### Every size these two vectors can give (computed, all in `selftest.py`)

| step-1 pair | template | product |
|---|---|---|
| **v2 adaptor** (the workbook pair) | **cloned lentiCRISPRv2 #73179** | **287 bp** ← expected |
| **v2 adaptor** | **lentiCRISPRv2 + stuffer** | **2,148 bp** |
| v2 adaptor | cloned lentiGuide-Puro #73178 | 556 bp |
| v2 adaptor | lentiGuide-Puro + stuffer | 2,417 bp |
| Shalem F1/R1 | cloned lentiCRISPRv2 | none (12.7 kb wrap-around only) |
| Shalem F1/R1 | cloned lentiGuide-Puro | 316 bp |
| Shalem F1/R1 | lentiGuide-Puro + stuffer | 2,177 bp |

(287 not 288 because these use a G-initial spacer, so no +1 G is appended. The stuffer
figures use the empty Addgene maps directly; both primer sites lie outside the excised
filler, so the amplicon shrinks by exactly as much as the plasmid does — 1,861 bp.)

### So: 🔴 ~1 kb is not any of them

The four reachable sizes are **287 / 556 / 2,148 / 2,417**. There is no ~1 kb product, and
this holds down to `min_anneal` 14 — relaxing primer stringency produces no additional
product of any size on either vector, cloned or stuffered. Three readings, most likely
first:

1. **Wrong reverse primer.** The two step-1 pairs share an *identical* forward primer and
   differ only in the reverse — the single easiest mix-up in this protocol. If the tube
   holds Shalem's `CTTTAGTTTGTATGTCTGTTGCTATTATGTCTACTATTCTTTCC` rather than the v2
   adaptor's `TCTACTATTCTTTCCCCTGCACTGTTGTGGGCGATGTGCGCTCTG`, then on lentiCRISPRv2
   **there is no specific product at all** and any band on the gel is spurious. Check the
   primer sequence before anything else; this costs nothing and explains the observation.
2. **Not leftover stuffer.** Uncut or re-ligated backbone would give **2,148 bp**, not
   1 kb, and a partially cloned prep would show a **doublet at 287 + 2,148** — not a
   single intermediate band. Note the §3 gel observation of a doublet in lanes 1–4; that
   is the pattern to re-examine, with this corrected pair of sizes rather than §3's
   lentiGuide-based ones.
3. **A different template.** If the band is real and specific, neither vector is the
   template: ~710 bp more than lentiCRISPRv2 between U6 and the graft site. Cut the band
   and Sanger it off `CRISPR_PCR1-F`; one read settles vector identity, cloning state and
   guide presence at once.

### The useful measurement

**Size step 1 precisely against a 100 bp ladder, not a 1 kb ladder.** The diagnostic
spread is 287 → 2,148 and the informative region is the bottom of the gel. A ~1 kb call
off a 1 kb ladder is also exactly the error mode that would make a 287 bp band
unresolvable from primer dimer, so the ladder used is worth confirming.
