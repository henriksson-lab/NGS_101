# Enzymes, buffers and cycling, per paper

Why this table exists: **GC bias is strongly enzyme- and buffer-dependent.** Polymerases
differ several-fold in how evenly they amplify GC-rich template (Quail *et al* 2012,
*Nat Methods* 9:10, compared KAPA HiFi, Phusion, AccuPrime and Taq-based mixes on exactly
this axis), and additives — DMSO, betaine, a dedicated GC buffer, 7-deaza-dGTP — change it
again. So before attributing a measured GC effect to "PCR", it is worth knowing which PCR.

Marking: 🟢 verbatim from the paper · 🟡 inferred · 🔴 not stated · ✅ computed/checked here.

## Summary — the amplification steps that can bias

| | **Schmierer 2017** | **Michlits 2017** | **CRISPR-StAR 2025** | **Behan 2019** (Project Score) | **CRISPR-MIP (ours)** |
|---|---|---|---|---|---|
| library oligo pool → insert | "double-stranded using outer primers" 🟢, enzyme 🔴 **not stated** | "ten cycles of PCR" 🟢, enzyme 🔴 **not stated** | **NEBNext High-Fidelity 2× Master Mix** 🟢 | library amplified in **bacteria**, not by PCR 🟢 (Addgene replication protocol) | — |
| cloning into vector | Gibson — **NEBuilder HiFi Assembly** 🟢 (no PCR) | Golden Gate, BbsI + T4 ligase 🟢 | Golden Gate, BsmBI 🟢 | BbsI, pKLV2 🟢 | — |
| capture / fill-in | — | — | — | — | **Phusion HF** (F530S, ThermoFisher), 4 U 🟢 |
| ligation | — | — | — | — | **Ampligase** (Lucigen), 5 U 🟢 |
| readout PCR | **KAPA HiFi HotStart** (KAPA Biosystems) 🟢 | 🔴 **not stated** | **KAPA HiFi HotStart ReadyMix** (Roche) 🟢 | **Q5 Hot Start High-Fidelity 2× Master Mix** (NEB) 🟢 — via Koike-Yusa 2014 | **KAPA HiFi HotStart ReadyMix** (Roche) 🟢 |
| cycles | 14 + 19 + 14 = **47** 🟢 | 🔴 not stated | 🔴 not stated | 🔴 not stated | determined per screen by qPCR 🟢 |
| 🟡 same lab elsewhere | — | — | — | **Phusion HF in GC buffer** for single-locus PCR; KAPA HiFi for adapter-ligation enrichment 🟢 | — |
| qPCR enzyme | — | **GoTaq** (Promega, A6001) 🟢 | — | — | **KAPA HiFi HotStart ReadyMix** + SYBR Green I 🟢 |
| additives (DMSO/betaine/GC buffer) | 🔴 none mentioned | 🔴 none mentioned | 🔴 none mentioned | 🔴 none mentioned | 🔴 none mentioned |
| sequencer | HiSeq 4000 🟢 | HiSeq 2500 🟢 | NextSeq 2000 P2 🟢 | **HiSeq 2000 v4**, 19-bp single-end 🟢 | 🟡 see preprint |

## How much DNA goes into the polymerase reaction

This is the other half of the question, and arguably the more important one: the **number of
template copies per guide in a single reaction** sets how stochastic the first cycles are —
which is precisely the mechanism our own preprint proposes for the GC effect. A guide present
at a handful of copies can lose several cycles to a slow first denaturation and never catch
up; at thousands of copies the average washes that out.

✅ Computed at 6.6 pg per diploid human genome, 6.0 pg mouse:

| | DNA per reaction | reactions | genome-equivalents **per reaction** | total | ✅ **copies per guide per reaction** |
|---|---|---|---|---|---|
| **Schmierer 2017** | **5 µg** 🟢 | 40 🟢 | 757,576 | 200 µg = 30M cells 🟢 | **≈ 33** |
| **CRISPR-StAR 2025** | **4 µg** 🟢 | 48 🟢 | 666,667 | ≈ 192 µg | **≈ 7** |
| **Michlits 2017** | 🔴 not stated | 200 × 50 µl 🟢 | — | 170M cells lysed per condition 🟢 | 🟡 see note |
| **Behan 2019** 🟡 | **1 µg** 🟢 (Koike-Yusa 2014) | 72 🟢 | 151,515 | 72 µg ≈ 1.1 × 10⁷ cells 🟢 | **≈ 1.7** ⚠ |
| *…or* 🟡 | **5 µg** in 50 µl (bench protocol, ⚠ possibly **pre-2014**) | — | 757,576 | — | **≈ 8.4** |
| *…its plasmid library* | **15 ng** 🟢 | 10 🟢 | — | 1.7 × 10¹⁰ molecules 🟢 | — |
| **CRISPR-MIP (ours)** | **1–10 µg** 🟢 per 20 µl | 1 | 151,515 – 1,515,152 | — | **≈ 7 – 65** |

🟡 **Michlits is the interesting case.** 170 million cells are lysed per condition, but the
cassette is enriched **10³–10⁴-fold by PacI digest and size selection before any PCR**. So the
DNA *mass* entering each reaction is small while the *cassette copy number* is preserved —
which is the whole point of the design, and the opposite of every other protocol here, where
cassette copies are diluted in a large excess of carrier genome. 🟢 The paper says so
explicitly: the enrichment "minimized the number of PCR cycles required for amplification,
**thereby reducing PCR amplification biases**."

⚠ **Our own method sits at the low end at 1 µg input** (≈7 copies per guide), which is the
most stochastic regime in the table — and CRISPR-MIP's capture is one-shot, so a molecule
lost at the fill-in is lost permanently rather than merely under-amplified. 🟡 If the GC
effect is driven by first-cycle stochasticity, input mass is a knob worth testing alongside
enzyme: 10 µg gives a ~10× deeper copy number per guide at identical chemistry.

⚠ **Behan/Project Score is the most stochastic regime in the table — at ~1.7 copies of each
guide per reaction.** 🟢 Koike-Yusa: "72 independent PCR reactions using 1 µg of the mouse ESC
library per reaction … These correspond to 1.7 × 10¹⁰ molecules of the plasmid DNA and
1.1 × 10⁷ ESCs in total". ✅ That arithmetic checks out (72 µg ÷ 6.0 pg = 1.2 × 10⁷ cells),
which also settles a typography trap: the PDF renders µ as `m`, so "1 mg" is **1 µg** — the
authors' own molecule count proves it, and the same substitution appears in "3 mg" of
lentivirus and "1 mg of genomic DNA (1.5 × 10⁵ cells)" (= 0.9 µg ✅).

⚠ **But the two "Yusa protocol" sources disagree 5-fold on input, and Behan does not say
which it followed.**

| source | per reaction | ✅ copies per guide | provenance |
|---|---|---|---|
| **Koike-Yusa 2014**, *Nat Biotechnol* 🟢 | **1 µg**, 72 reactions | **≈ 1.7** | the published method Behan reaches via Tzelepis 2016 |
| Yusa-lab **bench protocol** 🟡 | **5 µg in 50 µl** (= 100 ng/µl) | **≈ 8.4** | paper copy, supplied 2026-10-04; ⚠ **may predate the 2014 paper** |

🟡 **The citation chain favours 1 µg.** Behan → Tzelepis 2016 → Koike-Yusa 2014 is an explicit
published trail, and if the bench sheet is older than 2014 then the 1 µg figure **supersedes**
it rather than competing with it. On that reading Project Score ran at ≈1.7 copies per guide
per reaction — the most first-cycle-stochastic regime in the table.

⚠ But this is not settled, for a mundane reason: **labs run off bench sheets, not off their
own papers**, and a protocol sheet can stay in use unchanged for a decade. A 2019 screen
following a pre-2014 bench copy at 5 µg is entirely plausible, and ≈8.4 copies per guide is
unremarkable — it sits beside CRISPR-StAR (7.4) and our own 1 µg condition (6.5).

🔴 So: **do not assert "Project Score is where first-cycle stochasticity bites hardest"**. It
holds on the 1 µg reading only, and the 5-fold spread is unresolved. Dating the bench protocol
would settle it; so would asking the Sanger screening group directly, which is a single email
and worth sending if this line of argument is going into a response to reviewers.

⚠ Note also that 5 µg in 50 µl is **100 ng/µl** of genomic DNA in the reaction — high enough
that viscosity and incomplete denaturation of high-molecular-weight DNA become plausible, which
is the same accessibility argument the PacI-digest protocols make (see the purification
section).

### 🟡 The bench protocol is a TWO-step PCR, and round 2 uses KAPA

🟡 Further detail from the same bench copy: round-1 product is cleaned on **Qiagen PCR
purification columns**, eluted into **50 µl EB**, and **1 ng** of that is carried into a
**50 µl round-2 reaction using KAPA**.

⚠ **This undercuts the clean "Behan = Q5, ours = KAPA" enzyme contrast drawn above.** If
Project Score followed this sheet, its final amplification is **KAPA**, the same enzyme as
Schmierer, CRISPR-StAR and our own readout — making Q5 the round-1 enzyme at most, not the
enzyme that did most of the amplifying. 🔴 Which enzyme ran round 1 on this sheet is not
stated in what we have. So the enzyme difference between Project Score and our data is
**unresolved**, and the earlier statement that the preprint's external evidence "comes from a
Q5 screen" should be softened until the sheet is dated and its round-1 enzyme read off.

✅ The **1 ng carry-over is not a bottleneck**, which is worth knowing because it looks like
one. At a ~352 bp amplicon, 1 ng is **2.6 × 10⁹ molecules ≈ 29,000 copies per guide** —
about **3,500× more** per guide than the ~8.4 genome copies that entered round 1:

| | copies per guide |
|---|---|
| entering round 1 (5 µg gDNA) | **≈ 8.4** |
| entering round 2 (1 ng product) | ≈ 29,000 |

🟡 So library complexity is set entirely by **round 1's template input**, and round 2 cannot
recover diversity round 1 failed to sample — it only adds cycles. That is the general shape of
every two-step readout here, and it is why the round-1 input figure is the one worth arguing
about. It also means the *cycles* in round 2 contribute amplification bias without contributing
sampling noise, which is the opposite of round 1.

⚠ Minor: Qiagen PCR purification columns are guanidine-based, and the Broad protocol
separately warns that guanidine carry-over distorts A260 quantitation — relevant here because
the next step is measured in **ng**.

🟡 **Nobody reports a GC additive.** No DMSO, no betaine, no GC buffer, no 7-deaza-dGTP in
any of the four. If GC bias is enzyme/buffer-dependent, that is an untested axis across the
entire comparison set.

🟢 **Three of the five readouts use KAPA HiFi; Behan/Project Score uses Q5.** That is the one
real enzyme contrast in the set. KAPA HiFi is the most GC-even of the mainstream
high-fidelity enzymes; Q5 is a different buffer/enzyme system entirely. ✅ It matters for
interpreting our own result: the residual high-GC deficit measured on Schmierer
(≈**−0.031** log2, ~2 %) is what survives *best-in-class* chemistry. And it matters for the
external check — **the preprint's GC-bias evidence from Sanger SCORE comes from a Q5 screen,
while our own data are KAPA HiFi**, so those two numbers are not measuring the same chemistry
and should not be presented as if they were.

🟡 **The Yusa lab does use GC buffer — just not here.** 🟢 Koike-Yusa amplify individual
off-target loci "with Phusion High-Fidelity polymerase **in GC buffer** (Thermo Scientific)".
So the omission in the library readout is a choice, not an oversight, and it shows the
additive axis was available to them.

## 🟡 Also relevant: the Broad GPP protocol itself uses Ex Taq

🟢 The Broad GPP "PCR of sgRNAs for Illumina sequencing" protocol — the one most Brunello
users follow, and the basis of the workflow in `to_debug/` — specifies **Ex Taq DNA
polymerase (Clontech RR001A)**, a **single 28-cycle PCR** (95 °C 30 s / 53 °C 30 s /
72 °C 30 s), **up to 10 µg gDNA or 200 ng plasmid** per 100 µl reaction, ≥4 parallel
reactions, and 8 pooled staggered P5 primers for flowcell diversity.

⚠ Ex Taq is **Taq-based and non-proofreading** — the GC-biased end of the range. So the most
widely used pooled-screen readout in the field is also the most GC-bias-prone of any protocol
in this table, and it is the one real-world workflow here that uses neither KAPA HiFi nor Q5.
At 10 µg per reaction it is also the highest-input, i.e. the *least* first-cycle-stochastic —
the two effects push in opposite directions, which makes it an informative comparison point
rather than simply "the bad one".

## How the gDNA is purified, and what is done to it before the polymerase sees it

Extraction chemistry and any pre-PCR treatment, because both change what the polymerase is
handed: carry-over inhibitors, fragment length, and how accessible the target is inside a
3-Gb genome.

| | lysis / extraction | **pre-PCR treatment** |
|---|---|---|
| **Schmierer 2017** | Qiagen **Blood & Tissue Maxi** kit (silica column) 🟢 | 🔴 none stated |
| **Michlits 2017** | SDS lysis (10 mM Tris pH 8, 1 % SDS, 10 mM EDTA, 100 mM NaCl) + 1 mg/ml ProtK + RNase A → **phenol** → isopropanol 🟢 | ✅ **PacI digest + size-selective SpeedBeads precipitation** → cassette enriched **10³–10⁴×** |
| **CRISPR-StAR 2025** | same SDS/ProtK lysis, **48–72 h** → **phenol-chloroform** → isopropanol 🟢 | ✅ **PacI digest** — 🟢 "facilitated the **accessibility of the DNA** for PCR amplification" |
| **Behan 2019** | **QIAsymphony DSP DNA Midi** (automated) or Qiagen Blood & Cell Culture Maxi 🟢 | 🔴 none stated |
| *Yusa bench protocol* 🟡 | **DNeasy Blood & Tissue** kit, ⚠ **with RNase A** | 🔴 none stated |
| **Broad GPP protocol** | (user-supplied) | 🔴 none; ⚠ warns that **guanidine isothiocyanate** carry-over from gel extraction obscures quantitation (A230) and advises isopropanol precipitation |
| **CRISPR-MIP (ours)** | 50 mM Tris, 50 mM EDTA, 1 % SDS pH 8 + ProtK 100 µg/ml **overnight 55 °C** → RNase A 50 µg/ml → 7.5 M NH₄OAc → isopropanol → EtOH wash → Qubit 🟢 | ✅ **HindIII/BamHI digest**, and 🟢 **benchmarked against undigested** |

### 🟡 Restriction digestion before amplification is a real, deliberate variable

Three of the six cut the gDNA first, and two state why in terms that bear directly on the
GC hypothesis:

- 🟢 **CRISPR-StAR**: PacI digestion "facilitated the **accessibility of the DNA** for PCR
  amplification."
- 🟢 **Michlits**: PacI + size selection "minimized the number of PCR cycles required for
  amplification, **thereby reducing PCR amplification biases**."

⚠ Both are describing the same mechanism our own preprint proposes — that a target embedded
in long, structured, high-molecular-weight genomic DNA is harder to melt and prime than the
same target on a short free fragment, and that this is where sequence composition bites. 🟡 If
GC bias is partly a *template accessibility* effect rather than a pure polymerase property,
then **digestion is as much a candidate knob as enzyme choice**, and it is one the field
already varies without treating it as a variable.

### ✅ Does anyone test whether the digest actually helps?

Asked directly of each paper that includes one. The answer differs sharply:

| | digest | **tested?** | what was measured |
|---|---|---|---|
| **Michlits 2017** | PacI + size selection | ✅ **yes** | 🟢 Fig 2b: qPCR of the **589-bp target amplicon** against a **control amplicon on a 7.7-kb PacI fragment**, across gDNA / digested / fraction 1 (large) / fraction 2 (small). Gives the 10³–10⁴× enrichment figure. Technical triplicate, "repeated 3 times with different samples" |
| **CRISPR-StAR 2025** | PacI | 🔴 **no** | ⚠ **one sentence, no data.** "PacI restriction digestion facilitated the accessibility of the DNA for PCR amplification" — ✅ and `PacI` occurs **exactly once in the entire paper**. No figure, no supplementary panel, no citation supporting it |
| **CRISPR-MIP (ours)** | HindIII/BamHI | ✅ **yes, twice** | 🟢 Fig 1e: **unique UMIs after deduplication, digested vs undigested** at matched input (150,000 cell equivalents). 🟢 Supplementary Fig 1a: **qPCR amplification curves across different input amounts** for digested (pink) vs undigested (turquoise), plus Supp Fig 1b "analysis of fragment amplification based on qPCR" |

🟡 **So the field's position is weaker than the uniformity of practice suggests.** Michlits
measured *enrichment* — which is what their PacI + size-selection step is for, and it works —
but enrichment is not the same claim as "amplification is less biased"; that part is inferred
from needing fewer cycles, not measured. CRISPR-StAR simply asserts the accessibility benefit
and carries the step over from Michlits without testing it.

✅ **Our own Supplementary Fig 1a is the strongest evidence anyone has**, and it is stronger
than we have been giving it credit for: it crosses **digest state × input amount**, which is
exactly the two-factor design the accessibility hypothesis predicts. If digestion works by
making a target in long genomic DNA reachable, the benefit should be **largest at high input**
(where viscosity and incomplete denaturation bite hardest) and should shrink at low input.
🟡 That interaction is already in the figure.

🔴 **What nobody has tested is whether digestion changes GC bias specifically.** All three
measurements are about yield or enrichment, in aggregate, with no per-guide or
composition-resolved readout. That is the gap, and our own digested/undigested pair is the
material to close it with.

✅ **We are the only group that benchmarked the digest functionally.** The preprint runs
"an equivalent of 0.15 × 10⁶ cells (1.0 µg) of **either digested (BamHI/HindIII) or
undigested gDNA**" through one CRISPR-MIP reaction and compares unique-UMI recovery
(Fig 1e). That comparison is the closest thing in this whole table to a controlled test of a
pre-PCR variable — 🟡 but it was framed as a yield/efficiency check, and as far as the
preprint goes it was **not analysed for GC dependence**. Re-reading those two samples with
the GC machinery in [`../README.md`](../README.md) would cost nothing extra and tests the
accessibility mechanism directly.

### ⚠ The RNase detail is not cosmetic — it changes what "µg of DNA" means

🟡 The Yusa bench protocol specifies **DNeasy Blood & Tissue** and explicitly recommends
**RNase A during the prep**. DNeasy does not remove RNA on its own, and residual RNA absorbs
at 260 nm. A prep quantified by A260 without RNase therefore **over-states DNA**, so the true
template going into the reaction is *lower* than the nominal figure — and by an amount that
varies with cell type and harvest state rather than consistently.

⚠ This propagates straight into the only number that matters for the first-cycle argument:
**copies of each guide per reaction**. An input nominally 5 µg but really 3 µg drops ≈8.4
copies per guide to ≈5. So between the 1 µg / 5 µg ambiguity above and quantitation method,
the effective template depth for Project Score is uncertain by something like **5-fold**.
🟡 Anyone comparing GC bias across these datasets should treat per-reaction copy number as a
poorly determined quantity, not a known one — and should say how DNA was quantified (A260 vs
a dye method like Qubit, which does not see RNA).

✅ For contrast, our own protocol uses RNase A (50 µg/ml) and quantifies by **Qubit**, so its
µg figures are the most trustworthy in the table.

🟡 Note also the extraction chemistries split cleanly: **column/kit** (Schmierer, Behan, the
Yusa bench protocol) vs **phenol or phenol-chloroform** (Michlits, StAR, ours). Silica columns
size-bias against very large fragments and can carry chaotropic salts; phenol preps retain
longer DNA. Nobody reports fragment-size QC, so this is an uncontrolled difference across the
set.

## Per-paper detail

### Schmierer 2017 🟢
- **Library**: oligos on CustomArray; a 119-bp oligo carrying the 6-bp RSL and the i7 site is
  annealed to the array oligo and "double-stranded using outer primers". ⚠ Note this is a
  **fill-in, not an exponential PCR** — so the plasmid library has very little PCR history,
  and the GC effect we measure in the plasmid input is dominated by the readout, not by
  library construction.
- Assembly: NEBuilder HiFi DNA Assembly Master Mix (NEB); Endura electrocompetent *E. coli*;
  colonies grown **overnight at 30 °C** (low temperature, reduces recombination).
- **Readout**: 200 µg gDNA across **40 parallel PCR1 reactions** (5 µg each), KAPA HiFi
  HotStart, **14 cycles**; pooled → PCR2 **19 cycles** → PCR3 **14 cycles**.
  🟡 The enzyme is named only for PCR1; PCR2/3 presumably the same but it is not re-stated.

### Michlits 2017 🔴 — the gap
- **Library**: oligos on CustomArray, subsets amplified "with specific flanking primers with
  **ten cycles of PCR**"; enzyme 🔴 not stated. BbsI digest, Golden Gate with T4 ligase.
- **Enrichment before readout**: PacI digest + size selection on SpeedBeads, enriching the
  cassette **10³–10⁴-fold**. 🟢 Stated purpose: "This enrichment minimized the number of PCR
  cycles required for amplification, **thereby reducing PCR amplification biases**." So this
  paper is explicitly engineering *against* the effect we are measuring.
- **Readout**: "Each Sample was PCR-amplified in **200 individual 50-µl PCR reactions**" —
  🔴 **no polymerase, no buffer and no cycle number given**. Only the qPCR enzyme is named
  (GoTaq, Promega), and GoTaq is Taq-based, i.e. at the GC-biased end — but it is used for
  quantification, not for the library.
- 🟡 A step-by-step protocol is cited at Protocol Exchange (ref 51); the enzyme may be there.
  Worth retrieving before using this dataset for an enzyme comparison.

### CRISPR-StAR 2025 🟢
- **Library**: each oligo subpool amplified with **NEBNext High-Fidelity 2× PCR Master Mix**
  (Q5-based), then Golden Gate into the CRISPR-StAR vector or pLenti-UMI.
- **Readout**: **48 reactions of 50 µl, 4 µg DNA each**, KAPA HiFi HotStart ReadyMix (Roche);
  pooled per sample, quantified on a fragment analyzer, pooled by concentration.
- 🟡 The authors were already thinking about amplification bias, though about a different
  axis: "the polymerase chain reaction (PCR) strategy should exhibit minimal bias toward
  amplifying one conformation over the other" (active vs inactive sgRNA).

### Behan 2019 — Project Score / Sanger SCORE 🟢

The dataset the CRISPR-MIP preprint uses for its external GC-bias check
(`score.depmap.sanger.ac.uk`, Release 1, `raw_sgrnas_counts.zip`).

- **Scale**: 324 human cancer cell lines, 30 cancer types.
- **Library**: Human CRISPR Library **v1.0** (Addgene 67989), 18,009 genes / **90,709 sgRNAs**;
  v1.1 adds 1,004 non-targeting plus 5 extra guides against 1,876 genes. 🟢 All analyses use
  the **90,709 sgRNAs common to both**, so a GC analysis should use that intersection too.
- **Screen**: 3.3 × 10⁷ cells transduced at 30 % efficiency = **100× library coverage**;
  **500× coverage** maintained at each passage; ~2.5 × 10⁷ cells harvested.
- **gDNA**: QIAsymphony DSP DNA Midi Kit (Qiagen 937255) or Blood & Cell Culture DNA Maxi Kit
  (13362).
- **Sequencing**: **19-bp single-end** on **HiSeq 2000 v4** with a custom primer. ⚠ 19 bp is
  one base *shorter* than the 20-nt spacer, so the last base is never read — guides are
  identified on 19 of 20 bases, and any guide pair differing only at position 20 is
  indistinguishable. Worth checking against the library before a GC analysis.
- 🟢 Custom sequencing primer (from Tzelepis):
  `TCTTCCGATCTCTTGTGGAAAGGACGAAACACCG` — the same U6 anchor Schmierer uses.

#### 🔴 The amplification chemistry is not retrievable

A three-step citation chain that terminates behind a paywall:

| | says |
|---|---|
| **Behan 2019** | "PCR amplification, Illumina sequencing … were performed **as described previously¹⁰**" |
| **ref 10 = Tzelepis 2016**, *Cell Rep* 17:1193 (PMC5081405) | "Genomic DNA extraction and Illumina sequencing of gRNAs were conducted **as described previously (Koike-Yusa et al., 2014)**" |
| **Koike-Yusa 2014**, *Nat Biotechnol* 32:267 | 🔴 not open access, **no PMC record**; the Springer supplementary that is reachable (`MOESM2`) contains figures only |

✅ Checked: Addgene's protocol for the library (#67989,
`yusa_hcrispr_ko_replication_protocol.docx`) covers **bacterial replication of the plasmid
pool**, not the screen readout — useful for the library row above, useless for the polymerase.

> ⚠ So for the one dataset where we report a GC bias from someone else's data, **we do not
> know the polymerase, the buffer, the cycle count, or the DNA input.** If the reviewers'
> question is answered in terms of enzyme choice, that gap should be stated rather than
> glossed — or closed by writing to the authors.

### CRISPR-MIP (ours) 🟢
The only one of the four whose **capture** step is enzymatic rather than PCR, so it has an
extra enzyme in the chain:

| step | conditions |
|---|---|
| hybridisation | 1–10 µg gDNA, 0.2–0.002 µM probe, **1.5× Ampligase buffer**; 94 °C 5 min, ramp **−0.1 °C/s** to 60 °C, hybridise overnight |
| fill-in + ligation | **4 U Phusion HF** + **5 U Ampligase**, 0.2 mM dNTPs, **60 °C for 1 h** |
| exonuclease | 10 U Exo I + 50 U Exo III, 37 °C; inactivated 80 °C 20 min |
| PCR off the circle | **KAPA HiFi HotStart ReadyMix**; 95 °C 3 min, n × (98 °C 20 s, 60 °C 15 s, 72 °C 30 s) |
| conventional comparison | KAPA HiFi HotStart ReadyMix; PCR1 n × (98/20 s, 60/30 s, 72/30 s), PCR2 n × (98/20 s, 60/20 s, 72/30 s) |
| cycle number | set per screen by qPCR (KAPA HiFi + SYBR Green I) |

⚠ **Two things here are worth a second look given the hypothesis.**

1. 🟡 The gap-fill uses **Phusion**, not KAPA — and Phusion is the more GC-sensitive of the
   two. It runs **once per captured molecule**, so any bias is not compounded across cycles,
   which is the method's central advantage. But a single biased fill-in still sets which
   molecules circularise at all, and a molecule that fails to fill is lost entirely rather
   than merely under-amplified. That makes the fill-in an **all-or-nothing** GC filter, which
   is arguably worse per event than a mild per-cycle bias.
2. 🟡 The fill-in runs at **60 °C**, well below Phusion's 72 °C optimum, for a full hour. Low
   extension temperature slows the enzyme and is generally less forgiving on structured,
   GC-rich template. The temperature is set by the ligation (Ampligase needs 60 °C), so it is
   a real design constraint, not an oversight — but it is exactly the kind of choice the
   hypothesis predicts would matter.

🟢 The preprint already states a mechanism in these terms: "If naked plasmid DNA is more
accessible than gDNA … then the melting temperature can be the rate limiting factor for
amplification. We reason that low GC%, and thus low melting temperature, has an effect in the
first cycles of PCR."

⚠ Note this predicts a **stronger** GC effect on gDNA than on plasmid. ✅ Our re-analysis of
Schmierer found them **equal** once matched on reads-per-molecule (see
[`../README.md`](../README.md)) — but that is Schmierer's chemistry (KAPA HiFi, a 288-bp
amplicon), not ours, so it does not directly test the preprint's claim. Testing it properly
needs our own plasmid-vs-gDNA pair.

## What would actually test the hypothesis

🟡 In increasing order of effort:

1. **Same template, two enzymes.** Split one gDNA prep and run the readout with KAPA HiFi vs
   Phusion (or Q5, or a Taq-based mix). Everything else identical. This is the direct test
   and nothing in the public data substitutes for it.
2. **Same enzyme, ± additive.** KAPA HiFi with and without DMSO or betaine. Cheap, and it
   addresses the axis no published dataset covers.
3. **Fill-in enzyme for CRISPR-MIP.** Whether swapping Phusion for KAPA HiFi (or a
   GC-optimised fill-in) at the capture step changes the GC profile of what circularises.
   This is the one that could change the method rather than just characterise it.
4. **Cycle-number titration**, which also separates per-cycle bias from one-off capture bias.

🔴 None of the four papers varies enzyme, buffer or additive, so **the hypothesis cannot be
tested on public data from this set.** It needs bench work.
