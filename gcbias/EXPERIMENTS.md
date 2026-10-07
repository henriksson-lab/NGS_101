# Proposed experiments — GC bias in CRISPR screen readouts

Written 2026-10-04. Candidates for settling what the existing data cannot, with the
practical constraints that ruled out the obvious designs recorded alongside each.

Current state of the evidence is in [`README.md`](README.md) and
[`datasets/crispr-mip.md`](datasets/crispr-mip.md). In short: a GC-dependent
amplification bias is established (~3 % fewer reads per molecule at GC > 0.70, consistent
across two libraries, two platforms and four estimators), and UMI deduplication corrects a
real part of it (−0.021 log2 at GC > 0.70 after error-collapse, of which ~76 % survives
the collapse). What is **not** established is that the molecule count is *correct* — there
is no ground truth anywhere in the current data.

## ✅ Already possible on deposited data — do these first

Zero bench cost; both are listed here only so they are not forgotten.

- **Remaining matched conditions.** `normoxia_14d` and `hypoxia_3d` in E-MTAB-14379 carry
  both PCR and padlock arms, as `normoxia_3d` does. Three independent replications of the
  dedup effect rather than one.
- **Full depth.** The comparisons so far use 4M read pairs subsampled per sample. Full
  depth would tighten the replicate scatter, which is currently what limits every
  conclusion, and would let the collapse threshold be tuned rather than inherited.

## 🔴 Ruled out, and why — record these so they are not re-proposed

- **A no-amplification control.** The natural test of an amplification bias is to remove
  the amplification. ⚠ Not possible: without amplification there is nothing to sequence.
  Only *differences* in cycle number are achievable, never a zero point.
- **Naked oligo spike-ins.** Ordering a defined guide mixture at known ratios would give
  ground truth cheaply. ⚠ But synthetic oligos behave like short, fully accessible,
  digested DNA — not like a provirus embedded in chromatin-derived genomic DNA, which is
  the template whose behaviour is in question. A spike-in in the wrong physical state
  would answer a different question convincingly, which is worse than not answering it.

## B. Split one gDNA prep across both protocols at several inputs

**Question.** Does the gap between the read-based and molecule-based estimate scale with
template depth, as first-cycle stochasticity predicts?

**Why it is not already answered.** In E-MTAB-14379 `ug_dna` never varies *within* a
comparable group — the 0.2 / 1 / 7 / 10.1 µg spread tracks what each sample *is* (plasmid
library vs t0 gDNA vs screen gDNA), not a designed titration. The one input series that
does exist (Supplementary Fig 1a) is qPCR only and was never sequenced.

**Design.** One gDNA prep, split, then **1 / 3 / 10 µg** per reaction through **both** the
PCR and the CRISPR-MIP arms. Everything else held constant; three replicates.

**Prediction.** If the bias arises in the first cycles while template is sparse, the
read-vs-molecule gap should **shrink as input rises** (≈7 → ≈65 copies of each guide per
reaction over that range). A flat result would point instead at a per-cycle effect that
input cannot relieve.

**Cost.** 18 libraries, no new reagents, one gDNA prep. The cheapest of the bench options.

**Analysis.** `rarefy.py` to a common reads-per-UMI *before* comparing — input amount
changes molecules per guide directly, which changes reads-per-UMI, which changes the
compression on the statistic. This is the single comparison where that correction matters
most.

## C. Ground truth in the right physical state

**Question.** Is the molecule-based estimate actually *closer to truth*, as opposed to
merely *different* from the read-based one?

**Why nothing else answers it.** Every comparison run so far is one estimator against
another. The argument that molecules are less biased is mechanistic — `reads = molecules ×
amplification`, and molecules skip the amplification term — but it has never been checked
against a known answer. 🔴 **No amount of further internal comparison can license the word
"correct".**

**Design.** A cell line carrying a **known, integrated** guide set at known copy number,
spiked into screen gDNA at known ratios, with the spiked guides spanning the GC range.
Both arms read it out; whichever estimator recovers the known ratios is the better one.

**Why integrated rather than synthetic.** The template must be in the same physical state
as the real one — chromatinised, sheared or undigested as the protocol dictates, at the
same molecular crowding. That is precisely what a naked oligo spike-in cannot provide, and
why the cheap version of this experiment does not substitute.

**Cost.** By far the most work: a line to build or obtain, copy number to verify, ratios
to control. ⚠ But it is the only design in this file that licenses a claim of correctness
rather than of difference — so it is the one to do if the claim needs to survive a
determined referee.

**Lower-cost variant worth considering first.** If a suitable line already exists in-house
(e.g. the one-sgRNA `lentiCRISPRv2-GFPg1` model already used for the digest trial), mixing
two such lines at known cell ratios gives a two-point ground truth for a single GC value —
weaker than a spanning set, but a real check, and probably available now.

## A. Cycle-number titration — the strongest test that needs no ground truth

Kept here despite being listed as ruled-out above, because only the *zero* point is
impossible, not the design.

**Design.** Same gDNA, same library, **12 / 18 / 24 cycles** through the PCR arm, with the
MIP arm as the invariant reference.

**Prediction.** If reads carry an amplification bias that molecules do not, read-based GC
bias should grow roughly **linearly in cycle number** while the molecule-based estimate
stays flat. That is a sharp, directional prediction requiring no knowledge of the truth.

**Cost.** Two extra libraries beyond what a normal run produces.

🟡 Note the cycle count is currently chosen per screen by qPCR and is **not recorded in the
sample metadata**, so even retrospective use of existing runs is impossible. If nothing
else from this file is adopted, recording the cycle number per sample costs nothing and
makes this analysis available later.
