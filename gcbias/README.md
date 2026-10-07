# PCR GC bias, measured against a lineage UMI

Built to answer a CRISPR-MIP reviewer question: how much of the read-count variation
between sgRNAs is PCR GC bias rather than biology? A lineage UMI is the lever, because it
gives a molecule count that is independent of how hard each molecule was amplified.

Start with **plasmid libraries**: no selection, so a GC trend there is amplification rather
than differential growth. Dataset survey in [`datasets/`](datasets/README.md).

> **Data policy.** Nothing here downloads into the repo — see "Data policy" in the
> [root README](../README.md). `download/` holds the scripts, `$CHEM_DATA` (default
> `../_data/`, gitignored) holds the bytes.

## Workflow

```sh
# 0. the published guide libraries (small; fetch once)
python3 download/get_libraries.py --all

# 1. get a slice -- the whole file is 49 GB, and 25M reads answers most questions
python3 download/get_schmierer.py input --max-reads 25000000

# 2. FASTQ -> (guide, UMI) counts. --collapse 1 folds UMI sequencing errors
python3 tools/extract.py $CHEM_DATA/schmierer/input.25000000.fq.gz \
    --sample-index ATCACG --min-guide-reads 100 --collapse 1 \
    --out-prefix $CHEM_DATA/schmierer/input.collapsed

# then re-extract with guide error-correction against the PUBLISHED library.
# GC-rich guides lose more reads to basecalling error, which mimics amplification bias.
python3 tools/extract.py $CHEM_DATA/schmierer/input.25000000.fq.gz \
    --sample-index ATCACG --min-guide-reads 100 --collapse 1 \
    --whitelist $CHEM_DATA/libraries/schmierer.whitelist.tsv --correct-guide \
    --out-prefix $CHEM_DATA/schmierer/input.pub

# 3. LOOK FIRST -- saturation, error, and the GC/abundance confound
python3 tools/qualitative.py $CHEM_DATA/schmierer/input.collapsed.guide.tsv \
    --guide-umi $CHEM_DATA/schmierer/input.collapsed.guide_umi.tsv
Rscript R/plot_qualitative.R $CHEM_DATA/schmierer/input.collapsed.guide.tsv 6

# 4. only then quantify
python3 tools/gcbias.py $CHEM_DATA/schmierer/input.collapsed.guide.tsv \
    --guide-umi $CHEM_DATA/schmierer/input.collapsed.guide_umi.tsv

# model-free cross-check: holds the molecule count fixed, uses no estimator
python3 tools/stratified.py $CHEM_DATA/schmierer/input.collapsed.guide.tsv

# BEFORE comparing samples: equalise reads-per-UMI, not total reads
python3 tools/rarefy.py $CHEM_DATA/schmierer/d4.pub.guide_umi.tsv --target-rpu 1.649

# plasmid vs genomic, paired per guide. MATCH DEPTH FIRST (extract.py --max-reads)
python3 tools/paired.py $CHEM_DATA/schmierer/input.matched.guide.tsv \
    $CHEM_DATA/schmierer/d4.pub.guide.tsv

# control: is guide-read sequencing error GC-dependent? (yes -- correct it first)
python3 tools/error_control.py $CHEM_DATA/schmierer/input.25000000.fq.gz \
    $CHEM_DATA/schmierer/input.collapsed.guide.tsv
Rscript R/plot_gcbias.R $CHEM_DATA/schmierer/input.collapsed.gcbias.tsv

# expected vs actual read count, coloured by GC%  (5 panels + a confound-control panel)
Rscript R/plot_expected_vs_actual.R $CHEM_DATA/schmierer/input.collapsed.gcbias.tsv

# calibration: what the naive statistic does when there is NO GC effect
python3 tools/simulate.py
python3 tools/selftest.py        # 71 checks
```

## The model

In [`lib/umimodel.py`](lib/umimodel.py), with the derivation in its docstring. For one
guide, label space `M`, `N` lineages, reads per molecule ~ NB(mu, r):

```
P(a label gets 0 reads) = exp(-lam * (1 - p0)),   lam = N/M,  p0 = (r/(r+mu))^r
E[distinct UMIs]  D = M * (1 - exp(-lam * (1 - p0)))
E[total reads]    R = N * mu
```

`D` is what dedup gives you, `R` is what the read counter gives you, and the model says how
they relate when nothing GC-dependent is happening. Fit **one global** `(mu, r)` across all
guides, invert each guide's `D` to a molecule estimate, predict its reads, and regress
`log2(observed / predicted)` on GC%. The prediction uses only that guide's molecule count,
so the residual is the guide amplifying differently from the library average.

Two corrections matter and both get worse with abundance, which is why the naive
`reads / distinct UMI` ratio manufactures a GC trend:

1. **label collisions** — two lineages drawing the same UMI count once;
2. **unseen molecules** — a molecule that got zero reads is invisible.

### ⚠ Why not just deduplicate

`tools/simulate.py` simulates a library with **no GC effect on amplification at all**, but
with abundance correlated with GC (true of real libraries). The truth is a zero slope:

| lineages/guide | median occupancy | naive slope | model slope |
|---|---|---|---|
| 200 | 0.049 | −0.003 ± 0.018 | −0.046 ± 0.018 |
| 800 | 0.178 | **+0.165** ± 0.010 | +0.001 ± 0.010 |
| 2,000 | 0.383 | **+0.387** ± 0.007 | −0.001 ± 0.008 |
| 4,000 | 0.622 | **+0.703** ± 0.006 | −0.004 ± 0.007 |
| 12,000 | 0.947 | **+1.440** ± 0.004 | −0.022 ± 0.012 |

And injecting a *real* −1.0 slope, the naive statistic attenuates it to −0.29 at 62 %
occupancy while the model recovers −1.00. So the naive ratio is wrong in **both**
directions depending on regime — it invents signal where there is none and shrinks signal
where there is. Schmierer's ~3,355 lineages per guide in a 4,096-label space sit squarely in
the worst part of that table.

### The non-parametric route, which is what actually worked

The compound-NB fit is **rejected by the Schmierer plasmid input**: the read-per-UMI
spectrum is far more heavy-tailed than any NB that also reproduces the observed number of
distinct UMIs (at `m1 ≈ 2.2` the model wants occupancy ≈ 0.30; the data show 0.101).
`fit_moments` reports this rather than fitting anyway. PCR jackpotting is heavier-tailed
than a negative binomial, which is itself worth knowing.

So `umimodel` also provides an estimator with no amplification model at all:

```
chao1_labels(distinct, f1, f2)   # singletons/doubletons -> labels incl. unseen
labels_to_molecules(labels, M)   # -M * ln(1 - labels/M), undoes collisions
nonparametric_molecules(...)     # the two composed
good_turing_coverage(reads, f1)  # 1 - f1/R; low means too shallow to claim anything
```

✅ Validation of the collision step: 2,290 occupied labels in a 4,096 space inverts to
**3,354** molecules, against the paper's independently reported **3,355** lineages per
guide.

## What the data say so far

**Schmierer `RKO_NTU_INPUT`, first 25M reads (1 % of the run), UMIs collapsed at edit
distance 1.** 23,291 guides pass a ≥100-read threshold — against a library of 23,250, so
the threshold cleanly separates the real library from 1.34 M error guides.

🟢 Qualitative, and this is the part to trust:

| | |
|---|---|
| UMI occupancy | median **0.101** (max 0.153) — comfortably unsaturated at this depth |
| error inflation | **27 % of UMIs folded** by edit-distance-1 collapse |
| error threshold | 138,250 reads/guide; we are at ~1,070, so well clear |
| Good-Turing coverage | **0.71** — ⚠ shallow; only 71 % of reads are on UMIs seen ≥2× |
| GC vs abundance | Spearman **+0.07** — weak confound, but present |

✅ Quantitative, three estimators, `log2` per unit GC fraction:

| statistic | slope | per 10 % GC | Spearman |
|---|---|---|---|
| naive `log2(reads/distinct)` | −0.065 ± 0.010 | −0.0065 | −0.035 |
| NB model `log2(obs/pred)` | −0.078 ± 0.010 | −0.0078 | −0.047 |
| **Chao1 `log2(reads/molecule)`** | **−0.107 ± 0.018** | −0.0107 | −0.031 |

All three agree in sign: **GC-rich guides yield slightly fewer reads per molecule.** The
effect is small — about 0.7 % fewer reads per molecule per 10 % GC.

⚠ But it is **not linear**, and that is the substantive finding. A GAM in
`R/plot_gcbias.R` gives `edf 6.5/9, F = 27.3, p ≈ 0`, and the binned medians show why:

| GC bin | median bias (log2) |
|---|---|
| 0.45–0.50 | +0.014 |
| 0.50–0.55 | +0.014 |
| 0.55–0.60 | +0.008 |
| 0.60–0.65 | +0.003 |
| 0.65–0.70 | −0.010 |
| 0.70–0.75 | **−0.035** |

Flat through the middle, with the deviation concentrated above ~0.65 GC. 🟡 This matches the
prior expectation that GC bias is not a big global effect but does hit the extremes, and it
argues for reporting the end bins rather than a single slope — a linear slope on this data
reports `p = 5e-15` for a shape that is mostly flat.

⚠ Also worth noting: a large part of the *apparent* trend before error-collapse was error,
not amplification. Uncollapsed, the model slope was −0.110 with Spearman −0.22; after
collapsing it falls to −0.078 with Spearman −0.047, and the quintile pattern stops being
monotone. **Collapse UMIs before quoting anything.**

## Expected vs actual read count

`R/plot_expected_vs_actual.R`. Expected = Chao1 molecule estimate × one global
reads-per-molecule, so the totals match by construction and the **identity line is the
reference**: a guide off the line amplifies differently from the library average.

⚠ Read the panels in order, because the raw scatter cannot show the effect:

| panel | what it is |
|---|---|
| 1 | expected vs actual, log-log, points coloured by GC — the honest per-guide picture |
| 2 | the same split into GC bands with a fit per band |
| 3 | actual/expected vs expected — checks the prediction is unbiased in abundance |
| 4 | **binned medians per GC band** — the panel to actually read |
| 5 | offset from identity per GC band, median ± 95 % CI |
| 6 | GC effect *within* abundance quartiles (separate PDF) — the confound control |

✅ Per-guide noise is sd **0.298** log2 while the whole GC effect spans **0.088** log2, a
**3.4 : 1** noise-to-signal ratio, so panel 1 looks like an undifferentiated cloud however it
is coloured. That is a property of the data, not the plot. Binning is what makes GC legible,
and panel 4 shows the GC > 0.7 band sitting below identity across the whole abundance range.

✅ Deviation from expectation by GC band:

| GC band | n | log2(actual/expected) | 95 % CI |
|---|---|---|---|
| 0.15–0.45 | 4,954 | **−0.034** | −0.042, −0.027 |
| 0.45–0.55 | 6,868 | +0.001 | −0.005, +0.008 |
| 0.55–0.60 | 3,523 | −0.005 | −0.015, +0.004 |
| 0.60–0.70 | 5,856 | −0.025 | −0.033, −0.017 |
| **0.70–1.00** | 2,090 | **−0.087** | −0.101, −0.072 |

🟡 So the effect is **not a slope**: the middle sits on expectation and the high-GC end falls
clearly below it. A single linear slope misrepresents this, which is why the earlier
regression reported a confident `p` for a mostly-flat shape.

⚠ The apparent **low-GC** deficit in the table above does **not** survive the model-free test
below — it is an artefact of the molecule-estimation step, not a real effect. Only the
high-GC deficit is robust.

⚠ Confound control. GC correlates with depth (Spearman +0.148) and the ratio is not quite
flat in depth (−0.124), so the bands could partly be an abundance effect. Stratifying into
abundance quartiles, the pattern survives in all four — high-GC minus mid-GC is −0.181,
−0.042, −0.071, −0.032 from lowest to highest quartile — but the **magnitude is unstable**,
and it is largest in the lowest-abundance quartile, exactly where Chao1 is least reliable.
Treat the sign as established and the size as provisional until the deeper run.

🟡 Also worth stating plainly: the UMI-based expectation explains only part of the per-guide
read count — Pearson 0.52 on log10, Spearman 0.45. At this depth it is a weak per-guide
predictor and only useful aggregated.

## ⚠ Non-uniform guide abundance: what the model does and does not handle

Guide representation in a plasmid library is uneven, and every estimator here is a function
of abundance, so this needs answering explicitly rather than assuming.

| | handled? | why |
|---|---|---|
| **Guides differ in abundance from each other** | ✅ **yes, by construction** | `N_g` is a free per-guide parameter, estimated from that guide's own UMI data. The model never assumes equal representation — comparing reads to *that guide's* molecule count is the entire design. |
| **Lineages differ in size within a guide** | 🟡 partly | The global `r` absorbs an average overdispersion but assumes it is the same for every guide. Shape statistics do drift with GC (dispersion +0.037, singleton fraction +0.017 per unit GC), plausibly as a consequence of the mean shift rather than independently. |
| **Production bias vs PCR bias** | 🔴 **no — and this is the important one** | See below. |
| **Non-uniform UMI (RSL) pool** | 🔴 no | `labels_to_molecules` assumes a uniform UMI draw; the observed RSL distribution is C-skewed. Open item 5. |

### The model-free check — `tools/stratified.py`

Because the first three rows all involve abundance, the claim needs evidence that uses no
model at all. Guides are stratified into narrow bins of **observed distinct-UMI count** (a
direct molecule proxy — no Chao1, no NB, no collision correction), and within each stratum
`log2(reads)` is regressed on GC.

✅ All 10 strata give a negative slope, mean **−0.225** log2 per unit GC. By band, relative
to mid-GC guides *at the same molecule count*:

| stratum (distinct UMIs) | GC<0.45 | 0.45–0.55 | 0.55–0.60 | 0.60–0.70 | GC>0.70 |
|---|---|---|---|---|---|
| 35–369 | −0.010 | 0 | −0.010 | −0.021 | **−0.100** |
| 369–402 | +0.004 | 0 | −0.005 | −0.025 | −0.051 |
| 402–429 | +0.006 | 0 | −0.007 | −0.024 | −0.049 |
| 429–460 | −0.003 | 0 | −0.004 | −0.022 | −0.048 |
| 460–626 | −0.014 | 0 | −0.010 | −0.014 | −0.029 |
| **mean** | **−0.003** | 0 | −0.007 | −0.021 | **−0.055** |

🟢 **The high-GC deficit is real**: negative in 5/5 strata, ~−0.055 log2 (≈3.7 % fewer reads
per molecule). 🔴 **The low-GC deficit is not**: −0.003 on average, negative in only 3/5
strata. So the U-shape seen through the Chao1 statistic was half artefact — the low-GC arm
came from the estimator, not the chemistry.

✅ This test is also **conservative**: if GC-rich guides amplify less they lose reads *and*
lose detected UMIs, so at fixed distinct-UMI count they carry slightly more true molecules
and should show *more* reads. The measured deficit is a lower bound.

### 🔴 What "amplification" means here, and why it is not only PCR

The UMI labels a **cloning event**, not a plasmid molecule. Between that event and the
sequencer sit bacterial transformation, colony growth, plasmid copy number, miniprep, and
only then PCR. All of those are GC-sensitive. So the −0.055 log2 is **reads per cloning
event**, which bundles plasmid production together with PCR.

⚠ For a reviewer asking specifically about **PCR** GC bias, the plasmid input alone cannot
answer it. Separating them needs a contrast where production is held constant and
amplification varies:

1. the same plasmid pool amplified at **different cycle numbers** — production bias is common
   to both, PCR bias scales with cycles. Schmierer provides no cycle variation;
2. plasmid input vs the **day-4 genomic** sample — same library, different template and
   primers, and day 4 is early enough that selection is minimal;
3. Michlits' PacI enrichment exists precisely to cut cycle number 🟢 — the right design, but
   its UMI was never deposited.

Option 2 has now been **run** — see the next section.

## ✅ Plasmid input vs day-4 genomic — production bias is ruled out

The contrast that option 2 above describes, run on 25M reads from each of
`RKO_NTU_INPUT` and `RKO_NTU_Replicate1_D4`, identical settings, guide errors corrected,
UMIs collapsed, model-free stratified statistic (`tools/stratified.py`):

| | GC<0.45 | 0.45–0.55 | 0.55–0.60 | 0.60–0.70 | **GC>0.70** |
|---|---|---|---|---|---|
| **input** (plasmid) | −0.006 | 0 | −0.005 | −0.017 | **−0.049** |
| **day 4** (genomic) | +0.009 | 0 | −0.002 | −0.018 | **−0.054** |

🟢 The two are **the same within noise** (−0.049 vs −0.054 log2, ≈3.4 % vs 3.7 %), negative
in 5/5 strata in both, and the low-GC end is flat in both.

**Why that is informative.** In the plasmid the UMI labels a cloning event, and the reads
from it carry whatever plasmid copy-number and colony-size variation that clone accumulated
in bacteria. By day 4 each lineage is **one integrated provirus in a cell** — plasmid copy
number does not propagate. If the GC deficit came from plasmid production it should shrink
or vanish at day 4. It does not move.

🟡 So the deficit lives in a step **common to both samples**: the readout PCR (same primers,
same amplicon) and/or cluster amplification on the flowcell. ✅ And once the two are measured
at the same per-molecule operating point they are the same *size* as well — see below. 🔴 Those two **cannot** be
separated here — both are shared — so "PCR GC bias" remains the best-supported label but is
not proven against bridge amplification. Cycle-number variation is still the experiment that
would settle it.

⚠ Caveats: day 4 is early but not selection-free; both samples ran on the same HiSeq 4000
chemistry; the i5 filter discards 42 % of the day-4 reads (that file carries a second index,
`ATCCCG`, whose guide counts correlate with the dominant one at only Spearman +0.13, so it
is probably a different sample rather than miscalling); and only 1-mismatch substitution
errors in the guide are corrected.

### ✅ Use the published library, not one derived from the data

Both papers publish their library, and `download/get_libraries.py` fetches them:

| | source | guides |
|---|---|---|
| Schmierer | Dataset EV1 (`MSB-13-945-s002.csv`), via the EuropePMC `supplementaryFiles` endpoint for PMC5658704 | **23,279** unique (23,332 rows; 101 non-targeting controls) |
| Michlits | Supplementary Table 2 (`MOESM3`), via Springer static-content | **26,486** unique (28,565 rows) |

🟢 Schmierer's file is already in the `lib.csv` format the authors' own RSLC pipeline wants
(`GuideID,GuideSequence,TargetGene`, no header), and its 23,279 matches the paper's text
exactly. The script also pulls Michlits' primer table (`MOESM4`) and its **sample ↔ 6-bp
index map** (`MOESM5`), which names the original BAM file per sample — the only thing that
would make those BAMs usable.

⚠ **Deriving a whitelist from the data is circular.** Keeping every 20-mer above a read
threshold works well — a ≥100-read cut recovered 23,272 of the 23,279 — but it fails in a
GC-dependent direction: the 7 real guides it *missed* had median GC **0.70**, while the 19
false ones it *kept* had median GC **0.35**. Low-abundance guides are exactly the ones the
threshold drops, and abundance is what is being measured.

✅ In practice it changed nothing here — the high-GC deficit is **−0.0491** with the published
list against **−0.0493** with the derived one — but the result now rests on the real library
rather than on a threshold tuned against the data.

### ✅ Is the GC effect stronger on genomic DNA than on plasmid? **No** — corrected

An earlier version of this file said yes, ~1.5–1.7×. That was wrong: it compared samples
matched on **reads**, and the statistic is not comparable at matched reads. Corrected below.

#### The flaw: log(R/D) is a compressed function of the thing we want

The statistic is `log2(reads)` at fixed observed distinct UMIs, i.e. essentially `log(R/D)`.
But `R/D` is not the per-molecule amplification rate `mu` — molecules that got zero reads are
invisible, so what is observed is the **zero-truncated** mean `T(mu) = mu / (1 - p0(mu))`.
`T` is monotone in `mu` but compressed, and the compression depends on the operating point:

```
dlogT/dlogmu  ->  0   as mu -> 0     shallow: a real effect is squashed toward nothing
dlogT/dlogmu  ->  1   as mu -> inf   deep:    the effect appears in full
```

✅ Computed for these samples (Poisson; an overdispersed fit gives the same picture):

| | reads/UMI | implied `mu` | compression | measured | implied true |
|---|---|---|---|---|---|
| plasmid | 1.65 | 1.10 | **0.45** | −0.031 | −0.069 |
| day 4 rep 1 | 3.11 | 2.95 | **0.84** | −0.054 | −0.065 |
| day 4 rep 2 | 2.58 | 2.33 | 0.75 | −0.041 | −0.055 |
| day 28 | 4.53 | 4.48 | 0.95 | −0.006 | −0.006 |

⚠ Matching total read depth does **not** equalise reads-per-UMI, because reads-per-UMI also
depends on how many molecules a sample has — and lineage complexity collapses over a screen
(337 UMIs/guide in the plasmid, 177/214 at day 4, 122 by day 28). So the deeper-per-molecule
samples reported a bigger effect for the same chemistry.

#### The fix: match on molecules, not reads — `tools/rarefy.py`

Subsampling reads is exactly **binomial thinning** of each (guide, UMI) read count, so it can
be done on the table without re-reading FASTQ: draw `Binomial(reads, p)` per molecule and drop
molecules that reach zero, which correctly reproduces losing them to the detection limit.
`rarefy.py` bisects `p` until the median reads-per-UMI hits a target. Thin every sample to the
shallowest one's operating point, then compare.

#### ✅ The result, with every sample at reads/UMI = 1.649

| sample | **GC>0.70** |
|---|---|
| plasmid input | **−0.038** |
| genomic day 4, rep 1 | **−0.043** |
| genomic day 4, rep 2 | **−0.026** |
| genomic day 28 | −0.010 |

🟢 **Plasmid −0.038 vs day-4 mean −0.035** (difference +0.003, against a replicate spread of
~0.009). The difference is gone.

⚠ These figures are the **aggregate-then-log** values. An earlier version of this file used
the median of per-guide logs, which gave −0.031 / −0.037 / −0.024 / −0.008. ✅ Audited: at
these count scales (≈550 reads and 120–340 UMIs per guide) the Poisson transform bias is
−0.001 to −0.006, every figure moved by ≤0.007, and **every comparison was preserved** —
including this one. The rarefied day-28 sample is the only one that got close to trouble, at
60 reads per guide (bias −0.012). See "Counting safely" below. The GC effect is the same
on plasmid and on genomic DNA, which is what the earlier "it lives in a step common to both
samples" conclusion already implied — the magnitudes simply had to be measured at the same
operating point to see it.

🔴 Day 28 still reads lower, but it is the least trustworthy row: it was thinned hardest
(p = 0.105) and lost 2,200 guides to the ≥20-UMI floor. Do not read a time trend into it.

> **Rule for any future comparison:** equalise reads-per-UMI with `rarefy.py` before comparing
> samples. Equalising total reads is not sufficient and will inflate whichever sample has
> fewer molecules.

### 🟡 Where a real likelihood would earn its place — and why Schmierer is the case

Everything above is moments, closed forms and bisection, deliberately. But Schmierer's 6-nt
RSL saturates (82 % of the label space occupied at full depth), and that is exactly where the
cheap estimators give out:

> **The key point.** When `D → M`, the distinct-label count stops carrying information about
> the molecule count. But the **per-label read histogram still does** — a label carrying many
> reads is more likely to be multiply occupied. `labels_to_molecules` throws that away; a
> likelihood over the full compound distribution keeps it.

The generative model is a compound Poisson: each guide has `N` molecules thrown into `M`
labels, so a label carries `Poisson(N/M)` molecules, each contributing an amplification draw.
✅ The likelihood is cheap to evaluate exactly by **Panjer recursion** —

```
P(S=0) = exp(-lam')                       lam' = lam*(1-p0), g = amplification | X>=1
P(S=n) = (lam'/n) * sum_k k*g_k*P(S-n-k)
```

— so this needs `scipy.optimize`, not Stan. Stan would only pay for a hierarchical version
pooling guides or samples.

✅ **Validated on simulation** (M = 4096, fitting `D` *and* the full read histogram):

| truth | occupancy | fitted family | N̂ error | AIC |
|---|---|---|---|---|
| NB, N=2,000 | 0.35 | NB | **+1.6 %** | 11,747 |
| NB, N=9,000 | **0.85** | NB | **+3.1 %** | 24,004 |
| NB, N=9,000 | 0.85 | Poisson | −10.0 % | 24,215 |
| **Poisson-lognormal**, N=2,000 | 0.29 | **NB** | 🔴 **+531 %** | 11,016 |
| Poisson-lognormal, N=2,000 | 0.29 | Poisson-lognormal | **−1.6 %** | **10,982** |

🟢 **At 85 % occupancy the full likelihood still recovers N to 3 %**, which the `D`-based
inversion cannot do at all. That is the argument for doing this on Schmierer specifically.

⚠ **But the family matters more than the saturation does.** Fitting a negative binomial to
Poisson-lognormal truth gave N̂ **6× too high**. PCR is multiplicative, so a lognormal
amplification factor is the more physically natural choice — and recall the moment fit was
*rejected* by the Schmierer data precisely because the read spectrum was heavier-tailed than
any NB that also explained the occupancy. ✅ AIC separates the families cleanly, so this is
testable rather than a matter of taste.

#### ✅ Run on the real Schmierer data — Poisson-lognormal wins, and NB is pathological

40 guides from the full-depth input, mean occupancy 0.097, 399 distinct labels per guide.
Shared amplification parameters across guides, per-guide λ profiled out:

| family | fitted params | AIC | mean N̂ per guide |
|---|---|---|---|
| Poisson | μ = 1.55 | 160,802 | 573 |
| negative binomial | μ = 0.014, r = 0.005 | 150,584 | 🔴 **65,418** |
| **Poisson-lognormal** | μ = 0.441, σ = **1.79** | **149,904** | **2,160** |

🟢 **Poisson-lognormal wins by ~680 AIC units.** And the NB fit does not merely lose — it goes
**degenerate**, driving μ and r to the corner of the space and returning 65,418 molecules per
guide, which exceeds anything the library can contain. That is the same pathology the
simulation predicted (+531 % when NB is fitted to lognormal truth), seen on real data.

✅ The winning fit is also *sane*: 2,160 molecules per guide against the paper's independently
reported ~3,355 lineages per guide, from a sample where only 399 labels are observed. And
σ = 1.79 is a very heavy tail — strong PCR jackpotting, exactly what a multiplicative process
gives and what no NB can represent while also explaining the occupancy.

🟡 ⚠ This run was at **0.097 occupancy**, so it tests *family selection*, not saturation
handling. The saturation claim rests on the simulation above (3 % recovery at 0.85). Running
it at full library depth is the next step.

🟡 Remaining families worth adding: **zero-inflated NB** (decouples "molecule never amplified"
from dispersion — the identifiability problem the moment fit hit), and an **effective-M** or
Dirichlet concentration for **non-uniform label usage** — the observed RSL pool is C-skewed,
and non-uniform labels saturate faster than uniform, so a fitted `M_eff < 4096` is likely and
would otherwise be absorbed into N̂.

### ⚠ Counting safely: aggregate, then take the log

A second rule, learned the hard way on the CRISPR-MIP data. Taking the **median of per-guide
logs** is biased, because for Poisson `E[log2 X] = log2(lam) − 1/(2·lam·ln2)`:

| mean count | bias in log2 |
|---|---|
| 3 | −0.240 |
| 5 | −0.144 |
| 10 | −0.072 |
| 55 | −0.013 |
| 550 | −0.001 |

The bias only matters when it **differs between the groups being compared**, i.e. when their
counts differ. That makes it harmless at Schmierer's depth and serious at CRISPR-MIP's, where
UMIs/guide is ~5 and the statistic carried ~10× the bias of the read statistic.

> **Rule:** sum the counts within a group *first*, then take one log. Per-guide ratios
> (`reads/UMI`) are safe either way, because the bias largely cancels between numerator and
> denominator — which is why that statistic is the one to quote.

✅ Everything in this file has been re-run both ways. The only conclusion that changed was on
the CRISPR-MIP data (see [`datasets/crispr-mip.md`](datasets/crispr-mip.md)); all Schmierer
results and the probe-concentration null are unchanged.

### ⚠ Guide-read sequencing error is itself strongly GC-dependent

Found while checking the above, and worth its own note — `tools/error_control.py`:

| GC band | median fraction of reads 1 mismatch off |
|---|---|
| <0.45 | 0.047 |
| 0.45–0.55 | 0.053 |
| 0.55–0.60 | 0.062 |
| 0.60–0.70 | 0.070 |
| >0.70 | **0.083** |

✅ Spearman **+0.59** against GC — far stronger than any amplification effect here. Reads
whose guide is miscalled do not match the library and are silently discarded, so GC-rich
guides lose more reads, which mimics reduced amplification.

🟡 Its actual contribution is smaller than that table suggests: correcting the errors
(`extract.py --correct-guide`, recovering 6.6 % of reads) moves the high-GC deficit only from
**−0.055 to −0.049**, about **11 %** of the effect. The reason is that an error removes a
read *and* often the only read of its UMI, so both sides of "reads per molecule" shrink
together. The raw error-rate table is an upper bound on the confound, not an estimate of it.

> **Practical rule:** error-correct guides before any GC analysis. Not doing so inflates the
> apparent GC effect, and the inflation is itself GC-dependent, so it cannot be absorbed into
> a global normalisation.

## Open items

1. 🔴 **Depth.** Good-Turing coverage of 0.71 is too low for a confident molecule count —
   Chao1 puts the median at 880 molecules/guide against the paper's 3,355, as expected of a
   lower bound on shallow data. Re-run at 250M–500M reads
   (`--max-reads 250000000`, ~6 GB) and check that the GC shape is stable. Occupancy at full
   detection would be ~0.56, which the simulations handle to ±2 %, so there is a workable
   regime; we are simply below it.
2. 🔴 **A heavier-tailed amplification model.** NB is rejected. A Poisson-lognormal or an
   explicit branching-process amplification distribution would restore a parametric
   estimate; R (`MASS`, `mgcv`) is the natural place to fit candidates against the
   `*.umi_spectrum.tsv` output.
3. ✅ **Done** — the day-4 contrast is reported above. What remains is separating readout PCR
   from flowcell bridge amplification, which needs cycle-number variation that no public
   CRISPR-UMI dataset provides.
4. 🟡 **Michlits** would be the better dataset (10-nt UMI, 0.3 % occupancy) but its UMI was
   never deposited; see `get_michlits.py --show-bams`.
5. 🟡 **Non-uniform UMI pool.** `labels_to_molecules` assumes uniform UMI draw, and the
   observed RSL distribution is C-skewed. Non-uniformity reduces the effective label space
   and biases molecule counts upward; worth quantifying from the input sample directly.

## Layout

```
gcbias/
├── download/
│   ├── get_libraries.py   published guide libraries + Michlits' primer and index tables
│   ├── get_schmierer.py   ENA PRJEB18436 -- author-submitted FASTQ (RSL in the header)
│   └── get_michlits.py    SRA PRJNA383356 -- guide-level only, UMI was not deposited
├── tools/
│   ├── extract.py         FASTQ -> (guide, UMI) counts; --collapse folds UMI errors
│   ├── qualitative.py     saturation, error, confound, extremes -- run this first
│   ├── gcbias.py          global fit, per-guide molecule estimate, GC regression
│   ├── stratified.py      model-free GC test at fixed molecule count (no estimator)
│   ├── paired.py          per-guide plasmid-vs-genomic comparison, depth/selection aware
│   ├── rarefy.py          thin to a common reads-per-UMI -- REQUIRED before comparing
│   ├── error_control.py   is guide-read sequencing error GC-dependent? (it is)
│   ├── simulate.py        calibration: what the naive statistic does with no true effect
│   └── selftest.py        71 checks, all simulation-based
└── R/
    ├── plot_qualitative.R occupancy, confound, UMI spectrum
    └── plot_gcbias.R      model vs naive, binned medians, GAM non-linearity test
```
